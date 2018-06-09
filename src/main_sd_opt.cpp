#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <algorithm>

#include "Eigen/Dense"
#include "json.hpp"
#include "polynomial.h"
#include "cost.h"
#include "utils.h"
#include "spline.h"
#include "frenet.h"
#include "car.h"
#include "waypoint.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

template<typename T>
T norm(const T dx, const T dy) {
    return sqrt(dx * dx + dy * dy);
}

template<typename T>
bool validStep(const T dx, const T dy) {
    return norm(dx, dy) < mphToMps(50) * 0.02;
}

int main() {

    /* The default type used for all data */
    using T = double;

    /* Load the map */
    const string map_file_ = "../data/highway_map.csv";
    vector<T> map_x_, map_y_, map_s_, map_dx_, map_dy_;
    std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);
    string line;
    while (getline(in_map_, line)) {
        std::istringstream iss(line);
        map_x_.push_back(next<T>(iss));
        map_y_.push_back(next<T>(iss));
        map_s_.push_back(next<T>(iss));
        map_dx_.push_back(next<T>(iss));
        map_dy_.push_back(next<T>(iss));
    }

    /* Refine the map using splines */
    vector<T> map_x_refined, map_y_refined, map_s_refined;
    Refiner<T> refiner(map_x_, map_y_, map_s_);
    std::tie(map_x_refined, map_y_refined, map_s_refined) = refiner.linSpaced(map_x_.size() * 20);

//    /* Now, with this refined map, we need to generated the Frenet coordinates */
//    vector<T> map_s_refined(map_x_refined.size());
//    for (int i = 0; i < map_x_refined.size(); ++i) {
//        const auto sd = getFrenet(map_x_refined[i], map_y_refined[i], static_cast<T>(0), map_x_refined, map_y_refined);
//        map_s_refined[i] = sd[0];
//    }

    /* Print the original and refined maps */
    const bool printing = true;
    if (printing) {
        cout << "Original maps" << endl << "[" << endl;
        for (size_t i = 0; i < map_x_.size(); ++i)
            cout << map_x_[i] << "," << map_y_[i] << "," << map_s_[i] << endl;
        cout << "]" << endl << endl;
        cout << "Refined maps" << endl << "[" << endl;
        for (size_t i = 0; i < map_x_refined.size(); ++i)
            cout << map_x_refined[i] << "," << map_y_refined[i] << "," << map_s_refined[i] << endl;
        cout << "]" << endl << endl;
    }

//    const auto map_x = map_x_;
//    const auto map_y = map_y_;
//    const auto map_s = map_s_;
    const auto map_x = map_x_refined;
    const auto map_y = map_y_refined;
    const auto map_s = map_s_refined;

    /* The max s value before wrapping around the track back to 0 */
    const T max_s = 6945.554;

    /* Create the vector of times at which we will evaluate the cost */
    const T dt = 0.02;
    const size_t N = 50;
    using Times = Matrix<T, N, 1>;
    const Times times = Times::LinSpaced(dt, dt * N);

    /* Create the vector of times at which to compute the waypoints */
    const T waypoint_dt = dt;
    const size_t waypoint_N = N;
    using WaypointTimes = Matrix<T, waypoint_N, 1>;
    const WaypointTimes waypoint_times = WaypointTimes::LinSpaced(waypoint_dt, waypoint_dt * waypoint_N);
    cout << "waypoint_times: " << waypoint_times.transpose() << endl;

    /* Define the S and D trajectory parameterizations */
    using S_PARAM = S<T>;
    using D_PARAM = D<T>;
    S_PARAM s_parameterization;
    D_PARAM d_parameterization;

    /* Create the initial guess */
    Matrix<T, Dynamic, 1> solution = Matrix<T, Dynamic, 1>::Zero(
            s_parameterization.N_PARAMS + d_parameterization.N_PARAMS, 1);
    s_parameterization.goodInitialGuess(solution.data());
    d_parameterization.goodInitialGuess(solution.data() + s_parameterization.N_PARAMS);

    /* A list of the previously predicted trajectory times */
    std::queue<Waypoint<T>> waypoints;
    for (size_t i = 0; i < N; ++i) {
        waypoints.push(Waypoint<T>(waypoint_times[i]));
    }

    /* The uWebSocket and the lambda that will run whenever a new message is received */
    uWS::Hub h;
    h.onMessage([&solution,
                        &s_parameterization,
                        &d_parameterization,
                        &waypoints,
                        &times,
                        &waypoint_times,
                        &map_x,
                        &map_y,
                        &map_s](
            uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
            uWS::OpCode opCode) {

        /* "42" at the start of the message means there's a websocket message event.
         * The 4 signifies a websocket message
         * The 2 signifies a websocket event */
        if (length && length > 2 && data[0] == '4' && data[1] == '2') {
            auto s = hasData(data);

            if (s != "") {
                auto j = nlohmann::json::parse(s);
                if (j[0].get<string>() == "telemetry") {

                    /* Main car's localization data */
                    const T car_x = j[1]["x"];
                    const T car_y = j[1]["y"];
                    const T car_s = j[1]["s"];
                    const T car_d = j[1]["d"];
                    const T car_yaw = j[1]["yaw"];
                    const T car_speed_mph = j[1]["speed"];
                    const T car_speed = speedToMetersPerSecond(car_speed_mph);

                    const T car_x_dot = car_speed * cos(car_yaw);
                    const T car_y_dot = car_speed * sin(car_yaw);
                    auto car_sd_dot = getFrenet(car_x_dot, car_y_dot, car_yaw, map_x, map_y);
                    const T car_s_dot = car_sd_dot[0];
                    const T car_d_dot = car_sd_dot[1];

                    /* Show our car's data */
                    if (verbose)
                        cout << "(id,x,y,vx,vy,s,d)" << endl
                             << "("
                             << -1 << ", "
                             << car_x << ", "
                             << car_y << ", "
                             << car_speed * cos(car_yaw) << ", "
                             << car_speed * sin(car_yaw) << ", "
                             << car_s << ", "
                             << car_d << ")" << endl;

                    /* Previous path data given to the Planner */
                    const auto previous_path_x = j[1]["previous_path_x"];
                    const auto previous_path_y = j[1]["previous_path_y"];

                    /* Previous path's end s and d values */
                    const T end_path_s = j[1]["end_path_s"];
                    const T end_path_d = j[1]["end_path_d"];

                    /* Sensor fusion data, which includes a list of all other cars on the same side of the road. */
                    const auto sensor_fusion = j[1]["sensor_fusion"];

                    /* Extract the coefficients from the previous solution */
                    Map<Matrix<T, S_PARAM::N_PARAMS, 1>> s_coeffs(
                            solution.template topRows<S_PARAM::N_PARAMS>().data());
                    Map<Matrix<T, D_PARAM::N_PARAMS, 1>> d_coeffs(
                            solution.template bottomRows<D_PARAM::N_PARAMS>().data());

                    /* Now, we are going to append to the end of the previous solution.
                     * To make sure that we have a nice smooth continuous path,
                     * we want to know what the position and higher derivatives are
                     * so that we can match them.
                     *
                     * To accomplish this, we first pop off all of the "used" waypoints
                     */
                    for (size_t i = 0; i < waypoint_times.size() - previous_path_x.size(); ++i)
                        waypoints.pop();

                    if (waypoints.empty()) {

                        s_parameterization.x0 = car_s;
                        s_parameterization.v0 = car_speed;
                        s_parameterization.a0 = 0;
                        // s_parameterization.j0 = 0;

                        d_parameterization.x0 = car_d;
                        d_parameterization.v0 = 0;
                        d_parameterization.a0 = 0;
                        // d_parameterization.j0 = 0;

                    } else {

                        s_parameterization.x0 = waypoints.back().x;
                        s_parameterization.v0 = waypoints.back().dx;
                        s_parameterization.a0 = waypoints.back().ddx;
                        // s_parameterization.j0 = waypoints.back().dddx;

                        d_parameterization.x0 = waypoints.back().y;
                        d_parameterization.v0 = waypoints.back().dy;
                        d_parameterization.a0 = waypoints.back().ddy;
                        // d_parameterization.j0 = waypoints.back().dddy;
                    }

                    /* Create a list of all of the other cars. We will use this info to prevent collisions.
                     * Note that the parameterizations for our car assume that we are starting from t = 0.
                     * Therefore, when interpolating the positions of the other cars, we
                     * would say that their initial condition data that we received are
                     * 0.02 * number of waypoints in the past. */
                    using Cars = std::vector<Car<T>>;
                    Cars cars;
                    const T t0 = -0.02 * previous_path_x.size();
                    for (size_t i = 0; i < sensor_fusion.size(); ++i) {
                        const T s = sensor_fusion[i][5];
                        const T d = sensor_fusion[i][6];
                        const T vx = sensor_fusion[i][3];
                        const T vy = sensor_fusion[i][4];
                        cars.push_back(Car<T>(t0, s, d, sqrt(vx * vx + vy * vy)));
                        if (verbose)
                            cout << sensor_fusion[i] << endl;
                    }
                    cout << endl;

                    /* Create the new initial guess */
                    s_parameterization.goodInitialGuess(solution.data());
                    d_parameterization.goodInitialGuess(solution.data() + s_parameterization.N_PARAMS);

                    /* Define the cost function */
                    using Cost = TotalCost<T, Times, Cars, S_PARAM, D_PARAM>;
                    Cost total_cost(times, cars, s_parameterization, d_parameterization);
                    Eigen::NumericalDiff<Cost> numDiff(total_cost);
                    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Cost>, T> lm(numDiff);
                    lm.parameters.maxfev = 500;
                    lm.parameters.xtol = 1.0e-10;
                    lm.parameters.ftol = 1.0e-10;

                    /* Optimize */
                    cout << "Solution before: " << solution.transpose() << endl;
                    cout << "initial cost: " << total_cost.eval(solution) << endl;

                    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                    int ret = lm.minimize(solution);
                    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

                    cout << "Solution after:  " << solution.transpose() << endl;
                    cout << "final cost: " << total_cost.eval(solution) << endl;

                    cout << "Time (s): " << duration / 1e3 << endl;
                    cout << "Iterations: " << lm.iter << endl;
                    cout << "Return code: " << ret << endl;

                    /* We are going to append waypoints to the end of the previous path.
                     * But before we do, make sure that there are not any
                     * lingering waypoints that would violate the speed limit. */
                    vector<T> next_x_vals = previous_path_x;
                    vector<T> next_y_vals = previous_path_y;
                    T max_dist = 0;
                    T last_x = next_x_vals.size() > 0 ? next_x_vals[0] : car_x;
                    T last_y = next_y_vals.size() > 0 ? next_y_vals[0] : car_y;
                    for (size_t i = 1; i < next_x_vals.size(); ++i) {
                        const auto dist = norm(next_x_vals[i] - last_x, next_y_vals[i] - last_y);
                        if (dist > max_dist) {
                            max_dist = dist;
                        }
                        last_x = next_x_vals[i];
                        last_y = next_y_vals[i];
                    }
                    cout << "Max distance (m): " << max_dist << ", equivalent mph: " << mpsToMph(max_dist / 0.02)
                         << endl;

                    /* Generate the waypoints based on the s and d trajectories.
                     * We just want to append to the previous waypoint paths
                     * until the list of waypoints is the desired length. */
                    vector<T> next_s, next_d;
                    for (size_t i = 0; i < waypoint_times.size() - previous_path_x.size(); ++i) {

                        /* Save the time. Although it shouldn't be possible, it seems that sometimes
                         * the velocity bound is being exceeded. So we wil sometimes need to decrease the time. */
                        auto t = waypoint_times[i];
                        const T prev_t = i > 0 ? waypoint_times[i - 1] : 0;
                        auto s = s_parameterization.template eval<0>(s_coeffs, t);
                        auto d = d_parameterization.template eval<0>(d_coeffs, t);
                        auto xy = getXY(s, d, map_s, map_x, map_y);
                        const auto last_x = next_x_vals.empty() ? car_x : next_x_vals.back();
                        const auto last_y = next_y_vals.empty() ? car_y : next_y_vals.back();
                        while (!validStep(xy[0] - last_x, xy[1] - last_y)) {
                            cout << "Too fast. Truncating timestep from " << t - prev_t << " to ";
                            t -= 0.01 * (t - prev_t);
                            cout << t - prev_t << endl;
                            s = s_parameterization.template eval<0>(s_coeffs, t);
                            d = d_parameterization.template eval<0>(d_coeffs, t);
                            xy = getXY(s, d, map_s, map_x, map_y);
                        }

                        /* Evaluate the trajectories at the specified times */
                        const auto ds = s_parameterization.template eval<1>(s_coeffs, t);
                        const auto dd = d_parameterization.template eval<1>(d_coeffs, t);
                        const auto dds = s_parameterization.template eval<2>(s_coeffs, t);
                        const auto ddd = d_parameterization.template eval<2>(d_coeffs, t);
                        const auto ddds = s_parameterization.template eval<3>(s_coeffs, t);
                        const auto dddd = d_parameterization.template eval<3>(d_coeffs, t);
                        next_s.push_back(s);
                        next_d.push_back(d);

                        /* Save the new waypoints */
                        waypoints.emplace(Waypoint<T>(0, s, d, ds, dd, dds, ddd, ddds, dddd));

                        /* Convert from Frenet to Cartesian coordinates. */
                        next_x_vals.push_back(xy[0]);
                        next_y_vals.push_back(xy[1]);
                    }

                    if (verbose) {
                        cout << "s: " << Map<Matrix<T, Dynamic, 1>>(next_s.data(), next_s.size()).transpose() << endl;
                        cout << "d: " << Map<Matrix<T, Dynamic, 1>>(next_d.data(), next_d.size()).transpose() << endl;
                        cout << "x: " << next_x_vals.size() << " = "
                             << Map<WaypointTimes>(next_x_vals.data()).transpose()
                             << endl;
                        cout << "y: " << Map<WaypointTimes>(next_y_vals.data()).transpose() << endl;
                        cout << endl << "------------------" << endl;
                    }

                    /* Show the maximum distance between waypoints */
                    max_dist = 0;
                    size_t max_index = 0;
                    last_x = next_x_vals[0];
                    last_y = next_y_vals[0];
                    for (size_t i = 1; i < next_x_vals.size(); ++i) {
                        const auto dist = norm(next_x_vals[i] - last_x, next_y_vals[i] - last_y);
                        if (dist > max_dist) {
                            max_dist = dist;
                            max_index = i;
                        }
                        last_x = next_x_vals[i];
                        last_y = next_y_vals[i];
                    }
                    cout << "Index : " << max_index << ", Previous size: " << previous_path_x.size()
                         << ", max distance (m): " << max_dist << ", equivalent mph: "
                         << mpsToMph(max_dist / 0.02)
                         << endl;

                    /* Send the next set of waypoints */
                    nlohmann::json
                            msg_json;
                    msg_json["next_x"] = next_x_vals;
                    msg_json["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msg_json.dump() + "]";

                    std::this_thread::sleep_for(std::chrono::milliseconds(20));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                }
            } else {
                // Manual driving
                string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    /* We don't need this since we're not using HTTP but if it's removed the program
     * doesn't compile :-( */
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            /* I guess this should be done more gracefully? */
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        cout << "Connected!!!" << endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        cout << "Disconnected" << endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        cout << "Listening to port " << port << endl;
    } else {
        std::cerr << "Failed to listen to port" << endl;
        return -1;
    }
    h.run();
}
