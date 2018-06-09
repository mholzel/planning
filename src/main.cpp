#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "Eigen/Dense"
#include "json.hpp"
#include "polynomial.h"
#include "cost.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

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
    const auto map_x = map_x_;
    const auto map_y = map_y_;
    const auto map_s = map_s_;

    /* The max s value before wrapping around the track back to 0 */
    const T max_s = 6945.554;

    /* Create the vector of times at which we will evaluate the cost */
    const T dt = 0.02;
    const size_t N = 50;
    using Times = Matrix<T, N, 1>;
    const Times times = Times::LinSpaced(dt, dt * N);

    /* Create the vector of times at which we will evaluate the cost */
    const T waypoint_dt = dt;
    const size_t waypoint_N = N;
    using WaypointTimes = Matrix<T, waypoint_N, 1>;
    const WaypointTimes waypoint_times = WaypointTimes::LinSpaced(waypoint_dt, waypoint_dt * waypoint_N);
    cout << "waypoint_times: " << waypoint_times.transpose() << endl;

    /* Define the X and Y trajectory parameterizations */
    using X_PARAM = S<T>;
    using Y_PARAM = S<T>;
    X_PARAM x_parameterization;
    Y_PARAM y_parameterization;

    /* Create the initial guess */
    Matrix<T, Dynamic, 1> solution = Matrix<T, Dynamic, 1>::Zero(
            x_parameterization.N_PARAMS + y_parameterization.N_PARAMS, 1);
    x_parameterization.goodInitialGuess(solution.data());
    y_parameterization.goodInitialGuess(solution.data() + x_parameterization.N_PARAMS);

    size_t its;

    /* The uWebSocket and the lambda that will run whenever a new message is received */
    uWS::Hub h;
    h.onMessage([&its, &solution,
                        &x_parameterization,
                        &y_parameterization,
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

                    /* Show the car's data */
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
                    using NEARBY = Matrix<T, Dynamic, 1>;
                    NEARBY x_nearby = NEARBY::Zero(sensor_fusion.size());
                    NEARBY y_nearby = NEARBY::Zero(sensor_fusion.size());
                    for (size_t i = 0; i < sensor_fusion.size(); ++i) {
                        x_nearby(i) = sensor_fusion[i][5];
                        y_nearby(i) = sensor_fusion[i][6];
                        cout << sensor_fusion[i] << endl;
                    }
                    // cout << "x: " << x_nearby.transpose() << endl;
                    // cout << "y: " << y_nearby.transpose() << endl;
                    cout << endl;

                    vector<T> next_x_vals(waypoint_times.size());
                    vector<T> next_y_vals(waypoint_times.size());
                    if (true) {

                        /* Extract the coefficients from the previous solution */
                        Map<Matrix<T, X_PARAM::N_PARAMS, 1>> x_coeffs(
                                solution.template topRows<X_PARAM::N_PARAMS>().data());
                        Map<Matrix<T, Y_PARAM::N_PARAMS, 1>> y_coeffs(
                                solution.template bottomRows<Y_PARAM::N_PARAMS>().data());

                        /* First, set the initial position to match the car's and
                         * the velocity, acceleration, and jerk to match the old trajectory.
                         * We do that by evaluating the old parameterizations at the
                         * elapsed time. */
                        const size_t np = previous_path_x.size();
                        const T t = np == 0 ? 0 : 0.02 * (waypoint_times.size() - np);

                        cout << "Elapsed time: " << t << endl;

                        const T x0 = car_x;
                        const T dx0 = false ? x_parameterization.template eval<1>(x_coeffs, t)
                                            : car_x_dot;
                        const T ddx0 = x_parameterization.template eval<2>(x_coeffs, t);
                        x_parameterization.x0 = x0;
                        x_parameterization.v0 = dx0;
                        x_parameterization.a0 = ddx0;

                        const T y0 = car_y;
                        const T dy0 = false ? y_parameterization.template eval<1>(y_coeffs, t)
                                            : car_y_dot;
                        const T ddy0 = y_parameterization.template eval<2>(y_coeffs, t);
                        y_parameterization.x0 = y0;
                        y_parameterization.v0 = dy0;
                        y_parameterization.a0 = ddy0;

                        cout << x_parameterization.template eval<1>(x_coeffs, t) << ", "
                             << car_x_dot << ", "
                             << y_parameterization.template eval<1>(y_coeffs, t) << ", "
                             << car_y_dot << endl;
                        cout << "X initial conditions: " << x0 << ", " << dx0 << ", " << ddx0 << endl;
                        cout << "Y initial conditions: " << y0 << ", " << dy0 << ", " << ddy0 << endl;

                        /* Create the new initial guess */
                        x_parameterization.goodInitialGuess(solution.data());
                        y_parameterization.goodInitialGuess(solution.data() + x_parameterization.N_PARAMS);

                        /* Define the cost function */
                        using Cost = TotalCostXY<T, Times, NEARBY, NEARBY, X_PARAM, Y_PARAM>;
                        Cost total_cost(car_yaw,
                                        times,
                                        x_nearby,
                                        y_nearby,
                                        x_parameterization,
                                        y_parameterization,
                                        map_x,
                                        map_y);
                        Eigen::NumericalDiff<Cost> numDiff(total_cost);
                        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Cost>, T> lm(numDiff);
                        lm.parameters.maxfev = 1000;
                        lm.parameters.xtol = 1.0e-10;
                        lm.parameters.ftol = 1.0e-10;

                        /* Optimize */
                        cout << "Solution before: " << solution.transpose() << endl;
                        cout << "x_coeffs: " << x_coeffs.transpose() << endl;
                        cout << "y_coeffs: " << y_coeffs.transpose() << endl;

                        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                        int ret = lm.minimize(solution);
                        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

                        cout << "Solution after:  " << solution.transpose() << endl;
                        cout << "x_coeffs: " << x_coeffs.transpose() << endl;
                        cout << "y_coeffs: " << y_coeffs.transpose() << endl;

                        cout << "Time (s): " << duration / 1e6 << endl;
                        cout << "Iterations: " << lm.iter << endl;
                        cout << "Return code: " << ret << endl;

                        /* Generate the waypoints based on the x and y trajectories */
                        for (size_t i = 0; i < waypoint_times.size(); ++i) {

                            /* Evaluate the trajectories at the specified times */
                            const auto x = x_parameterization.template eval<0>(x_coeffs, waypoint_times[i]);
                            const auto y = y_parameterization.template eval<0>(y_coeffs, waypoint_times[i]);
                            next_x_vals[i] = x;
                            next_y_vals[i] = y;
                        }

                        cout << "x: " << Eigen::Map<WaypointTimes>(next_x_vals.data()).transpose() << endl;
                        cout << "y: " << Eigen::Map<WaypointTimes>(next_y_vals.data()).transpose() << endl;
                        cout << endl << "------------------" << endl;
                    }

                    /* Send the next set of waypoints */
                    nlohmann::json msg_json;
                    msg_json["next_x"] = next_x_vals;
                    msg_json["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msg_json.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    if (++its < 3)
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
