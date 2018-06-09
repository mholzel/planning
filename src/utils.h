#ifndef UTILS_HEADER
#define UTILS_HEADER

#include "spline.h"

const double pi = M_PI;

template<typename T>
T deg2rad(T x) { return x * pi / 180; }

template<typename T>
T rad2deg(T x) { return x * 180 / pi; }

template<typename T>
T speedToMetersPerSecond(T mph) {
    return mph * 0.44704;
}

/** Meters per second to miles per hour. */
template<typename T>
T mpsToMph(T mps) {
    return mps / 0.44704;
}

/** miles per hour to meters per second. */
template<typename T>
T mphToMps(T mph) {
    return mph * 0.44704;
}

/** Rotate the coordinates (x,y) by "yaw" radians.  */
template<typename T>
std::tuple<T, T> rotate(T x, T y, T yaw) {
    return std::make_tuple(cos(yaw) * x - sin(yaw) * y, sin(yaw) * x + cos(yaw) * y);
};

template<typename T>
T distance(T x1, T y1, T x2, T y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/**
 * Checks if the SocketIO event has JSON data.
 * If there is data the JSON object in string format will be returned,
 * else the empty string "" will be returned.
 */
std::string hasData(std::string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != std::string::npos) {
        return "";
    } else if (b1 != std::string::npos && b2 != std::string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

template<typename X, typename Y, typename VX, typename VY>
int closestWaypoint(X x, Y y, const std::vector<VX> &maps_x, const std::vector<VY> &maps_y) {

    if (maps_x.size() == 0) {
        return 0;
    } else {
        int closest_waypoint = 0;
        auto closestLen = distance(x, y, maps_x[0], maps_y[0]);
        for (int i = 1; i < maps_x.size(); i++) {
            auto dist = distance(x, y, maps_x[i], maps_y[i]);
            if (dist < closestLen) {
                closestLen = dist;
                closest_waypoint = i;
            }
        }
        return closest_waypoint;
    }
}

template<typename X, typename Y, typename T, typename VX, typename VY>
int nextWaypoint(X x, Y y, T theta, const std::vector<VX> &maps_x, const std::vector<VY> &maps_y) {

    int closest_waypoint = closestWaypoint(x, y, maps_x, maps_y);

    VX map_x = maps_x[closest_waypoint];
    VY map_y = maps_y[closest_waypoint];

    auto heading = atan2((map_y - y), (map_x - x));

    auto angle = fabs(theta - heading);
    angle = std::min(2 * pi - angle, angle);

    if (angle > pi / 4) {
        closest_waypoint++;
        if (closest_waypoint == maps_x.size()) {
            closest_waypoint = 0;
        }
    }

    return closest_waypoint;
}

/**
 * Transform from Cartesian x,y coordinates to Frenet s,d coordinates
 */
template<typename T, typename VX, typename VY>
std::vector<T> getFrenet(T x, T y, T theta, const std::vector<VX> &maps_x, const std::vector<VY> &maps_y) {

    int next_wp = closestWaypoint(x, y, maps_x, maps_y) + 1;
    if (next_wp == maps_x.size()) {
        next_wp = 0;
    }
    int prev_wp = next_wp - 1;
    if (prev_wp < 0) {
        prev_wp = maps_x.size() - 1;
    }

    VX n_x = maps_x[next_wp] - maps_x[prev_wp];
    VY n_y = maps_y[next_wp] - maps_y[prev_wp];
    T x_x = x - maps_x[prev_wp];
    T x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    T proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    T proj_x = proj_norm * n_x;
    T proj_y = proj_norm * n_y;

    auto frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    VX center_x = 1000 - maps_x[prev_wp];
    VY center_y = 2000 - maps_y[prev_wp];
    auto centerToPos = distance(center_x, center_y, x_x, x_y);
    auto centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef) {
        frenet_d *= -1;
    }

    // calculate s value
    T frenet_s = 0;
    for (int i = 0; i < prev_wp; i++) {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(static_cast<T>(0), static_cast<T>(0), proj_x, proj_y);

    return {frenet_s, frenet_d};

}

/**
 * Transform from Cartesian x,y coordinates to Frenet s,d coordinates
 */
template<typename T>
struct FrenetSpline {

    /** This spline is used simply to generate the map's y-value given a value of x */
    tk::spline y_spline;

    /** This spline tells us the equivalent s-coordinate of an x-value */
    tk::spline s_spline;

    FrenetSpline(const std::vector<T> &map_x,
                 const std::vector<T> &map_y) {
        y_spline.set_points(map_x, map_y);
    }

    std::vector<T> getFrenet(T x, T y) {

        const T dx = 1;
        const T y_minus_dy = y_spline(x - dx);
        const T n_y = y_spline(x) - y_minus_dy;
        const T dy = y - y_minus_dy;

        // find the projection of x onto n
        const T proj_norm = (dx * dx + dy * n_y) / (dx * dx + n_y * n_y);
        const T proj_x = proj_norm * dx;
        const T proj_y = proj_norm * n_y;

        auto frenet_d = distance(dx, dy, proj_x, proj_y);

        //see if d value is positive or negative by comparing it to a center point

        T center_x = 1000 - x + dx;
        T center_y = 2000 - y + dy;
        auto centerToPos = distance(center_x, center_y, dx, dy);
        auto centerToRef = distance(center_x, center_y, proj_x, proj_y);

        if (centerToPos <= centerToRef) {
            frenet_d *= -1;
        }

        // calculate s value
        T frenet_s = 0;
        // TODO prev_wp
//        for (int i = 0; i < prev_wp; i++) {
//            frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
//        }
        frenet_s += distance(static_cast<T>(0), static_cast<T>(0), proj_x, proj_y);

        return {frenet_s, frenet_d};
    }
};

/**
 * Transform from Frenet s,d coordinates to Cartesian x,y
 */
template<typename T, typename V>
std::vector<T>
getXY(T s, T d, const std::vector<V> &maps_s, const std::vector<V> &maps_x, const std::vector<V> &maps_y) {

    int prev_wp = -1;
    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    auto heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    auto seg_s = (s - maps_s[prev_wp]);

    auto seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    auto seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    auto perp_heading = heading - pi / 2;

    T x = seg_x + d * cos(perp_heading);
    T y = seg_y + d * sin(perp_heading);

    return {x, y};

}

/** Return the next value of type T from the specified istringstream */
template<typename T>
T next(std::istringstream &iss) {
    T value;
    iss >> value;
    return value;
}

/** Convert a trajectory to waypoints */
template<typename Poly, typename T, typename V>
std::tuple<std::vector<typename Poly::Scalar>, std::vector<typename Poly::Scalar>>
generateWaypoints(const Poly &s_coeffs,
                  const Poly &d_coeffs,
                  const T &times,
                  const std::vector<V> &maps_s,
                  const std::vector<V> &maps_x,
                  const std::vector<V> &maps_y) {

    std::vector<typename Poly::Scalar> x(times.size()), y(times.size());
    for (size_t i = 0; i < times.size(); ++i) {

        /* Evaluate the trajectories at the specified times */
        const auto s = Polynomial::eval<0>(s_coeffs, times[i]);
        const auto d = Polynomial::eval<0>(d_coeffs, times[i]);

        /* Convert from Frenet to Cartesian coordinates. */
        const auto xy = getXY(s, d, maps_s, maps_x, maps_y);
        x[i] = xy[0];
        y[i] = xy[1];
    }
    return std::make_tuple(x, y);
};

#endif /* UTILS_HEADER */
