#ifndef WAYPOINT_HEADER
#define WAYPOINT_HEADER

/** This class saves the time and x,y position of each waypoint.
 * Additionally, some schemas may call for higher derivatives.
 * This class can hold everything up to the jerk derivative.
 * Note that everything that you don't specify in the constructor
 * will be set to 0. */
template<typename T>
struct Waypoint {

    T t, x, y, dx, dy, ddx, ddy, dddx, dddy;

    Waypoint(T t, T x, T y, T dx, T dy, T ddx, T ddy, T dddx, T dddy)
            : t(t), x(x), y(y), dx(dx), dy(dy), ddx(ddx), ddy(ddy), dddx(dddx), dddy(dddy) {
    }

    Waypoint(T t, T x, T y, T dx, T dy, T ddx, T ddy)
            : Waypoint(t, x, y, dx, dy, ddx, ddy, 0, 0) {
    }

    Waypoint(T t, T x, T y, T dx, T dy)
            : Waypoint(t, x, y, dx, dy, 0, 0) {
    }

    Waypoint(T t, T x, T y)
            : Waypoint(t, x, y, 0, 0) {
    }

    Waypoint(T t)
            : Waypoint(t, 0, 0) {
    }

    Waypoint()
            : Waypoint(0) {
    }
};

#endif /* WAYPOINT_HEADER */
