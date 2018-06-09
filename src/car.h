#ifndef CAR_HEADER
#define CAR_HEADER

#include <tuple>

/** This class is used to predict a car's (s,d) position at a future time. */
template<typename T>
struct Car {

    /* The car's position and velocities at time t0 */
    const T t0, s0, d0, ds0, dd0;

    /* Constructors */
    Car(const T t0, const T s0, const T d0, const T ds0, const T dd0)
            : t0(t0), s0(s0), d0(d0), ds0(ds0), dd0(dd0) {
    }

    Car(const T t0, const T s0, const T d0, const T v0)
            : Car(t0, s0, d0, v0, 0) {
    }

    /* Get the (s,d) position at a future time. */
    template<typename Time>
    std::tuple<T, T> operator()(const Time t) const {
        return std::make_tuple(s0 + (t - t0) * ds0, d0 + (t - t0) * dd0);
    }
};

#endif /* CAR_HEADER */
