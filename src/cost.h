#ifndef COST_HEADER
#define COST_HEADER

#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"
#include "polynomial.h"
#include "utils.h"

using std::abs;
using std::pow;
using std::cos;
using std::sin;
using std::sqrt;
using std::exp;

using Eigen::Dynamic;
using Eigen::Map;
using Eigen::Matrix;

using std::cout;
using std::endl;

const bool verbose = false;

/**
 * Generic functor
 * See http://eigen.tuxfamily.org/index.php?title=Functors
 * C++ version of a function pointer that stores meta-data about the function
 */
template<typename _Scalar, int NX = Dynamic, int NY = Dynamic>
struct Functor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}

    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }

    int values() const { return m_values; }

};

/* Constants */
template<typename T>
constexpr T ALPHA() {
    return 2;
}

template<typename T>
constexpr T GAMMA() {
    return 0.95;
}

template<typename T>
constexpr T LOW_SCALE() {
    return 1 / GAMMA<T>();
}

template<typename T>
constexpr T HIGH_SCALE() {
    return (ALPHA<T>() - 1) / (1 - GAMMA<T>());
}

template<typename T>
constexpr T HIGH_INTERCEPT() {
    return (1 - GAMMA<T>() * ALPHA<T>()) / (1 - GAMMA<T>());
}

template<typename T>
constexpr T LANES() {
    return 3;
}

template<typename T>
constexpr T LANE_WIDTH() {
    return 4;
}

template<typename T>
constexpr T ROAD_WIDTH() {
    return LANES<T>() * LANE_WIDTH<T>();
}

template<typename T>
constexpr T TIGHT_ROAD_WIDTH() {
    return ROAD_WIDTH<T>() - 1.5;
}

template<typename T>
constexpr T SAFE_S_DISTANCE() {
    return 30;
}

template<typename T>
constexpr T S_BAR() {
    return SAFE_S_DISTANCE<T>() / sqrt(2 - 2 * GAMMA<T>());
}

template<typename T>
constexpr T SAFE_D_DISTANCE() {
    return 3;
}

template<typename T>
constexpr T D_BAR() {
    return SAFE_D_DISTANCE<T>() / sqrt(2 - 2 * GAMMA<T>());
}

template<typename T>
constexpr T SIGMA() {
    return 10;
}

template<typename T>
constexpr T REFERENCE_VELOCITY() {
    return 21;
}

template<typename T>
constexpr T REFERENCE_VELOCITY_SQ() {
    return REFERENCE_VELOCITY<T>() * REFERENCE_VELOCITY<T>();
}

/**
 * The prototypical cost function that we use for the critical constraints.
 */
template<typename T>
constexpr T MAX_POWER() {
    return 2.5;
}

/** The prototypical cost function that we use for the critical constraints. */
template<typename T, typename Z>
constexpr T proto(const Z &z) {
    return pow(ALPHA<T>(), (z - GAMMA<T>()) / (1 - GAMMA<T>()));
//    return z < 0 ? 0
//                 : z < GAMMA<T>() ? z * LOW_SCALE<T>()
//                                  : z * HIGH_SCALE<T>() + HIGH_INTERCEPT<T>();
//    return z < GAMMA<T>() ? exp(z - GAMMA<T>())
//                          : z * HIGH_SCALE<T>() + HIGH_INTERCEPT<T>();
}

/** A constraint meant to ensure that we stay on the correct side of the road */
template<typename T, typename S, typename D>
constexpr T road(const S &s, const D &d) {
    return proto<T>(d / TIGHT_ROAD_WIDTH<T>()) + proto<T>(1 - d / TIGHT_ROAD_WIDTH<T>())
           + abs((d - ROAD_WIDTH<T>() / 2) / (ROAD_WIDTH<T>() / 2));
}

/** Critical constraints for the velocity, acceleration, and jerk. */
template<typename T, typename S_DOT, typename D_DOT>
T velocity(const S_DOT &s_dot, const D_DOT &d_dot) {
    const auto v_sq = pow(s_dot, 2) + pow(d_dot, 2);
    const auto cost = 10 * exp(-7 * v_sq / REFERENCE_VELOCITY_SQ<T>());
    return cost + (v_sq > 420) ?: proto<T>(v_sq / REFERENCE_VELOCITY_SQ<T>());
}

template<typename T, typename S_DDOT, typename D_DDOT>
constexpr T acceleration(const S_DDOT &s_ddot, const D_DDOT &d_ddot) {
    return proto<T>((pow(s_ddot, 2) + pow(d_ddot, 2)) / 100);
}

template<typename T, typename S_DDDOT, typename D_DDDOT>
constexpr T jerk(const S_DDDOT &s_dddot, const D_DDDOT &d_dddot) {
    return proto<T>((pow(s_dddot, 2) + pow(d_dddot, 2)) / 100);
}

/** Collision avoidance cost */
template<typename T, typename S, typename D, typename S_NEARBY, typename D_NEARBY>
T collisionAvoidance(const S &s, const D &d, const S_NEARBY &s_nearby, const D_NEARBY &d_nearby) {
    T cost = 0;
    for (size_t i = 0; i < s_nearby.size(); ++i)
        cost += pow(proto<T>(1 - pow((s_nearby[i] - s) / S_BAR<T>(), 2))
                    * proto<T>(1 - pow((d_nearby[i] - d) / D_BAR<T>(), 2)), 3);
    return cost;
}

/** Lane-keeping cost */
template<typename T, typename S, typename D>
constexpr T lane(const S &s, const D &d) {
//    return pow(cos(M_PI * d / LANE_WIDTH<T>()), 2);
    return sqrt(abs(cos(M_PI * d / LANE_WIDTH<T>())));
}

/** Reference velocity cost */
template<typename T, typename S_DOT, typename D_DOT>
constexpr T referenceVelocity(const S_DOT &s_dot, const D_DOT &d_dot) {
    return 1 - exp(-pow((sqrt(pow(s_dot, 2) + pow(d_dot, 2)) - REFERENCE_VELOCITY<T>()) / SIGMA<T>(), 2));
}

/**
 * The parameterization for the S trajectory.
 */
template<typename T>
struct S {

    static constexpr size_t N_PARAMS = 1;
    T x0 = 0, v0 = 0, a0 = 0;
    const T min_v0 = 5;
    const T default_b = 1;
    const T default_g = 2;

    /** Fill x with a reasonable initial guess for the parameters */
    void goodInitialGuess(T *x) {
        x[0] = v0 > min_v0 ? v0 : min_v0;
        if (N_PARAMS > 1)
            x[1] = default_b;
        if (N_PARAMS > 2)
            x[2] = default_g;
    }

    template<size_t deriv, typename Params>
    T eval(const Params &params, T t) const {
        auto v = abs(params(0));
        auto b = N_PARAMS > 1 ? abs(params(1)) : default_b;
        auto g = N_PARAMS > 2 ? abs(params(2)) : default_g;
        b = b > 10 ? 10 : b;
        g = g > 10 ? 10 : g;
        v = v > REFERENCE_VELOCITY<T>() ? REFERENCE_VELOCITY<T>() : v;

        auto gamma = -(a0 - b * v + b * v0) / (g * (b - g));
        auto beta = -(v0 - v + gamma * g) / b;
        auto delta = x0 - beta - gamma;
        if (deriv == 0) {
            return delta + v * t + beta * exp(-b * t) + gamma * exp(-g * t);
        } else if (deriv == 1) {
            return v - b * beta * exp(-b * t) - g * gamma * exp(-g * t);
        } else {
            return pow(-1, deriv) *
                   (pow(b, deriv) * beta * exp(-b * t) + pow(g, deriv) * gamma * exp(-g * t));
        }
        return 0;
    }
};

/**
 * The parameterization for the D trajectory
 */
template<typename T>
struct D {

    static constexpr size_t N_PARAMS = 1;
    T x0 = 2, v0 = 0, a0 = 0;
    const T default_b = 1.3;
    const T default_g = 1.6;

    /** Fill x with a reasonable initial guess for the parameters */
    void goodInitialGuess(T *x) {
        x[0] = x0;
        if (N_PARAMS > 1)
            x[1] = default_b;
        if (N_PARAMS > 2)
            x[2] = default_g;
    }

    template<size_t deriv, typename Params>
    T eval(const Params &params, T t) const {
        auto x = params(0);
        auto b = N_PARAMS > 1 ? abs(params(1)) : default_b;
        auto g = N_PARAMS > 2 ? abs(params(2)) : default_g;
        auto delta = -(a0 + v0 * b + v0 * g - x * b * g + x0 * b * g) /
                     ((b - 1) * (g - 1));
        auto beta = -(a0 + delta - delta * g + v0 * g) / (b * (b - g));
        auto gamma = (a0 + delta - delta * b + v0 * b) / (g * (b - g));
        if (verbose)
            cout << "delta=" << delta << ", beta=" << beta << ", gamma=" << gamma << ", b=" << b
                 << ", g=" << g << endl;
        if (deriv == 0) {
            return x0 - delta * (exp(-t) - 1) - beta * (exp(-t * b) - 1) - gamma * (exp(-t * g) - 1);
        } else {
            return -pow(-1, deriv) * (
                    delta * exp(-t) + beta * pow(b, deriv) * exp(-t * b) +
                    gamma * pow(g, deriv) * exp(-t * g));
        }
    }
};

/**
* A parameterization for a X or Y trajectory
*/
template<typename T>
struct XY {

    static constexpr size_t N_PARAMS = 3;
    T x0 = 0, v0 = 0, a0 = 0;

    /** Fill x with a reasonable initial guess for the parameters */
    void goodInitialGuess(T *x) {
        x[0] = v0;
        x[1] = 2;
        x[2] = 3;
    }

    template<size_t deriv, typename Params>
    T eval(const Params &params, T t) const {
        auto v = params(0);
        auto b = abs(params(1));
        auto g = abs(params(2));

        auto gamma = -(a0 - b * v + b * v0) / (g * (b - g));
        auto beta = -(v0 - v + gamma * g) / b;

        if (deriv == 0) {
            return x0 - beta - gamma + v * t + beta * exp(-b * t) + gamma * exp(-g * t);
        } else if (deriv == 1) {
            return v - b * beta * exp(-b * t) - g * gamma * exp(-g * t);
        } else {
            return pow(-1, deriv) * (
                    beta * pow(b, deriv) * exp(-t * b) +
                    gamma * pow(g, deriv) * exp(-t * g));
        }
    }
};

/** The total cost function */
template<typename T, typename Times, typename Cars, typename S_PARAM, typename D_PARAM>
struct TotalCost : Functor<T> {

    const Times times;
    Cars cars;
    S_PARAM s_parameterization;
    D_PARAM d_parameterization;

    TotalCost(const Times &times,
              const Cars &cars,
              const S_PARAM &s_parameterization,
              const D_PARAM &d_parameterization)
            : Functor<T>(S_PARAM::N_PARAMS + D_PARAM::N_PARAMS, 7 * times.size()),
              times(times),
              cars(cars),
              s_parameterization(s_parameterization),
              d_parameterization(d_parameterization) {
    }

    int operator()(const Matrix<T, Dynamic, 1> &x, Matrix<T, Dynamic, 1> &f) const {

        /* The first degree + 1 terms correspond to the s polynomial.
         * The remaining terms correspond to the d polynomial */
        const Matrix<T, S_PARAM::N_PARAMS, 1> s_coeffs = x.template topRows<S_PARAM::N_PARAMS>();
        const Matrix<T, D_PARAM::N_PARAMS, 1> d_coeffs = x.template bottomRows<D_PARAM::N_PARAMS>();
        for (size_t i = 0; i < times.size(); ++i) {

            /* Evaluate the parameterizations and their derivatives at t */
            const auto t = times(i);
            const T s = s_parameterization.template eval<0>(s_coeffs, t);
            const T d = d_parameterization.template eval<0>(d_coeffs, t);
            const T s_dot = s_parameterization.template eval<1>(s_coeffs, t);
            const T d_dot = d_parameterization.template eval<1>(d_coeffs, t);
            const T s_ddot = s_parameterization.template eval<2>(s_coeffs, t);
            const T d_ddot = d_parameterization.template eval<2>(d_coeffs, t);
            const T s_dddot = s_parameterization.template eval<3>(s_coeffs, t);
            const T d_dddot = d_parameterization.template eval<3>(d_coeffs, t);

            /* Evaluate the positions of the other vehicles at this time */
            std::vector<T> s_nearby(cars.size()), d_nearby(cars.size());
            for (size_t i = 0; i < cars.size(); ++i) {
                std::tie(s_nearby[i], d_nearby[i]) = cars[i](t);
            }

            f(7 * i + 0) = road<T>(s, d);
            f(7 * i + 1) = velocity<T>(s_dot, d_dot);
            f(7 * i + 2) = acceleration<T>(s_ddot, d_ddot);
            f(7 * i + 3) = jerk<T>(s_dddot, d_dddot);
            f(7 * i + 4) = collisionAvoidance<T>(s, d, s_nearby, d_nearby);
            f(7 * i + 5) = lane<T>(s, d);
            f(7 * i + 6) = 0; //ALPHA<T>() * referenceVelocity<T>(s_dot, d_dot);
        }
//        cout << "cost: " << f.array().square().sum() << endl;
//        cout << Map<Matrix<T, Dynamic, Dynamic>>(f.data(), 7, times.size()).rowwise().sum().transpose() << endl;
        return 0;
    }

    T eval(const Matrix<T, Dynamic, 1> &x) {
        Matrix<T, Dynamic, 1> f(this->values());
        this->operator()(x, f);
        return f.array().square().sum();
    }
};

/** The total cost function */
template<typename T, typename Times, typename S_NEARBY, typename D_NEARBY, typename X_PARAM, typename Y_PARAM>
struct TotalCostXY : Functor<T> {

    const T theta;
    const Times times;
    const std::vector<T> maps_x;
    const std::vector<T> maps_y;
    S_NEARBY s_nearby;
    D_NEARBY d_nearby;
    X_PARAM x_parameterization;
    Y_PARAM y_parameterization;

    TotalCostXY(const T theta,
                const Times &times,
                const S_NEARBY &s_nearby,
                const D_NEARBY &d_nearby,
                const X_PARAM &x_parameterization,
                const Y_PARAM &y_parameterization,
                const std::vector<T> &maps_x,
                const std::vector<T> &maps_y)
            : Functor<T>(X_PARAM::N_PARAMS + Y_PARAM::N_PARAMS, 7 * times.size()),
              theta(theta),
              times(times),
              s_nearby(s_nearby),
              d_nearby(d_nearby),
              x_parameterization(x_parameterization),
              y_parameterization(y_parameterization),
              maps_x(maps_x),
              maps_y(maps_y) {
    }

    int operator()(const Matrix<T, Dynamic, 1> &solution, Matrix<T, Dynamic, 1> &f) const {

        const Matrix<T, X_PARAM::N_PARAMS, 1> x_coeffs = solution.template topRows<X_PARAM::N_PARAMS>();
        const Matrix<T, Y_PARAM::N_PARAMS, 1> y_coeffs = solution.template bottomRows<Y_PARAM::N_PARAMS>();
        T x, y, dx, dy, s, d;
        f = Matrix<T, Dynamic, 1>::Zero(7 * times.size());
        for (size_t i = 0; i < times.size(); ++i) {

            /* Evaluate the parameterizations and their derivatives at t */
            const auto t = times(i);
            x = x_parameterization.template eval<0>(x_coeffs, t);
            y = y_parameterization.template eval<0>(y_coeffs, t);
            dx = x_parameterization.template eval<1>(x_coeffs, t);
            dy = y_parameterization.template eval<1>(y_coeffs, t);
            const T ddx = x_parameterization.template eval<2>(x_coeffs, t);
            const T ddy = y_parameterization.template eval<2>(y_coeffs, t);
            const T dddx = x_parameterization.template eval<3>(x_coeffs, t);
            const T dddy = y_parameterization.template eval<3>(y_coeffs, t);

            /* These can be defined in terms of sd or xy */
            f(7 * i + 0) = referenceVelocity<T>(dx, dy);
            f(7 * i + 1) = velocity<T>(dx, dy);
            f(7 * i + 2) = acceleration<T>(ddx, ddy);
            f(7 * i + 3) = jerk<T>(dddx, dddy);

            /* These are tied to the Frenet coordinates, so we will need to convert them. */
            // TODO Theta?
            auto sd = getFrenet(x, y, static_cast<T>(pi), maps_x, maps_y);
            s = sd[0];
            d = sd[1];
//            f(7 * i + 4) = collisionAvoidance<T>(s, d, s_nearby, d_nearby);
            f(7 * i + 5) = lane<T>(s, d);
//            f(7 * i + 6) = road<T>(s, d);
        }
        cout << "s: " << s << ", d: " << d << ", dx: " << dx << ", dy: " << dy << endl;
        cout << "cost: " << f.array().square().sum() << endl;
        cout << "evaluating at : " << solution.transpose() << endl;
        cout << Map<Matrix<T, Dynamic, Dynamic>>(f.data(), 7, times.size()).rowwise().mean().transpose() << endl;
        return 0;
    }
};

#endif /* COST_HEADER */