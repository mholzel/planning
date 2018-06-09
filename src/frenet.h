#ifndef FRENET_HEADER
#define FRENET_HEADER

#include <iostream>
#include <vector>
#include "spline.h"

/* Given a map of (x,y) coordinates, this
 * class uses splines to create a smooth polar representation
 * where you can query for the (x,y) coordinates of a
 * point at a given angle along the curve,
 * or for the s-distance at that angle.
 * This is primarily useful for refining a coarse map.
 * For instance, suppose that the map is a circle, but we only
 * have a very small number of points on this circle.
 * Then it is hard to interpolate between points.
 * So we could refine the map by using a spline to generate
 * points for all of the angles -pi...pi with a smaller step
 * than we originally had. */
template<typename T>
struct Refiner {

    Refiner(const std::vector<T> &x, const std::vector<T> y, const std::vector<T> s) {

        // TODO Check that x, y, and s are the same length

        /* First, compute the mean x and y coordinates */
        x0 = 0, y0 = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            x0 += x[i] / x.size();
            y0 += y[i] / x.size();
        }

        /* Now, find the map point with the smallest s value. When computing the angles of
         * each point, we want that to have the smallest angle.
         * 2 * pi should be added to all of the smaller angles. */
        size_t ind = 0;
        T min = 1e10;
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] < min) {
                min = s[i];
                ind = i;
            }
        }
        theta0 = std::atan2(y[ind] - y0, x[ind] - x0);

        /* Now compute the angle from the (x,y) data
         * and construct a vector of tuples containing (angle,x,y,s) */
        std::vector<std::tuple<T, T, T, T>> polar(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            const T dx = x[i] - x0;
            const T dy = y[i] - y0;
            T theta = std::atan2(dy, dx);
            if (theta < theta0)
                theta += 2 * pi;
            polar[i] = std::make_tuple(theta, x[i], y[i], s[i]);
        }

        /* Sort the vector of tuples by angle */
        std::sort(polar.begin(), polar.end());

        /* Now let's add n points from the end to the beginning, and vice versa,
         * so that we don't have discontinuities at the intersections. To do that, first
         * get the angles associated with the first and last points */
        const size_t n = 1;
        decltype(polar) back(n), front(n);
        auto first = polar.front();
        auto last = polar.back();
        using std::get;
        for (size_t i = 0; i < n; ++i) {

            const size_t I = n - 1 - i;
            back[i] = polar[i];
            front[I] = polar[polar.size() - 1 - i];

            /* Fix the angles of the points we are inserting */
            std::get<0>(back[i]) += 2 * pi;
            std::get<0>(front[I]) -= 2 * pi;

            /* Fix the s-values */
            auto dx = std::get<1>(back[i]) - std::get<1>(last);
            auto dy = std::get<2>(back[i]) - std::get<2>(last);
            std::get<3>(back[i]) = std::get<3>(last) + std::sqrt(dx * dx + dy * dy);
            last = back[i];

            dx = std::get<1>(front[I]) - std::get<1>(first);
            dy = std::get<2>(front[I]) - std::get<2>(first);
            std::get<3>(front[I]) = std::get<3>(first) - std::sqrt(dx * dx + dy * dy);
            first = front[I];
        }

        /* Copy polar to the back of the front. Then reassign */
        for (const auto &i : polar)
            front.push_back(i);
        polar = front;

        /* Now push the new back onto polar. */
        for (const auto &i : back)
            polar.push_back(i);

        /* Now split the sorted polar coordinates into theta, r, and s vectors */
        std::vector<T> theta(polar.size()), r(polar.size()), s_map(polar.size());
        for (size_t i = 0; i < polar.size(); ++i) {
            theta[i] = std::get<0>(polar[i]);
            const T dx = std::get<1>(polar[i]) - x0;
            const T dy = std::get<2>(polar[i]) - y0;
            r[i] = std::sqrt(dx * dx + dy * dy);
            s_map[i] = std::get<3>(polar[i]);
        }

        /* To be able to interpolate, we need to convert to polar */
        f.set_points(theta, r);
        g.set_points(theta, s_map);
    }

    /** Get the (x,y) and s coordinates of the path at a specific angle */
    std::tuple<T, T, T> operator()(T theta) {

        /* First, convert the angle to the range [-pi,pi] */
        theta = std::atan2(std::sin(theta), std::cos(theta));

        /* Now, if the angle is less than theta0, add 2 pi since our splines map theta
         * in the [theta0, theta0 + 2*pi] */
        if (theta < theta0)
            theta += 2 * pi;
        const T r = f(theta);
        const T s = g(theta);
        return std::make_tuple(x0 + r * cos(theta), y0 + r * sin(theta), s);
    }

    /** Get the (x,y) and s coordinates of the path at the specified angles */
    std::tuple<std::vector<T>, std::vector<T>> operator()(const std::vector<T> &theta) {
        std::vector<T> x(theta.size()), y(theta.size()), s(theta.size());
        for (size_t i = 0; i < theta.size(); ++i)
            std::tie(x[i], y[i], s[i]) = this->operator()(theta[i]);
        return std::make_tuple(x, y, s);
    }

    /** Generate the specified number of (x,y) and s coordinates which are linearly spaced according to angle. */
    std::tuple<std::vector<T>, std::vector<T>, std::vector<T>> linSpaced(size_t n) {
        std::vector<T> x(n), y(n), s(n);
        for (size_t i = 0; i < n; ++i)
            std::tie(x[i], y[i], s[i]) = this->operator()(theta0 + 2 * pi * i / (n - 1));
        return std::make_tuple(x, y, s);
    }

private:

    tk::spline f, g;
    double x0, y0, theta0;

};

#endif /* FRENET_HEADER */
