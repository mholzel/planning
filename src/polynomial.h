#ifndef POLYNOMIAL_HEADER
#define POLYNOMIAL_HEADER

#include "Eigen/Dense"
#include <iostream>

namespace {

    template<size_t N>
    struct Factorial {
        static const size_t value = N * Factorial<N - 1>::value;
    };

    template<>
    struct Factorial<0> {
        static const size_t value = 1;
    };

    template<size_t low, size_t high, typename Enabled = void>
    struct Product {
        static const size_t value = high * Product<low, high - 1>::value;
    };

    template<size_t low, size_t high>
    struct Product<low, high, typename std::enable_if<high < low>::type> {
        static const size_t value = 1;
    };

    template<size_t low>
    struct Product<low, low, void> {
        static const size_t value = low;
    };

    template<typename C>
    struct Size {
        static const size_t value = C::SizeAtCompileTime;
    };

    template<typename C>
    struct Degree {
        static const size_t value = C::SizeAtCompileTime - 1;
    };

    /** Return the coefficient of the derivative polynomial at order i.
     * For instance, the coefficient vector, coeffs = [ a, b, c ] represents the polynomial
     *
     * a x^2 + b x + c
     *
     * So for a derivative of 0, the polynomial is unchanged. The coefficient at order
     *
     * i = 0 is c
     * i = 1 is b
     * i = 2 is a
     * i > 2 is 0
     *
     * For a derivative of 1, the polynomial is a 2 a x + b. The coefficient at order
     *
     * i = 0 is b
     * i = 1 is 2a
     * i > 1 is 0
     */
    template<size_t order, size_t derivative, typename C>
    constexpr typename C::Scalar derivativeCoefficient(const C &coeffs) {
        return order > Degree<C>::value - derivative ? 0
                                                     : coeffs(Degree<C>::value - derivative - order) *
                                                       Product<order + 1, order + derivative>::value;
    }

    template<size_t i>
    struct Int {
    };

    template<size_t order, size_t derivative, typename C, typename X>
    constexpr auto
    evalImpl(const C &coeffs, const X &x) -> typename std::enable_if<order >= Degree<C>::value - derivative, X>::type {
        return derivativeCoefficient<order, derivative, C>(coeffs);
    }

    template<size_t order, size_t derivative, typename C, typename X>

    auto
    constexpr
    evalImpl(const C &coeffs, const X &x) -> typename std::enable_if<order<Degree<C>::value - derivative, X>::type {
        return derivativeCoefficient<order, derivative, C>(coeffs) +
               x * evalImpl<order + 1, derivative, C, X>(coeffs, x);
    }
}

namespace Polynomial {

    template<typename C>
    constexpr auto
    toMonic(const C &coeffs, bool is_monic = false) -> Eigen::Matrix<typename C::Scalar, 1, C::SizeAtCompileTime> {
        return is_monic ? coeffs : coeffs / coeffs(0, 0);
    }

    template<typename C>
    constexpr auto roots(const C &coeffs, bool is_monic = false)
    -> decltype(Eigen::Matrix<typename C::Scalar,
            C::SizeAtCompileTime - 1, C::SizeAtCompileTime - 1>::Ones().eigenvalues()) {

        using T = typename C::Scalar;
        const size_t n = C::SizeAtCompileTime - 1;
        Eigen::Matrix<T, n, n> A = Eigen::Matrix<T, n, n>::Zero();
        A.diagonal(1) = Eigen::Matrix<T, 1, n - 1>::Ones();
        A.row(n - 1) = -toMonic(coeffs, is_monic).template rightCols<n>().reverse();
        return A.eigenvalues();
    }

    /** If passed an object with a "coeffs" field (like this class) see if we can compute the roots of those.
     * Note that SNIFAE will ignore this option for things like Eigen matrices since they don't have a "coeffs" field. */
    template<typename P>
    constexpr auto roots(const P &polynomial) -> decltype(roots(polynomial.coeffs)) {
        return roots(polynomial.coeffs, P::is_monic);
    }

    /** Use Horner's method to evaluate the polynomial in a numerically stable way. */
    template<size_t derivative, typename C, typename X>
    constexpr X eval(const C &coeffs, const X &x) {
        return evalImpl<0, derivative, C, X>(coeffs, x);
    }

    template<size_t derivative, typename P, typename X>
    constexpr auto
    eval(const P &polynomial, const X &x) -> decltype(eval<derivative>(polynomial.coeffs, x)) {
        return eval<derivative>(polynomial.coeffs, x);
    }

    template<typename C>
    constexpr auto derivative(const C &coeffs) -> Eigen::Matrix<typename C::Scalar, 1, C::SizeAtCompileTime - 1> {
        return coeffs.template leftCols<C::SizeAtCompileTime - 1>().array() *
               Eigen::Array<typename C::Scalar, 1, C::SizeAtCompileTime - 1>::LinSpaced(1, C::SizeAtCompileTime -
                                                                                           1).reverse();
    }

    template<typename T, size_t i, template<typename, size_t> class P>
    constexpr auto derivative(const P<T, i> &polynomial) -> decltype(P<T, i - 1>(derivative(polynomial.coeffs))) {
        return P<T, i - 1>(derivative(polynomial.coeffs));
    }
};

template<typename T, size_t n>
struct Poly {

    using Scalar = T;
    static const size_t degree = n;

    const Eigen::Matrix<T, 1, n + 1> coeffs;
    const bool is_monic;

    template<typename C>
    Poly(const C &coeffs, bool is_monic = false)
            : coeffs(coeffs), is_monic(is_monic) {
    }

    Poly(const Poly<T, n> &p)
            : Poly(p.coeffs, p.is_monic) {
    }

    Poly<T, n> &operator=(const Poly<T, n> &p) {
        coeffs = p.coeffs;
        is_monic = p.is_monic;
        return *this;
    }

    template<size_t derivative, typename X>
    constexpr X eval(X x) const {
        return Polynomial::eval<derivative>(coeffs, x);
    }

    /** Used for cout */
    friend std::ostream &operator<<(std::ostream &o, const Poly<T, n> &p) {
        o << p.coeffs;
        return o;
    }

    /** Used for cin */
    friend std::istream &operator>>(std::istream &i, const Poly<T, n> &p) {
        i >> p.coeffs;
        return i;
    }
};

template<typename T, size_t n>
struct PolyWithRoots : Poly<T, n> {

    const Eigen::Matrix<std::complex<T>, 1, n> roots;

    template<typename C>
    PolyWithRoots(const C &coeffs, bool is_monic = false)
            : Poly<T, n>(coeffs, is_monic), roots(Polynomial::roots(coeffs, is_monic)) {
    }

    PolyWithRoots(const Poly<T, n> &p)
            : PolyWithRoots(p.coeffs, p.is_monic) {
    }

    /** Used for cout */
    friend std::ostream &operator<<(std::ostream &o, const PolyWithRoots<T, n> &p) {
        o << p.coeffs << std::endl << p.roots;
        return o;
    }
};

template<typename T, size_t n>
struct MonicPoly : Poly<T, n> {

    template<typename C>
    MonicPoly(const C &coeffs)
            : Poly<T, n>(Polynomial::toMonic(coeffs), true) {
    }

    MonicPoly(const Poly<T, n> &p)
            : MonicPoly(p.coeffs, p.is_monic) {
    }
};

template<typename T, size_t n>
struct MonicPolyWithRoots : PolyWithRoots<T, n> {

    template<typename C>
    MonicPolyWithRoots(const C &coeffs)
            : PolyWithRoots<T, n>(Polynomial::toMonic(coeffs), true) {
    }

    MonicPolyWithRoots(const Poly<T, n> &p)
            : MonicPolyWithRoots(p.coeffs, p.is_monic) {
    }
};

#endif /* POLYNOMIAL_HEADER */