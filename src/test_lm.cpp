#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"
#include "polynomial.h"
#include "cost.h"

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

int main(int argc, char *argv[]) {

    using T = double;

    /* Create the vector of times at which we will evaluate the cost */
    const size_t time_horizon = 3;
    const size_t N = 50;
    using Times = Matrix<T, N, 1>;
    const Times times = Times::LinSpaced(time_horizon / N, time_horizon);

    /* Create the initial nearby cars */
    using NEARBY = Matrix<T, Dynamic, 1>;
    NEARBY s_nearby = NEARBY::Zero(0);
    NEARBY d_nearby = NEARBY::Zero(0);

    /* Define the S and D trajectory parameterizations */
    S<T> s;
    D<T> d;

    /* Create the initial guess */
    const size_t degree = 3;
    Matrix<T, Dynamic, 1> x = Matrix<T, Dynamic, 1>::Zero(decltype(s)::N_PARAMS + decltype(d)::N_PARAMS,
                                                          1);
    s.goodInitialGuess(x.data());
    d.goodInitialGuess(x.data() + s.N_PARAMS);

    /* Define the problem */
    using Cost = TotalCost<T, Times, NEARBY, NEARBY>;
    Cost total_cost(times, s_nearby, d_nearby, s, d);
    Eigen::NumericalDiff<Cost> numDiff(total_cost);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Cost>, T> lm(numDiff);
    lm.parameters.maxfev = 100;
    lm.parameters.xtol = 1.0e-5;
    lm.parameters.ftol = 1.0e-5;

    /* Optimize */
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    int ret;
    for (int i = 0; i < 1; ++i) {
        ret = lm.minimize(x);
    }

    /* Display the results */
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    cout << "Time (s): " << duration / 1e6 << endl;
    cout << "Iterations: " << lm.iter << endl;
    cout << "Return code: " << ret << endl;
    cout << "F: " << lm.fvec.transpose() << endl;
//    cout << "dF: " << lm.fjac << endl;

    cout << "x that minimizes the function: " << endl << x.transpose() << endl;
    return 0;
}