if False:
    from autograd import numpy as np
else:
    import numpy as np
from autograd import jacobian
import lm
import sys

import scipy

import numdifftools


def ALPHA():
    return 100


def GAMMA():
    return 0.99


def LOW_SCALE():
    return 1 / GAMMA()


def HIGH_SCALE():
    return (ALPHA() - 1) / (1 - GAMMA())


def HIGH_INTERCEPT():
    return (1 - GAMMA() * ALPHA()) / (1 - GAMMA())


def LANES():
    return 3


def LANE_WIDTH():
    return 4


def ROAD_WIDTH():
    return LANES() * LANE_WIDTH()


def SAFE_S_DISTANCE():
    return 20


def S_BAR():
    return SAFE_S_DISTANCE() / np.sqrt(2 - 2 * GAMMA())


def SAFE_D_DISTANCE():
    return 2


def D_BAR():
    return SAFE_D_DISTANCE() / np.sqrt(2 - 2 * GAMMA())


def SIGMA():
    return 10


def REFERENCE_VELOCITY():
    return 45


# The prototypical cost function that we use for the critical constraints.
def MAX_POWER():
    return 1.5


def proto(z):
    if False:
        power = (z - GAMMA()) / (1 - GAMMA())
        return np.power(ALPHA(), min(power, MAX_POWER()))
    else:
        if z < 0:
            return 0
        elif z < GAMMA():
            return z * LOW_SCALE()
        else:
            return z * HIGH_SCALE() + HIGH_INTERCEPT()


# A constraint meant to ensure that we stay on the correct side of the road
def road(s, d):
    return proto(d / ROAD_WIDTH()) + proto(1 - d / ROAD_WIDTH())


# Critical constraints for the velocity, acceleration, and jerk.
def velocity(s_dot, d_dot):
    return proto((np.square(s_dot) + np.square(d_dot)) / 500)


def acceleration(s_ddot, d_ddot):
    return proto((np.square(s_ddot) + np.square(d_ddot)) / 100)


def jerk(s_dddot, d_dddot):
    return proto((np.square(s_dddot) + np.square(d_dddot)) / 100)


# Collision avoidance cost
def collisionAvoidance(s, d, s_nearby, d_nearby):
    cost = 0
    for (s_i, d_i) in zip(s_nearby, d_nearby):
        cost += proto(1 - np.square((s_i - s) / S_BAR())) \
                * proto(1 - np.square((d_i - d) / D_BAR()))
    return cost


# Lane-keeping cost
def lane(s, d):
    return np.square(np.cos(np.pi * d / LANE_WIDTH()))


# Reference velocity cost
def referenceVelocity(s_dot, d_dot):
    return 1 - np.exp(-np.square((s_dot - REFERENCE_VELOCITY()) / SIGMA()))


def polyder(c):
    return c[:-1] * range(c.size - 1, 0, -1)


def polyval(c, t):
    return np.dot(c, (t ** np.arange(c.size - 1, -1, -1)))


class S:
    n_params = 2

    @staticmethod
    def eval(params, t, deriv):
        ''' Evaluate the specified derivative of S at t using the specified parameters '''
        for i in range(deriv):
            params = polyder(params)
        # return polyval(params, t)
        return 0


class D:
    n_params = 3

    @staticmethod
    def eval(params, t, deriv):
        ''' Evaluate the specified derivative of D at t using the specified parameters '''
        d0, v0, a0 = 0, 0, 0
        if params.size == 2:
            x, y = params
            delta = 0
        else:
            x, y, d = params
            delta = -(a0 + v0 * x + v0 * y - d * x * y + d0 * x * y) / ((x - 1) * (y - 1))
        alpha = d0
        beta = -(a0 + delta - delta * y + v0 * y) / (x * (x - y))
        gamma = (a0 + delta - delta * x + v0 * x) / (y * (x - y))
        if deriv == 0:
            return alpha - delta * (np.exp(-t) - 1) - beta * (np.exp(-t * x) - 1) - gamma * (np.exp(-t * y) - 1)
        else:
            return -(-1) ** deriv * (
                    delta * np.exp(-t) + beta * x ** deriv * np.exp(-t * x) + gamma * y ** deriv * np.exp(-t * y))


# The total cost function
def totalCost(times, s_nearby, d_nearby, x):
    d_params = x[:D.n_params]
    s_params = x[D.n_params:]
    f = ()
    for i in range(times.size):
        # Evaluate s and d (as well as their derivatives) at t
        t = times[i]
        s = S.eval(s_params, t, 0)
        d = D.eval(d_params, t, 0)
        s_dot = S.eval(s_params, t, 1)
        d_dot = D.eval(d_params, t, 1)
        s_ddot = S.eval(s_params, t, 2)
        d_ddot = D.eval(d_params, t, 2)
        s_dddot = S.eval(s_params, t, 3)
        d_dddot = D.eval(d_params, t, 3)

        f = np.hstack((f, (road(s, d),
                           velocity(s_dot, d_dot),
                           acceleration(s_ddot, d_ddot),
                           jerk(s_dddot, d_dddot),
                           collisionAvoidance(s, d, s_nearby, d_nearby),
                           lane(s, d),
                           referenceVelocity(s_dot, d_dot))))
    # print("x: ", x.transpose())
    # print("sum(f**2): ", np.square(f).sum())
    return f


# Create the list of times at which to evaluate the cost
time_horizon = 3
N = 10
times = np.linspace(0, time_horizon, N)

# Create the initial nearby cars
s_nearby = []
d_nearby = []

# Create the initial guess
d_params = np.array([2, 3, 4])
degree = 6
s_params = np.zeros((degree + 1), )
x0 = np.hstack((d_params, s_params))

# Define the problem
f = lambda x: totalCost(times, s_nearby, d_nearby, x)
# df = jacobian(f)
df = numdifftools.core.Jacobian(f)

print("f(x0)", f(x0).sum())
for i in range(2, 6):
    print("f(x0+1e-" + str(i) + ")", f(x0 + 10 ** -i).sum(), (f(x0).sum() - f(x0 + 10 ** -i).sum()))

print()

# The stop condition for the optimization routine
maxIts = 1000
dcostTol = 1e-5
dxTol = 1e-5


def stop(x, fx, dfdx, cost, it, step, dcost, dx, scales, updated):
    if dcost is None or dx is None or np.isnan(dcost) or np.isinf(dcost):
        done = False
    else:
        done = it >= maxIts
        if not done and updated:
            done = np.abs(dcost) < dcostTol or np.linalg.norm(dx) < dxTol
    return done


# Optimize
tryLM = True
if tryLM:
    x, fx, dfdx, cost, it, step, dcost, dx, scales, updated = lm.LM(np).solve(f, df, x0, stop=stop, scaling=False)
    print("it", it)
    print("x", x)
    print("cost", cost)
    print("dcost", dcost)
    print("np.linalg.norm(dx)", np.linalg.norm(dx))
else:
    F = lambda x: f(x).sum()
    sol = scipy.optimize.minimize(F, x0, tol=1e-20)
    print(sol)
