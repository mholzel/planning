import matplotlib.pyplot as plt
import numpy

# Specify the time vector
tn = numpy.linspace(-5, 0, 100)
tp = numpy.linspace(0, +5, 100)

# The method that will handle plotting the specified derivatives of the given polynomials
def plot(p, s, derivative, filename):
    # Compute the derivatives of the polynomials
    dp = numpy.polyder(p, m=derivative)
    ds = numpy.polyder(s, m=derivative)

    fig = plt.figure()
    fig.set_size_inches(5, 2)
    if derivative > 0:
        pre = "\\"
        pre += "d" * derivative
        pre += "ot{"
        post = "}"
    else:
        pre = ""
        post = ""
    plabel = "$\,\,\, \\underbar \!\!\!" + pre + "p" + post + "(t)$"
    slabel = "$" + pre + "p" + post + "(t)$"
    print("Labels:", plabel, slabel)
    ph, = plt.plot(tn, numpy.polyval(dp, tn), label=plabel)
    sh, = plt.plot(tp, numpy.polyval(ds, tp), label=slabel)
    plt.legend(handles=[ph, sh])
    plt.grid()
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig(filename)


# Discontinuous in the position
p = [0]
s = [1]
plot(p, s, 0, "discontinuous_position.svg")

# Discontinuous in the velocity
s = [1, 0]
plot(p, s, 0, "discontinuous_velocity_position.svg")
plot(p, s, 1, "discontinuous_velocity_velocity.svg")

# Discontinuous in the acceleration
s = [1, 0, 0]
plot(p, s, 0, "discontinuous_acceleration_position.svg")
plot(p, s, 1, "discontinuous_acceleration_velocity.svg")
plot(p, s, 2, "discontinuous_acceleration_acceleration.svg")
