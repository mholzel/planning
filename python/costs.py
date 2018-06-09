import matplotlib.pyplot as plt
from matplotlib import cm
import numpy
import sys
from mpl_toolkits.mplot3d.axes3d import Axes3D

# Protypical cost function 
proto = lambda alpha,gamma,z: alpha**( ( z - gamma ) / ( 1 - gamma ) )

if True:

    # Protoype
    def plotCost( params, z, f, xlabel, filename ):
        fig = plt.figure()
        fig.set_size_inches(5, 4.5)
        handles = []
        for (alpha,gamma) in params:
            label = "$(\\alpha,\gamma) = (" + str(alpha) + "," + str(gamma) + ")$"
            h, = plt.plot( z, f(alpha,gamma,z), label = label)
            handles.append(h)
        plt.grid()
        plt.legend(handles=handles)
        plt.xlabel(xlabel)
        plt.ylabel("cost")
        plt.ylim((0,10))
        plt.savefig(filename)

    params = ((10,.9),(10,.95),(10,.99))
    params2 = ((10,.5),(10,.7),(10,.9))

    # Staying on the correct side of the road 
    f = lambda alpha,gamma,d: proto(alpha,gamma,d/12) + proto(alpha,gamma,1-d/12)

    plotCost( params, numpy.linspace( 0, 12, 1000 ), f, "$d$ (m)", "correct_lane.svg" )
    plotCost( params, numpy.linspace( 0, 2, 1000 ), f, "$d$ (m)", "correct_lane_left.svg" )
    plotCost( params, numpy.linspace( 10, 12, 1000 ), f, "$d$ (m)", "correct_lane_right.svg" )
    
    # Velocity, acceleration, and jerk
    f = lambda alpha,gamma,x: proto(alpha,gamma,x**2 / 500)
    plotCost( params, numpy.linspace( 20, numpy.sqrt(500), 1000 ), f, "$velocity$ (m/s)", "velocity.svg" )

    f = lambda alpha,gamma,x: proto(alpha,gamma,x**2 / 100)
    plotCost( params2, numpy.linspace( 0, 10, 1000 ), f, "$acceleration$ (m/s^2)", "acceleration.svg" )

    f = lambda alpha,gamma,x: proto(alpha,gamma,x**2 / 100)
    plotCost( params2, numpy.linspace( 0, 10, 1000 ), f, "$jerk$ (m/s^3)", "jerk.svg" )

# Collision avoidance
gamma  = 0.95
s_i    = 0
d_i    = 0

s_star = 20
s_bar  = s_star / numpy.sqrt(1-gamma) / numpy.sqrt(2)

d_star = 2
d_bar  = d_star / numpy.sqrt(1-gamma) / numpy.sqrt(2)

f      = lambda alpha,gamma,s,d: proto(alpha,gamma,1 - ((s-s_i) / s_bar)**2) * proto(alpha,gamma,1 - ((d-d_i) / d_bar)**2)

# First, evaluate at an alpha of 10
s      = numpy.linspace( -20, 20, 35 )
d      = numpy.linspace(  -2,  2, 35 )
S,D    = numpy.meshgrid( s, d )
G      = f(10,gamma,S,D)

fig    = plt.figure()
fig.set_size_inches(6,6)
fig.subplots_adjust(left=0, bottom=0, right=.925, top=1, wspace=0, hspace=0)
rect = fig.add_subplot(1, 1, 1).get_position()
plt.axis('off')
ax     = Axes3D(fig, rect)
surf   = ax.plot_surface(S, D, G, cmap=plt.get_cmap('plasma'))
plt.grid()
ax.set_zlim( (0,100) )
ax.set_xlabel('$s - s_i$ (m)')
ax.set_ylabel('$d - d_i$ (m)')
ax.set_zlabel('cost')
ax.set_yticks(list(range(-2,3, 1)))
ax.set_xticks(list(range(-20,21,10)))
ax.view_init(15,-45)
plt.savefig("collision_3d.svg")

s      = numpy.linspace( -30, 30, 1000 )
d      = numpy.linspace(  -3,  3, 1000 )
S,D    = numpy.meshgrid( s, d )
G      = f(10,gamma,S,D)

levels = [ 0.01, 1, 10, 40 ]
fmt = {}
for lev in levels:
    fmt[lev] = str(lev)
fig    = plt.figure()
fig.set_size_inches(5, 5)
CS = plt.contour(D,S,G,levels)
plt.clabel(CS, inline=1, fmt=fmt)
plt.ylabel('$s - s_i$ (m)')
plt.xlabel('$d - d_i$ (m)')
plt.grid()
plt.savefig("collision_contours.svg")
     

sys.exit()
    
# Speed limit cost
def speedLimitCost(params):
    filename = "speed_limit.svg"
    x = numpy.linspace( 15, 22.352, 1000 )
    fig = plt.figure()
    fig.set_size_inches(5, 4.5)
    handles = []
    for (a,xstar) in params:
        label = "$(\\alpha_{v},v^\star) = (" + str(a) + "," + str(xstar) + ")$"
        h, = plt.plot( x, numpy.exp( a * (x**2 - xstar**2) ), label = label)
        handles.append(h)
    plt.grid()
    plt.legend(handles=handles)
    plt.xlabel("$\sqrt{\\dot{s} + \\dot{d}}$")
    plt.ylabel("cost")
    plt.ylim((0,10))
    plt.savefig(filename)

speedLimitCost( ( (.01,19), (.03,20), (.05,21), ) )

# Lane positions 
def laneCost(powers):
    filename = "lane_cost.svg"
    d = numpy.linspace( 0, 12, 100 )
    fig = plt.figure()
    fig.set_size_inches(5, 4)
    handles = []
    for power in powers:
        label = "$\\beta = " + str(power) + "$"
        h, = plt.plot(d, numpy.cos( numpy.pi * d / 4 )**power, label = label)
        handles.append(h)
    plt.grid()
    plt.legend(handles=handles)
    plt.xlabel("d")
    plt.ylabel("cost")
    plt.savefig(filename)

laneCost([2,6,20])

# Reference velocity cost
def referenceVelocityCost(v, stds):
    filename = "reference_velocity_costs.svg"
    max_std = numpy.max(stds) 
    s_dot = numpy.linspace( 0, v+3*max_std, 1000 )
    fig = plt.figure()
    fig.set_size_inches(5, 4)
    handles = []
    for std in stds:
        h, = plt.plot( s_dot, 1 - numpy.exp( -( s_dot - v )**2 / std**2 ), label = "$\\sigma = " + str(std) + "$")
        handles.append(h)
    plt.grid()
    plt.legend(handles=handles)
    plt.xlabel("$\dot{s}$")
    plt.ylabel("cost")
    plt.savefig(filename)

referenceVelocityCost( 45, [1, 5, 10] )