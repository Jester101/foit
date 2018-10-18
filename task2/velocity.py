import scipy.constants as const
import numpy as np
import scipy.special as sp


def_power = 10**9
def_eps = 5.7
def_radius = 0.009
def_accuracy = -6
def_start_freq = 2.0 * np.pi
def_max_root = 10


def NewtonMethod(edges, accuracy):
    a = edges[0]
    b = edges[1]
    if sp.jv(0, a)*sp.jvp(0, a, 2) > 0:
        x0 = a
    else:
        x0 = b
    xn = x0 - sp.jv(0, x0)/sp.jvp(0, x0)
    while np.fabs(xn - x0) > accuracy:
        x0 = xn
        xn = x0 - sp.jv(0, x0)/sp.jvp(0, x0)
    return xn


def getVelocity(bessel_arg, freq, radius, eps):
    x = radius**2 * freq**2
    y = bessel_arg**2 * const.c**2
    vel = np.sqrt(np.fabs(x / (eps * x - y)))
    return vel


def VelocitySort(total, max_root):
    velocity = list()
    for k in range(max_root):
        velocity.append([[], []])
    for unit in total:
        fr = unit[0]
        vel = unit[1]
        for i in range(len(vel)):
            velocity[i][0].append(fr)
            velocity[i][1].append(vel[i])
    return velocity


def buildGraph(start_freq=def_start_freq, max_root=def_max_root, eps=def_eps,
               radius=def_radius, accuracy=10**def_accuracy, number=0):
    freq = start_freq * def_power
    print(freq)
    step = True
    total = list()
    votn = np.linspace(1/np.sqrt(eps), 4+1/eps, 1000)
    while step:
        arg = freq/const.c*np.sqrt(np.fabs(eps - 1/votn**2))*radius
        y = sp.jv(0, arg)
        edges = list()
        root = 0
        for l in range(1, len(y)):
            if (y[l]*y[l-1]) < 0:
                root += 1
                if arg[l-1] < accuracy:
                    arg[l-1] = 2
                edges.append([arg[l-1], arg[l]])
        if (root > 0) and (root <= max_root):
            answer = list()
            for edge in edges:
                answer.append(getVelocity(NewtonMethod(edge, accuracy), freq, radius, eps))
            total.append([freq, answer])
        if root > max_root:
            step = False
        freq += 0.1*def_power
    points = VelocitySort(total, max_root)
    resullt = points[number]
    for i in range(len(resullt[0])):
        resullt[0][i] = resullt[0][i]/(2*np.pi*def_power)
    return resullt
