from matplotlib import rcParams
import matplotlib.pyplot as plt
import scipy.constants as const
import numpy as np
import scipy.special as sp

rcParams['font.family'] = 'Times New Roman', 'Arial', 'Tahoma'
rcParams['font.fantasy'] = 'Times New Roman'

# Изменение параметров рисования (смена чёрного по белому на белое по чёрному)
# facecolor = 'k'
#
# rcParams['figure.edgecolor'] = facecolor
# rcParams['figure.facecolor'] = facecolor
# rcParams['axes.facecolor'] = facecolor
# #rcParams['axes.edgecolor'] = 'k'
# rcParams['grid.color'] = 'w'
# rcParams['xtick.color'] = 'w'
# rcParams['ytick.color'] = 'w'
# rcParams['axes.labelcolor'] = 'w'


eps = 5.7
radio = 0.009
power = 10**9
accuracy = 10**-6
start_freq = 2.0 * np.pi

def NewtonMethod(edges):
    a = edges[0]
    b = edges[1]
    if sp.jv(0,a)*sp.jvp(0,a,2) > 0:
        x0 = a
    else:
        x0 = b
    xn = x0 - sp.jv(0,x0)/sp.jvp(0,x0)
    while np.fabs(xn - x0) > accuracy:
        x0 = xn
        xn = x0 - sp.jv(0,x0)/sp.jvp(0,x0)
    return xn

def getVelocity(bassel_arg, freq):
    x = radio**2 * freq**2
    y = bassel_arg**2 * const.c**2
    vel = np.sqrt(x /(eps* x - y))
    return vel

def VelocitySort(total):
    velocity = list()
    for k in range(10):
        velocity.append([[],[]])
    for unit in total:
        fr = unit[0]
        vel = unit[1]
        for i in range(len(vel)):
            velocity[i][0].append(fr)
            velocity[i][1].append(vel[i])
    return velocity

def buildGraph():
    freq = np.arange(2*np.pi*power,70*2*np.pi*power,0.1*power)
    max_root = 0
    total = list()
    votn = np.linspace(1/np.sqrt(eps),4+1/eps,1000)
    for i in range (len(freq)):
        arg = freq[i]/const.c*np.sqrt(eps - 1/votn**2)*radio
        y = sp.jv(0,arg)
        edges = list()
        root = 0
        for l in range(1,len(y)):
            if (y[l]*y[l-1]) < 0:
                root +=1
                # если нижняя граница меньше 2, то Ньютон выдаёт отрицательные числа
                if arg[l-1] < accuracy:
                    arg[l-1] = 2
                edges.append([arg[l-1],arg[l]])
                # print('votn: %(1)f, arg: %(3)f y: %(2)f' %{'1':votn[i],'2':y[i],'3':arg[i]})
        # print(edges)
        if root > 0:
            answer = list()
            if root > max_root:
                print(root)
                max_root = root
            ar = list()
            for edge in edges:
                # ar.append(NewtonMethod(edge))
                answer.append(getVelocity(NewtonMethod(edge),freq[i]))
            total.append([freq[i],answer])
            # print(sp.jv(0,ar))
            # print(answer)
            # if i == len(freq)-1:
            # for k in answer:
            #     # print('arg: ',k)
            #     arr = getVelocity(k,freq[i])
            #     # print('vel: ',arr)
            #     # print(freq[i]/const.c*np.sqrt(eps - 1/arr**2)*radio)
            #     print(sp.jv(0,freq[i]/const.c*np.sqrt(eps - 1/arr**2)*radio))
            #     # print(sp.jv(0,k))
            # # print('x: ', answer)
            # print('y: ',sp.jv(0,answer))\
    points = VelocitySort(total)
    fig = plt.figure()
    for i in range(10):
        plt.plot(points[i][0],points[i][1],'blue')
        # plt.xlim(np.min(points[i][0]),np.max(points[i][0]))
        # plt.ylim(np.min(points[i][1]),np.max(points[i][1]))
        plt.grid(True, color='w')
    plt.xlim(0,np.max(points[i][0]))
    plt.ylim(0,np.max(points[i][1]))
    plt.show()

if __name__ == "__main__":
    buildGraph()
