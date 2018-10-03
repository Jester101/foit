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


def_eps = 5.7
def_radius = 0.009
def_power = 10**9
def_accuracy = -6
def_start_freq = 2.0 * np.pi
def_max_root = 10

def NewtonMethod(edges,accuracy):
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

def getVelocity(bassel_arg, freq, radius, eps):
    x = radius**2 * freq**2
    y = bassel_arg**2 * const.c**2
    vel = np.sqrt(np.fabs(x /(eps* x - y)))
    return vel

def VelocitySort(total, max_root):
    velocity = list()
    for k in range(max_root):
        velocity.append([[],[]])
    for unit in total:
        fr = unit[0]
        vel = unit[1]
        for i in range(len(vel)):
            velocity[i][0].append(fr)
            velocity[i][1].append(vel[i])
    return velocity

def buildGraph(start_freq, max_root, eps, radius, accuracy):
    # freq = np.arange(2*np.pi*power,70*2*np.pi*power,0.1*power)
    print("Подождите немного, идёт обработка данных...")
    freq = start_freq * def_power
    step = True
    total = list()
    votn = np.linspace(1/np.sqrt(eps),4+1/eps,1000)
    # print(votn)
    # print(np.sqrt(eps - 1/votn**2))
    while step:
        # print(freq)
        try:
            arg = freq/const.c*np.sqrt(np.fabs(eps - 1/votn**2))*radius
        except Exception as e:
            print(arg)
            print(np.sqrt(eps - 1/votn**2))
        y = sp.jv(0,arg)
        edges = list()
        root = 0
        for l in range(1,len(y)):
            if (y[l]*y[l-1]) < 0:
                root +=1
                # Менять ничего тут нельзя. вообще
                if arg[l-1] < accuracy:
                    arg[l-1] = 2
                edges.append([arg[l-1],arg[l]])
        if (root > 0) and (root <= max_root):
            answer = list()
            for edge in edges:
                answer.append(getVelocity(NewtonMethod(edge,accuracy),freq,radius,eps))
            total.append([freq,answer])
        if root > max_root:
            step = False
        freq += 0.1*def_power
    points = VelocitySort(total,max_root)
    group = list()
    for i in range(len(points)):
        group.append([[],[]])
        for k in range(len(points[i][0])-1):
            group[i][1].append((points[i][0][k+1]-points[i][0][k])/
            (points[i][0][k+1]/points[i][1][k+1] - points[i][0][k]/points[i][1][k]))
            group[i][0].append(points[i][0][k])
        # group_velocity.append(group_velocity[k])
        # print(group_velocity)
        # group.append(group_velocity)

    # print(group)
    # print(points[1])
    fig = plt.figure()
    for i in range(max_root):
        # print(i)
        plt.subplot(1,3,3)
        plt.plot(points[i][0],points[i][1],'blue')
        plt.grid(True, color='w')
    for i in range(max_root):
            plt.subplot(1,3,3)
            plt.plot(group[i][0],group[i][1],'red')
            # plt.plot(group[i][0],group[i][1],'red')
            plt.grid(True, color='w')
    plt.yscale('logit')
    plt.xlabel(u'Частота')
    plt.ylabel(u'Относительная фазовая скорость')
    for i in range(max_root):
        # print(i)
        plt.subplot(1,3,1)
        plt.plot(points[i][0],points[i][1],'blue')
        plt.grid(True, color='w')
    plt.xlabel(u'Частота')
    plt.ylabel(u'Относительная фазовая скорость')
    for i in range(max_root):
            plt.subplot(1,3,2)
            plt.plot(group[i][0],group[i][1])
            plt.grid(True, color='w')
    plt.xlabel(u'Частота')
    plt.ylabel(u'Относительная фазовая скорость')
    plt.show()


def initData():
    print("Введите Eps или '0' для значения по умолчанию \n" +
    "значения по умолчанию = %f: " %(def_eps) )
    eps = float(input()) or def_eps
    print("Введите Радиус в метрах или '0' для значения по умолчанию \n" +
    "значения по умолчанию = %f: " %(def_radius))
    radius = float(input()) or def_radius
    print("Введите порядок точности вычислений или '0' для значения по умолчанию \n" +
    "значения по умолчанию = %f: " %(10**def_accuracy))
    accuracy = int(input()) or def_accuracy
    accuracy = 10**accuracy
    max_root = def_max_root
    print("Введите начальную частоту в степени 2*pi*10^9 или '0' для значения по умолчанию \n" +
    "значения по умолчанию = %f * %e: " %(def_start_freq, def_power))
    start_freq = float(input()) or def_start_freq
    print("Ваши данные:")
    print ("Eps = %f \n"
            "Радиус = %f \n"
            "Точность вычисления = %f \n"
            "Количество кривых = %f \n"
            "Начальная частота = %f * %e \n"
            %(eps, radius, accuracy, max_root, start_freq, def_power))
    print("Нажмите Enter, чтобы начать")
    input()
    buildGraph(start_freq, max_root,eps,radius,accuracy)
if __name__ == "__main__":
    initData()
