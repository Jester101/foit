import numpy as np
import scipy.constants as spc

POWER = 10**9


def signal_function_time(time_array, z, board_freq, last_freq, arguments, spectrum, velocity_function):
    signal = list()
    for time in time_array:
        sum = 0
        for k in range(len(arguments)):
            if board_freq <= arguments[k] <= last_freq:
                sum += spectrum[k] * np.cos(arguments[k]*2*np.pi *
                        (time - z*POWER /(velocity_function(arguments[k])*spc.c)))
        signal.append(sum)
    return signal


def signal_function_dist(time, dist_array, board_freq, last_freq, arguments, spectrum, velocity_function):
    signal = list()
    for dist in dist_array:
        sum = 0
        for k in range(len(arguments)):
            if board_freq <= arguments[k] <= last_freq:
                sum += spectrum[k] * np.cos(arguments[k]*2*np.pi * time -
                    (arguments[k]*2*np.pi*dist*POWER /(velocity_function(arguments[k])*spc.c)))
        signal.append(sum)
    return signal
