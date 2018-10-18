import numpy


def stage_function(argument):
    if argument < 0:
        return 0
    elif argument > 0:
        return 1
    else:
        return 0.5


def create_impulse(time_array, time_diff, impulse_freq):
    impulse = list()
    for time in time_array:
        impulse.append(numpy.sin(impulse_freq*2*numpy.pi*time)*
                       (stage_function(time-time_diff)-stage_function(time)))
    return impulse
