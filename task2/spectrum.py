import numpy


def custom_sin(x):
    if x == 0:
        return 1
    return numpy.sin(x)/x


def create_spectrum(array, impulse_freq, time1_diff):
    spectrum = list()
    for argument in array:
        spec = custom_sin((argument - impulse_freq)*numpy.pi*time1_diff)*time1_diff
        spectrum.append(spec)
    return spectrum
