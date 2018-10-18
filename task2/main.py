import scipy.fftpack as spfft
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plot
from velocity import buildGraph
from impulse import create_impulse
from spectrum import create_spectrum
from signal import signal_function_dist, signal_function_time

POWER = 10**9
SIGNAL_FREQ = 13
SIZE = 100
STEP = 0.01
TIME_1 = 0.7692
TIME_2 = 1.538
TIME_3 = 0.07692
LIMIT = 0.5
N_COUNT = 2**10
N1_COUNT = 60
NUMBER = 1000

def signal_transmission():
    basic_velocity = buildGraph(max_root=2)
    border_freq = np.min(basic_velocity[0])
    velocity_function = interp1d(basic_velocity[0], basic_velocity[1])
    freq_array = list()
    freq_array.append(border_freq)
    for i in range(1, SIZE):
        freq_array.append(freq_array[i-1]*1.01)
    last_freq = np.max(basic_velocity[0])
    print(border_freq)
    time_array = [TIME_2/(N_COUNT-1)*i for i in range(N_COUNT)]
    impulse_array = create_impulse(time_array, TIME_1, SIGNAL_FREQ)
    fast_fourier = spfft.fft(impulse_array)
    fast_fourier_print = fast_fourier/np.max(fast_fourier)
    k = [i/TIME_2 for i in range(N1_COUNT+1)]
    spectrum = create_spectrum(k, SIGNAL_FREQ, TIME_1)
    spectrum2 = np.fabs(spectrum)
    spectrum2 = spectrum2/np.max(spectrum2)
    distance = 0
    time = [4*TIME_2/(NUMBER-1) * i for i in range(NUMBER)]
    time_modify = time
    for i in range(len(time_modify)):
        time_modify[i] = time_modify[i]-TIME_1/2
    signa_time = signal_function_time(time_modify, distance, border_freq, last_freq, k, spectrum, velocity_function)

    distance = [i for i in range(SIZE)]
    signal_dist = signal_function_dist(0+TIME_1/2, distance, border_freq, last_freq, k, spectrum, velocity_function)


    fig = plot.figure()
    plot.subplot(2, 3, 1)
    plot.plot(basic_velocity[0], basic_velocity[1], 'blue')
    plot.grid(True, color='w')
    plot.title("Повторимый оригинал")

    plot.subplot(2, 3, 4)
    plot.plot(freq_array, velocity_function(freq_array), 'red')
    plot.grid(True, color='w')
    plot.title("Достойная пародия")

    plot.subplot(2, 3, 2)
    plot.plot(time_array, impulse_array, 'green')
    plot.xlim(0, TIME_2)
    plot.grid(True, color='w')
    plot.title("Импульс")

    plot.subplot(2, 3, 5)
    plot.plot(fast_fourier_print, 'red')
    plot.plot(spectrum2, 'blue')
    plot.xlim(0, np.max(k))
    plot.grid(True, color='w')
    plot.title("Спектр")

    plot.subplot(2, 3, 3)
    plot.plot(time, signa_time, 'red')
    plot.grid(True, color='w')
    plot.title("Сигнал")

    plot.subplot(2, 3, 6)
    plot.plot(distance, signal_dist, 'blue')
    plot.grid(True, color='w')
    plot.title("Сигнал")

    print("end")
    plot.show()


if __name__ == "__main__":
    signal_transmission()