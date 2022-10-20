import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Freq_1,Phase_1=np.loadtxt('Ph_34.txt',unpack=True)

Phase_1=np.unwrap(Phase_1*(np.pi/180))-(6*2*np.pi)
Freq_1 = Freq_1/10**9

def f(freq,m,k):
    return m*freq+k

popt,pcov=curve_fit(f,Freq_1,Phase_1)

m = popt[0]
k = popt[1]

print('m=',m,'k=',k)

Time_delay = Phase_1/(Freq_1*2*np.pi)

plt.plot(Freq_1,Time_delay)
plt.xlabel('Frequency [GHz]')
plt.ylabel('Time Delay [ns]')
plt.title('Nominal Power = -34 dBm')

slope = (Phase_1[-1]-Phase_1[0])/(Freq_1[-1]-Freq_1[0]) # Radians/GHz


plt.show()
