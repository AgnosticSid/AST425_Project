import sys
sys.path.insert(0, '/home/ubuntu/athena/vis/python')
import athena_read
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from glob import glob
import natsort
from numpy import diff
from scipy.optimize import curve_fit

files = glob("/home/ubuntu/athena/bin/disk.out1.*.athdf")
files = natsort.natsorted(files)

time,mass,mom_1,mom_2,mom_3,KE_1,KE_2,KE_3,tot_E = np.loadtxt('hst.txt',dtype=str,skiprows=2,unpack=True,usecols=(0,2,3,4,5,6,7,8,9))

radius,density = np.loadtxt('v_theta.txt',dtype=float,unpack=True,delimiter=',')

time = time.astype('float64')
mass = mass.astype('float64')
mom_1 = mom_1.astype('float64')
mom_2 = mom_2.astype('float64')
mom_3 = mom_3.astype('float64')
KE_1 = KE_1.astype('float64')
KE_2 = KE_2.astype('float64')
KE_3 = KE_3.astype('float64')
tot_E = tot_E.astype('float64')

KE_tot = KE_1+KE_2+KE_3
mom_tot = np.sqrt(mom_1**(2) + mom_2**(2) + mom_3**(2))

# print(np.max(density))
# print(np.min(density))

# plt.scatter(time[750:],tot_E[750:],s=5)
# plt.xlabel('Time')
# plt.ylabel('Total Energy')
# plt.axhline(y=9.787415,color='red')

plt.scatter(time,mass,s=5)
plt.xlabel('Time')
plt.ylabel('Total Mass')
plt.show()

plt.scatter(time,tot_E,s=5)
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.show()

plt.scatter(time,KE_tot,s=5)
plt.xlabel('Time')
plt.ylabel('Total KE')
plt.show()

print(mass[-1])
print(KE_tot[-1])

# plt.axhline(y=9.787415,color='red')

# plt.scatter(time[750:],mass[750:],s=5)
# plt.xlabel('Time')
# plt.ylabel('Mass')
# plt.axhline(y=20.9891,color='red')

# plt.scatter(time,tot_E,s=5)
# plt.xlabel('Time')
# plt.ylabel('Mass')
# plt.axhline(y=20.9891,color='red')

# print(athena_read.athdf(files[0])['x1f'])


# avg_E = (np.max(tot_E[750:])+np.min(tot_E[750:]))/2
#
# print(avg_E)
#
# print(np.abs(((np.max(tot_E[750:])-avg_E)/avg_E)*100))
# print(np.abs(((np.min(tot_E[750:])-avg_E)/avg_E)*100))

# avg_m = (np.max(mass[750:])+np.min(mass[750:]))/2
#
# print(avg_m)
#
# print(np.abs(((np.max(mass[750:])-avg_m)/avg_m)*100))
# print(np.abs(((np.min(mass[750:])-avg_m)/avg_m)*100))

