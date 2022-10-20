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

# mom_2 = athena_read.athdf(files[3])['mom2']
#
# print(mom_2[0][0])
#
# L = []
# t = []
#
# for i in range(len(files)):
#     data = athena_read.athdf(files[i])
#     r = data['x1f'][0:64]
#     mom2 = data['mom2'][0, :, :]
#     L.append(mom2)
#     t.append(data['Time'])
#
# L = np.array(L)
# t = np.array(t)
#
# print(np.sum(mom_2[0][0]))
#
# L_tot = []
# i = 0
# j = 0
# k = 0
# sum = 0
# while i < len(files):
#     while j < 64:
#         while k < 64:
#             sum = sum + L[i][j][k]
#             k = k + 1
#         j = j + 1
#     L_tot.append(sum)
#     i = i + 1
#
# L_tot = np.array(L_tot)

# print(L_tot)

# plt.plot(t,L_tot)

# data = athena_read.athdf(files[74])
# mom3 = data['mom3'][0, :, :]
#
# print(np.max(mom3))
# print(np.min(mom3))

L = []
t = []

for k in range(len(files)):
    data = athena_read.athdf(files[k])
    r = data['x1f'][0:64]
    theta = data['x2f'][0:64]
    t.append(data['Time'])
    sum = 0
    i = 0
    j = 0
    while i <= 63:
        while j <= 63:
            sum = sum + r[i] * data['mom2'][0, j, i]
            j = j + 1
        i = i + 1
    L.append(sum)

print(len(L))
print(len(t))

plt.plot(t,L)


plt.show()