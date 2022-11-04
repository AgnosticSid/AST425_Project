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

# data = athena_read.athdf(files[0])

# #
# fig, ax = plt.subplots(1,subplot_kw=dict(projection='polar'))
# line, = plt.plot([], [], lw=2)
#
# x1 = data['x1v']
# y1 = data['x2v']
# z3 = np.transpose(data['vel2'][0])
# plt.contourf(y1, x1, z3)

# r = []
# t = []
# j = 0
# for i in range(len(files)):
#     t.append(data['Time'])
#     data = athena_read.athdf(files[i])
#     while j < 64:
#         r_ind = np.min(data['x1f'][np.where(data['rho'][0,j,:]==np.max(data['rho'][0,j,:]))])
#         r.append(r_ind)
#         j = j + 1
#
# r = np.array(r)
# t = np.array(t)
#
# print(np.shape(r))

#
# w = w[0:39]
# t = t[0:39]
#
# dwdt = diff(w)/diff(t)

# print(w)
#
# plt.plot(t[1:],dwdt)

# r = data['x1f'][0:64]
# theta = data['x2f'][0:64]
# rad = np.tile(r,(64,1))
# rho = data['rho'][0, :, :]
# v_theta = data['vel2'][0, :, :]
#
# L = rho*rad*v_theta
#
# z3 = np.transpose(L)
# plt.contourf(theta, r, z3)
# fig.colorbar(plt.contourf(theta, r, z3))
#

# def f(arr):
#     L_tot0 = []
#     for i in range(arr):
#         for j in range(arr):
#             sum = 0
#             sum = sum+L[i][j]
#     L_tot0.append(sum)
#     return L_tot0

# print(f(64))

L = []
t = []
M = []
E = []

for k in range(len(files)):
    data = athena_read.athdf(files[k])
    r = data['x1f'][0:64]
    theta = data['x2f'][0:64]
    t.append(data['Time'])
    L_sum = 0
    E_sum = 0
    M_sum = 0
    i = 0
    j = 0
    for i in range(63):
        for j in range(63):
            L_sum = L_sum + ((r[i+1]+r[i])/2)**(2) * data['rho'][0, j, i] * data['vel2'][0, j, i] * (r[i+1]-r[i]) * (theta[j+1]-theta[j])
            M_sum = M_sum + data['rho'][0, j, i] * (r[i+1]-r[i]) * (theta[j+1]-theta[j]) * ((r[i+1]+r[i])/2)
            E_sum = E_sum + (1/2)*(data['rho'][0, j, i]*(r[i+1]-r[i])*(theta[j+1]-theta[j])*((r[i+1]+r[i])/2))*((data['vel1'][0, j, i]**(2)+data['vel2'][0, j, i]**(2)+data['vel3'][0, j, i]**(2)))
    L.append(L_sum)
    M.append(M_sum)
    E.append(E_sum)

plt.scatter(t,L)
plt.xlabel('time')
plt.ylabel('Angular Momentum')
plt.show()

plt.scatter(t,M)
plt.xlabel('time')
plt.ylabel('Mass')
plt.show()

plt.scatter(t,E)
plt.xlabel('time')
plt.ylabel('KE')
plt.show()

print(M[-1])
print(E[-1])

# for i in range(len(files)):
#     data = athena_read.athdf(files[i])
#     r = data['x1f'][0:64]
#     rho = data['rho'][0, :, :]
#     v_theta = data['vel2'][0, :, :]
#     L.append(rho * r * v_theta)
#     t.append(data['Time'])
#
# L = np.array(L)
# t = np.array(t)
#
# # print(np.shape(L))
# #
# L_tot = []
# i = 0
# while i < len(files):
#     for j in range(64):
#         for k in range(64):
#             sum = 0
#             sum = sum + L[i][j][k]
#     L_tot.append(sum)
#     i = i + 1
#
# L_tot = np.array(L_tot)

# plt.plot(t,L_tot)

#1.3773943




# def f(x,a,b,c):
#     return a*np.sin(b*x+c)
#
# popt,pcov = curve_fit(f,t[1:],dwdt)
#
# a = popt[0]
# b = popt[1]
# c = popt[2]
#
# t_arr = np.arange(t[1],t[-1],0.01)
#
# plt.plot(t_arr,f(t_arr,a,b,c))


# data = athena_read.athdf(files[0])
#
# rows, cols = np.where(data['rho'][0,:,:]==np.max(data['rho'][0,:,:]))
#
# print(data['rho'][0,:,:][rows,cols])
# print(data['x2f'][0:64][rows])


# dmass = density * area element
# dA = r*dr*dtheta
# Angular mom = r*dmass*v_tangential
