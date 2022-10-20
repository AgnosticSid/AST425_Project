# Creating animation
import sys
sys.path.insert(0, '/home/ubuntu/athena/vis/python')
import athena_read
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from glob import glob
import natsort

files = glob("/home/ubuntu/athena/bin/disk.out1.*.athdf")
files = natsort.natsorted(files)

data = athena_read.athdf(files[0])
x = data['x1v']
y = data['x2v']
z = np.transpose(data['rho'][0])
mycmap = plt.get_cmap('viridis') #plt.get_cmap('gist_heat')
cs = plt.contourf(y, x, z, cmap=mycmap)

z = np.transpose(data['rho'][0])
cs1 = plt.contourf(y, x, z, cmap=mycmap)
z = np.transpose(data['vel1'][0])
cs2 = plt.contourf(y, x, z, cmap=mycmap)
z = np.transpose(data['vel2'][0])
cs3 = plt.contourf(y, x, z, cmap=mycmap)
z = np.transpose(data['rho'][0] * data['x1v'] * data['vel2'][0])
cs4 = plt.contourf(y, x, z, cmap=mycmap)

fig, ax = plt.subplots(2,2,subplot_kw=dict(projection='polar'))
line, = ax[0,0].plot([],[],lw=2)

# ax[0,0].grid(False)
# ax[0,1].grid(False)
# ax[1,0].grid(False)
# ax[1,1].grid(False)

fig.colorbar(cs1, ax=ax[0, 0], pad = 0.2)
fig.colorbar(cs2, ax=ax[0, 1], pad = 0.2)
fig.colorbar(cs3, ax=ax[1, 0], pad = 0.2)
fig.colorbar(cs4, ax=ax[1, 1], pad = 0.2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

def animate(i):
    ax[0,0].cla()
    ax[0,1].cla()
    ax[1,0].cla()
    ax[1,1].cla()

    data = athena_read.athdf(files[i])
    r = data['x1f'][0:64]
    theta = data['x2f'][0:64]
    rho = data['rho'][0,:,:]
    v_r = data['vel1'][0,:,:]
    v_theta = data['vel2'][0, :, :]
    mycmap = plt.get_cmap('viridis')

    ax[0,0].contourf(theta, r, rho.T, cmap=mycmap)
    ax[0,1].contourf(theta, r, v_r.T, cmap=mycmap)
    ax[1,0].contourf(theta, r, v_theta.T, cmap=mycmap)
    ax[1,1].contourf(theta, r, (rho * r * v_theta).T, cmap=mycmap)
    # ax[0, 0].grid(False)
    # ax[0, 1].grid(False)
    # ax[1, 0].grid(False)
    # ax[1, 1].grid(False)
    time = data['Time']
    ax[0,0].set_title('Density')
    ax[0,1].set_title('v_r')
    ax[1,0].set_title('v_theta')
    ax[1,1].set_title('Angular Momentum')
    fig.suptitle('Time = '+str(time))
    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.8,wspace=0.6,hspace=0.6)
    print(i)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(files), interval=100, blit=True)
anim.save('Animation_velcyl.mp4', fps=1)

# phi,v_phi = np.loadtxt('v_theta.txt',float,delimiter=',',unpack=True)
# plt.scatter(phi,v_phi)
# plt.show()