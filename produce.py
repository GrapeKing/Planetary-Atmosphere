from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter

r = np.log(np.loadtxt("coordinate.txt"))

rho = np.log(np.loadtxt("density.txt"))

rho_ana = np.log(0.0005)+4/np.exp(r)-1/4;

fig, ax = plt.subplots()

line0, = ax.plot(r, rho_ana,label='Analytic')
line1, = ax.plot(r, rho[-1],label='Numerical')

ax.legend()

plt.xlabel('ln[r]')
plt.ylabel('ln[Rho]')

def animate(i):
	line1.set_ydata(rho[i])
	return line0, line1,

ani = anim.FuncAnimation(
	fig, animate, interval = 1, save_count = 20)

plt.show()