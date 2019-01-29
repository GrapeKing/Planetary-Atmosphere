import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter

r = np.log(np.loadtxt("coordinate.txt"))

rho = np.log(np.loadtxt("density.txt"))

fig, ax = plt.subplots()

line, = ax.plot(r, np.linspace(5,0,len(r)),label='Numerical')
line2, = ax.plot(r,4/np.exp(r),label='Analytic')
ax.legend()

plt.xlabel('ln[r]')
plt.ylabel('ln[Rho]')

def animate(i):
	line.set_ydata(rho[i])
	return line, line2,

ani = anim.FuncAnimation(
	fig, animate, interval = 1, save_count = 2000)

Writer = anim.writers['ffmpeg']
writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=1800)
ani.save("0_density.mp4", writer=writer)
