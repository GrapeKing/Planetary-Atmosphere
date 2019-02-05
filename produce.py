from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter
from scipy import optimize as opt
from scipy import integrate as intg

save = True

rb = 4
ms = rb
# B.g. density
roubg = 0.00057
rou0 = roubg*np.exp(-1/4)

# No. of steps
div = 128

# differentiation approx order
o = 3

# Core radius = 1
furthest_point = 16
upper_bound = np.log(furthest_point)

#Initializing Grid
ln_r = np.linspace(0,upper_bound,div)

#Analytic solution
ana = rb*np.exp(-ln_r)+np.log(rou0)
mass = np.ones(div)
for i in range(div):
	mass[i] = intg.quad(lambda r: rou0*4*np.pi*(r**2)*np.exp(rb/r), 1, np.exp(ln_r[i]))[0]

guess_rou = ana
guess_m = mass + ms
guess = np.append(guess_rou,guess_m)

def dif_matrix(x,j):
	#First order approximation
	if j == 1:
		l=len(x)
		dx = x[1]-x[0]
	
		#unit matrix with size l
		dif_n = np.diag(np.ones(l))
	
		#list of (-1)
		dif_pre_ones = np.ones(l-1)*(-1)
		
		#diagonal matrix of -1 shifted to the left by 1 element
		dif_pre = np.diag(dif_pre_ones,k=-1)
	
		dif = dif_n+dif_pre
		'''its gonna be sth look like this
		[[1,0,0,0,0]
		[-1,1,0,0,0]
		[0,-1,1,0,0]
		[0,0,-1,1,0]
		[0,0,0,-1,1]]'''
		dif[0,0]=-1
		dif[0,1]=1
	
		dif = dif/dx

	#Second order approximation
	elif j == 2:

		l=len(x)
		dx = x[1]-x[0]
	
		dif_n_ones = np.ones(l-1)
		dif_n = np.diag(dif_n_ones,k=1)

		dif_pre_ones = np.ones(l-1)*(-1)

		dif_pre = np.diag(dif_pre_ones,k=-1)
	
		dif = (dif_n+dif_pre)
		'''its gonna be sth look like this
		[[0,1,0,0,0]
		[-1,0,1,0,0]
		[0,-1,0,1,0]
		[0,0,-1,0,1]
		[0,0,0,-1,0]]'''
	
		dif[0,0]=-1
	
		dif[l-1,l-1]=1
	
		dif = dif/(2*dx)

		dif[0]=dif[0]*2

		dif[l-1]=dif[l-1]*2

	#Four point approximation
	elif j == 3:

		l=len(x)
		dx = x[1]-x[0]

		dif_8 = np.diag(np.ones(l-1)*8,k=1)

		dif_n8 = np.diag(np.ones(l-1)*-8,k=-1)

		dif_n1 = np.diag(np.ones(l-2)*-1,k=2)

		dif_1 = np.diag(np.ones(l-2),k=-2)

		dif = dif_8 + dif_n8 + dif_1 + dif_n1

		dif = dif /(12*dx)

		dif[0,0]=-1
		dif[0,1]=1
		dif[0,2]=0
		dif[0]=dif[0]/dx

		dif[1,0]=-1
		dif[1,1]=0
		dif[1,2]=1
		dif[1,3]=0
		dif[1]=dif[1]/(2*dx)

		dif[l-1,l-1]=1
		dif[l-1,l-2]=-1
		dif[l-1,l-3]=0
		dif[l-1]=dif[l-1]/dx

		dif[l-2,l-1]=1
		dif[l-2,l-2]=0
		dif[l-2,l-3]=-1
		dif[l-2,l-4]=0
		dif[l-2]=dif[l-2]/(2*dx)

	#Three point approximation
	elif j == 4:

		l=len(x)
		dx = x[1]-x[0]

		dif_1 = np.diag(np.ones(l-2),k=2)

		dif_3 = np.diag(np.ones(l)*3)

		dif_n4 = np.diag(np.ones(l-1)*-4,k=-1)

		dif = dif_1 + dif_3 + dif_n4

		dif = dif /(6*dx)

		dif[0,0]=-1
		dif[0,1]=1
		dif[0,2]=0
		dif[0]=dif[0]/dx

		dif[-1,-1]=1
		dif[-1,-2]=-1
		dif[-1,-3]=0
		dif[-1]=dif[-1]/dx

		dif[-2,-1]=1
		dif[-2,-2]=0
		dif[-2,-3]=-1
		dif[-2,-4]=0
		dif[-2]=dif[-2]/(2*dx)

	return dif

mat = dif_matrix(ln_r,o)

#Function to be reduced to zero
def zero(y):

	l = len(y)/2

	ln_rou = y[:l]
	m = y[l:]

	fun1 = mat.dot(ln_rou)+m*np.exp(-ln_r)
	fun2 = mat.dot(m)-4*np.pi*np.exp(3*ln_r)*np.exp(ln_rou)

	fun = np.append(fun1,fun2)
	fun = np.append(fun,ln_rou[-1]-np.log(roubg))
	fun = np.append(fun,m[0]-ms)
	return fun

#Solve for ODE
sol = opt.root(zero,guess,method='lm')

res = sol.x


res_rou = res[:div]
res_m = res[div:]

fig, ax = plt.subplots()

line0, = ax.plot(ln_r,res_rou,label="massive profile")

r = np.log(np.loadtxt("coordinate.txt"))

rho = np.log(np.loadtxt("density.txt"))

line1, = ax.plot(r, rho[-1],label='Numerical')
line2, = ax.plot(r, 4/np.exp(r)+np.log(rou0),label='massless')

ax.legend()

plt.xlabel('ln[r]')
plt.ylabel('ln[Rho]')

def animate(i):
	line1.set_ydata(rho[i])
	return line0, line1,

ani = anim.FuncAnimation(
	fig, animate, interval = 1, save_count = 1000)

if save == True:
	Writer = anim.writers['ffmpeg']
	writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=1800)
	ani.save("0_density.mp4", writer=writer)
else:
	plt.show()
