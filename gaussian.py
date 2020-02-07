import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import matplotlib.animation as animation
import cmath

pi = np.pi
h_x = 0.01
h_t = 0.5			

x0 = 0.0
xe = 300.0
L = abs(xe-x0)
N = int(abs(xe-x0)/h_x)

sigma = L/20.0

V0 = -0.50
lamda = 2.0*h_x**2/h_t

spread = 1.5
total_time = lamda*(sigma**2)*((spread**2 - 1.0)**0.5)/(4.0*h_x**2)

k = (N+1)*lamda/(8*h_x*total_time)#n*pi/L

xv = 130.0

x = np.arange(x0,xe,h_x)

#	Intial psi
psi = np.exp(1j*k*x)*np.exp(-(x-L/4)**2/(2.0*sigma**2))
psi[0] = complex(0.0,0.0)
psi[N-1] = complex(0.0,0.0)

#	Potential
V = []
for i in range(N):
	#if 	(x0 + i*h_x) >= xv:
	#	V.append(V0)
	if (x0 + i*h_x) >= xv and (x0 + i*h_x) <= xv + 50:
		V.append(V0)#((i*h_x-xv)/100)**0.5)
	else:
		V.append(0.0)
	#V.append(((150.0 - i*h_x)/37.5)**2)

# 	Omega
def om(psi):
	global V
	oma = [-psi[1]]
	for i in range(1,N):
		if i != N-1:
			oma1 = psi[i+1]
		else:
			oma1 = 0.0
		oma.append(-oma1 + (2.0 + h_x**2 * V[i] + 1j*lamda)*psi[i] - psi[i-1])
	z=np.arange(N)
	return oma
omega = om(psi)

def v1v2(psi,omega):
	global V
	v1 = [complex(0.0,0.0)]
	v2 = [complex(0.0,0.0)]
	v1.append(2.0 + h_x**2 * V[1] - 1j*lamda)
	v2.append(omega[1])
	v1[0] = 0.0
	v2[0] = v1[0]*(v2[1] - omega[1])
	for i in range(2,N):
		v1.append(2.0 + h_x**2 * V[i] - 1j*lamda - 1/v1[i-1])
		v2.append(omega[i] + v2[i-1]/v1[i-1])
	return v1,v2
v1,v2 = v1v2(psi,omega)

fig= figure()
fig.suptitle("Gaussian Wave Packet")

ax01 = subplot2grid((2,1), (0, 0))
ax03 = subplot2grid((2,1), (1, 0))

#ax01.set_title('Real $\psi$ vs x')
#ax02.set_title('Imaginary $\psi$ vs x')
#ax03.set_title('Probability vs x')

ax01.set_ylabel("Real part of $\psi$")
ax03.set_xlabel("x")
ax03.set_ylabel("Proability")

wave1, pot1 = ax01.plot(x,psi.real,x,V)
wave1.set_label('$\psi$')
pot1.set_label('V')
ax01.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
ax01.set_ylim([-2,2])

wave2, pot2 = ax03.plot(x,abs(psi)**2,'xkcd:black',x,V)
wave2.set_label('Prob.')
pot2.set_label('V')
ax03.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
ax03.set_ylim([V0 - 0.5,3.8])
#wave, = ax[0].plot(x,v1.real)

def animate(t):
	if t == 0:
		return wave1,wave2
	global psi,v1,v2,omega,x

	psi[0] = complex(0.0)
	psi[N-1] = complex(0.0)
	psi[N-2] = -v2[N-2]/v1[N-2]
	for i in range(N-3,0,-1):
		psi[i] = (psi[i+1] - v2[i])/v1[i]
	
	#print t

	omega = om(psi)

	v1,v2 = v1v2(psi,omega)
	
	wave1.set_data(x,psi.real)
	wave2.set_data(x,abs(psi)**2)

	return wave1, wave2

ani = animation.FuncAnimation(fig,animate,None, interval = 1)

plt.show()

