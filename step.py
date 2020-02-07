import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

E = 2.0
Vc = 1.5		
V_x = 0
m = 1.0
h_bar = 1.0
x_le = -10.0
x_re = 10.0

k1 = (2.0*m*E)**0.5/h_bar
if Vc > E:
	k2 = (2.0*m*(Vc-E))**0.5/h_bar
	w_exp = True
else:
	k2 = (2.0*m*(E-Vc))**0.5/h_bar
	w_exp = False

w1 = E
w2 = E

A = 1.0
if w_exp == False:
	B = A*(k1-k2)*np.exp(2.0*1j*k1*V_x)/(k1+k2)
	C = A*2.0*k1*np.exp(1j*(k1-k2)*V_x)/(k1 + k2)
else:
	B = A*(k1 - 1j * k2)*np.exp(2.0*1j*k1*V_x)/(k1 + 1j * k2)
	C = A*2.0*k1*np.exp((1j*k1+k2)*V_x)/(k1 + 1j * k2)

h = 0.001

x_l = np.arange(x_le,V_x,h)
x_r = np.arange(V_x,x_re,h)
x = np.arange(x_le,x_re,h)

n = len(x)


psi_i = A*np.exp(1j*k1*x_l)

psi_r = B*np.exp(-1j*k1*x_l)

if w_exp == False:
	psi_t = C*np.exp(1j*k2*x_r)
else:
	psi_t = C*np.exp(-k2*x_r)

psi_c = psi_i + psi_r

def potential(x):
	if x  >= V_x:
		return Vc
	else:
		return 0.0
	
V = []
for i in range(n):
	V.append(potential(x_le + i*h))

fig, ax = plt.subplots(2,1)
wave_c, wave_t, pot = ax[0].plot(x_l,psi_c.real, x_r,psi_t.real, x,V)

def animate(t):
	global psi_c
	global psi_t

	psi_c_t = psi_c*np.exp(-1j*w1*t/100)
	wave_c.set_ydata(psi_c_t.real)
	
	psi_t_t = psi_t*np.exp(-1j*w1*t/100)
	wave_t.set_ydata(psi_t_t.real)

	#print psi_c_t[n1-1]
	return wave_c,wave_t,

#def init():
#	psi_i.set_data([],[]) 
#	psi_r = np.linspace(0.0,0.0,int(n1)) + 1j * np.linspace(0.0,0.0,int(n1))
#	psi_t = np.linspace(0.0,0.0,int(n1)) + 1j * np.linspace(0.0,0.0,int(n1))
#	ax.vlines(10,0,5)
#	pot = ax.plot(x,V)
#	return wave_r,

ani = animation.FuncAnimation(fig, animate,None,interval=50)

plt.show()