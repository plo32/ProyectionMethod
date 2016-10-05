import numpy as np
from matplotlib import pyplot

# Operadores
# Primera Derivada Centrada
def Dx(u,i,j,hx):
	return (u[i+1,j] - u[i-1,j])/(2*hx)
def Dy(u,i,j,hy):
	return (u[i,j+1] - u[i,j-1])/(2*hy)
# Segunda Derivada Centrada
def DDx(u,i,j,hx):
	return (u[i-1,j] - 2*u[i,j] + u[i+1,j])/(hx**2)
def DDy(u,i,j,hy):
	return (u[i,j-1] - 2*u[i,j] + u[i,j-1])/(hy**2)
# Condicion Neumann Orden 2 (Cuidado con Frwd {+dx} o Bkwd {-dx})
def Nmn2(u_1,u_2,du,dx):
	return 4/3*u_1 - 1/3*u_2 + 2/3*du*dx

nu = 1.785*10**(-6)

N  = 20
xlim = [0,1]
ylim = [0,1]
tlim = [0,5]

X = np.linspace(xlim[0],xlim[1], N)
Y = np.linspace(ylim[0],ylim[1], N)

dx = X[1]-X[0]
dy = Y[1]-Y[0]
dt = 0.1			#dx**2/(2*nu)
Nt = int(tlim[1] * dt) + 1 
t = np.linspace(tlim[0],tlim[1], Nt)

# Obtenibles explicitamente
u = np.zeros((N,N))
v = np.zeros((N,N))
u_mid = np.zeros((N,N))
v_mid = np.zeros((N,N))

p = np.zeros(N**2)
p_tem = np.zeros(N**2)

# PREDICCION u*, v*
# Condiciones de contorno no deslizamiento (x=0,L):
u_mid[ 0 ,:] = 0
u_mid[N-1,:] = 0
v_mid[ 0 ,:] = 0
v_mid[N-1,:] = 0

# Euler Explplicito con diferencia centrada
	# CONDICION PERIODICA EN (Y=0,L)
for i in range(1,N-1):			# Puntos 0 y N-1 ya fueron calculados
	for j in range(-1,N-1):		# Se parte desde el ultimo punto para lograr periodicidad
		u_mid[i,j] =  u[i,j] - dt*( u[i,j]*Dx(u,i,j,dx) + v[i,j]*Dy(u,i,j,dy) )
		u_mid[i,j] -= dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) )

		v_mid[i,j] =  v[i,j] - dt*( u[i,j]*Dx(u,i,j,dx) + v[i,j]*Dy(u,i,j,dy) )
		v_mid[i,j] -= dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) )

def Pois_RHS(u_mid, v_mid,P, n, delta, dt, rho):
	C = delta*rho/(2*dt)
	k=0
	for i in range(1,n-1):
		for j in range(1,n-1):
			RHS[k] = (u_mid[i+1,j] - u_mid[u-1,j] + v_mid[i,j+1] - v_mid[u,j-1])*C
	return RHS

def Pois_Matrix(N):	#PONER CONDICIONES PERIODICAS !!!!1
	n = N-2
	M = np.zeros((n**2,n**2))
	k = 0
	for j in range(n):
		for i in range(n):
			M[k,k] = -4

			# Neumann + Periodica !!! Stencil orden 2 para bordes con Neumann
			if i==0 and j==0:
				M[k,k] += 4/3
				M[k,k+1] = 1 - 1/3
				M[k,k+n] = 1
				M[k,-1-n]  = 1 		##### Falta poner la condicion periodica
			elif i==0 and j==n-1:
				M[k,k] += 4/3
				M[k,k+1] = 1 - 1/3
				M[k,k-n] = 1
				M[k,-n+j]  = 1		##### Falta poner la condicion periodica
			elif i==n-1 and j==0:
				M[k,k] += 4/3
				M[k,k-1] = 1 - 1/3
				M[k,k+n] = 1
			elif i==n-1 and j==n-1:
				M[k,k] += 4/3
				M[k,k-1] = 1 - 1/3
				M[k,k-n] = 1

			# Solo Neumann
			elif i==0:
				M[k,k] += 4/3
				M[k,k+1] = 1
				M[k,k+n] = 1
				M[k,k-n] = 1
			elif i==n-1:
				M[k,k] += 4/3
				M[k,k-1] = 1
				M[k,k+n] = 1
				M[k,k-n] = 1

			# Solo Periodica
			elif j==0:
				M[k,k+1] = 1
				M[k,k-1] = 1
				M[k,k+n] = 1
			elif j==n-1:
				M[k,k+1] = 1
				M[k,k-1] = 1
				M[k,k-n] = 1
			else:
				M[k,k+1] = 1
				M[k,k-1] = 1
				M[k,k+n] = 1
				M[k,k-n] = 1
			k +=1
	return M

A = Pois_Matrix(5)
print A
# Ecuacion de Poisson con Map: 2D -> 1D
#k = 0
#for i in range(1,N-1):
#	for j in range(1,N-1):
		
