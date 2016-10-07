import numpy as np
from matplotlib import pyplot

# Operadores
# Primera Derivada Atrasada
def Dx_a(u,i,j,hx):
	return (u[i,j] - u[i-1,j])/(hx)
def Dy_a(u,i,j,hy):
	return (u[i,j] - u[i,j-1])/(hy)
# Primera Derivada Centrada
def Dx_c(u,i,j,hx):
	return (u[i+1,j] - u[i-1,j])/(2*hx)
def Dy_c(u,i,j,hy):
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
g = 9.81 # m/s**2
rho = 1  #

# Los limites y numero de puntos deben estar de tal forma que dx = dy !!!
nx = 41
ny = 81
xlim = [0,2]
ylim = [0,4]
tlim = [0,5]

X = np.linspace(xlim[0],xlim[1], nx)
Y = np.linspace(ylim[0],ylim[1], ny)

dx = X[1]-X[0]
dy = Y[1]-Y[0]
dt = 0.1			#dx**2/(2*nu)
Nt = int(tlim[1] / dt) + 1 
t = np.linspace(tlim[0],tlim[1], Nt)

# Obtenibles explicitamente
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))
u_mid = np.zeros((nx,ny))
v_mid = np.zeros((nx,ny))

p = np.zeros(nx*ny)
p_tem = np.zeros(nx*ny)

# PREDICCION u*, v*
# Condiciones de contorno no deslizamiento (x=0,L):
#u_mid[ 0  ,:] = 0
#u_mid[nx-1,:] = 0
#v_mid[ 0  ,:] = 0
#v_mid[nx-1,:] = 0

# Euler Explplicito con diferencia centrada
	# CONDICION PERIODICA EN (Y=0,L)
for i in range(1,nx-1):			# Puntos 0 y nx-1 ya fueron calculados
	for j in range(-1,ny-1):		# Se parte desde el ultimo punto para lograr periodicidad
		u_mid[i,j] =  u[i,j] - dt*( u[i,j]*Dx_a(u,i,j,dx) + v[i,j]*Dy_a(u,i,j,dy) )
		u_mid[i,j] += dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) )

		v_mid[i,j] =  v[i,j] - dt*( u[i,j]*Dx_a(u,i,j,dx) + v[i,j]*Dy_a(u,i,j,dy) )
		v_mid[i,j] += dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) ) - dt*g

def Pois_RHS(u_mid, v_mid,P, nx, ny, delta, dt, rho):
	C = -delta*rho/(2*dt)
	k=0
	for i in range(1,nx-1):
		for j in range(1,ny-1):
			# poner casos
			RHS[k] = (u_mid[i+1,j] - u_mid[i-1,j] + v_mid[i,j+1] - v_mid[u,j-1])*C
	return RHS

def Pois_Matrix(nx,ny):	#PONER CONDICIONES PERIODICAS !!!!1
	M = np.zeros((nx*ny,nx*ny))
	k = 0
	for j in range(ny):
		for i in range(nx):
			M[k,k] = -4

			# Neumann + Periodica !!! Stencil orden 2 para bordes con Neumann
			if i==0 and j==0:
				M[k,k]          += 4/3

				M[k,k+1]        = 1 - 1/3
				M[k,k+nx]       = 1
				M[k,nx*(ny-1)]  = 1 		##### Revisar
			elif i==0 and j==ny-1:
				M[k,k]    += 4/3

				M[k,k+1]  = 1 - 1/3
				M[k,k-nx] = 1
				M[k,0]    = 1		
			elif i==nx-1 and j==0:
				M[k,k]      += 4/3

				M[k,k-1]     = 1 - 1/3
				M[k,k+nx]    = 1
				M[k,nx*ny-1] = 1
			elif i==nx-1 and j==ny-1:
				M[k,k]   += 4/3

				M[k,k-1]  = 1 - 1/3
				M[k,k-nx] = 1
				M[k,nx-1] = 1
			# Solo Neumann
			elif i==0:
				M[k,k]   += 4/3

				M[k,k+1]  = 1 - 1/3##
				M[k,k+nx] = 1
				M[k,k-nx] = 1
			elif i==nx-1:
				M[k,k]   += 4/3

				M[k,k-1]  = 1 - 1/3##
				M[k,k+nx] = 1
				M[k,k-nx] = 1

			# Solo Periodica
			elif j==0:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,k+nx] = 1
				M[k,(ny-1)*nx+i] = 1
			elif j==ny-1:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,k-nx] = 1
				M[k,i]    = 1
			else:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,k+nx] = 1
				M[k,k-nx] = 1
			# Condicion para fijar Presion
			if j==(ny-1)/2 and i==(nx-1)/2:
				M[k,k] += -1
			k +=1
	return M

A = Pois_Matrix(nx,ny)
B = Pois_RHS(u_mid, v_mid, p, nx, ny, dx, dt, rho)
print np.shape(A)

# Ecuacion de Poisson con Map: 2D -> 1D
#k = 0
#for i in range(1,N-1):
#	for j in range(1,N-1):
		

