import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

import matplotlib.pyplot as plt
#from streamplot import streamplot

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
	return (u[i,j-1] - 2*u[i,j] + u[i,j+1])/(hy**2)

# Fucion lado Derecho
def Pois_RHS(u_mid, v_mid,P, nx, ny, delta, dt, rho):
	C = rho*delta/(2*dt)	# EL delta va por la matriz
	k=0
	RHS = np.zeros(nx*ny)

	for j in range(-1,ny-1):
		for i in range(nx):

			if i==0:
				RHS[k] = ( (-3*u_mid[0,j] + 4*u_mid[1,j] - u_mid[2,j] ) + ( v_mid[i,j+1] - v_mid[i,j-1] ) )*C
			elif i==nx-1:
				RHS[k] = ( ( 3*u_mid[nx-1,j] - 4*u_mid[nx-2,j] + u_mid[nx-3,j] ) + ( v_mid[i,j+1] - v_mid[i,j-1]) )*C
			else:
				RHS[k] = (u_mid[i+1,j] - u_mid[i-1,j] + v_mid[i,j+1] - v_mid[i,j-1])*C
			k +=1
	return RHS

# parte 3 resolucion ecuacion eliptica 
def Pois_Matrix(nx,ny):	#PONER CONDICIONES PERIODICAS !!!!1
	M = np.zeros((nx*ny,nx*ny))
	k = 0
	for j in range(ny):
		for i in range(nx):
			M[k,k] = -4

			# Neumann + Periodica !!! Diferencia centrada para determinar el punto P[-1,j] => P[-1,j] = P[1,j] 
			if i==0 and j==0:
				#M[k,k]    += 4/3

				M[k,k+1]  = 1 + 1
				M[k,k+nx] = 1
				M[k,-nx]   = 1
			elif i==0 and j==ny-1:
				#M[k,k]    += 4/3

				M[k,k+1]  = 1 + 1
				M[k,k-nx] = 1
				M[k,0]    = 1		
			elif i==nx-1 and j==0:
				#M[k,k]      += 4/3

				M[k,k-1]  = 1 + 1
				M[k,k+nx] = 1
				M[k,-1]   = 1
			elif i==nx-1 and j==ny-1:
				#M[k,k]   += 4/3

				M[k,k-1]  = 1 + 1
				M[k,nx-1] = 1
				M[k,k-nx] = 1
			# Solo Neumann
			elif i==0:
				#M[k,k]   += 4/3

				M[k,k+1]  = 1 + 1

				M[k,k+nx] = 1
				M[k,k-nx] = 1
			elif i==nx-1:
				#M[k,k]   += 4/3

				M[k,k-1]  = 1 + 1 
				M[k,k+nx] = 1
				M[k,k-nx] = 1

			# Solo Periodica
			elif j==0:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,k+nx] = 1
				M[k,-nx+i] = 1
			elif j==ny-1:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,i]    = 1
				M[k,k-nx] = 1

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


nu = 0.1
g = 9.81 # m/s**2
rho = 1  #

# Los limites y numero de puntos deben estar de tal forma que dx = dy !!!
nx = 81
ny = (nx-1)*2 + 1
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

# Campo de Velocidades
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))

# Velocidades Intermedias
u_mid = np.zeros((nx,ny))
v_mid = np.zeros((nx,ny))

# Velocidades Iteracion Anterior
u_n = np.zeros((nx,ny))
v_n = np.zeros((nx,ny))

# Presion como matriz
P = np.zeros((nx,ny))

# Presion como vector
P_vec = np.zeros(nx*ny)

# Crear Matriz del Sistema
A = Pois_Matrix(nx,ny)
A_s = sparse.csr_matrix(A)


error = 1
iteracion = 0
for paso in range(10):
	# Euler Explplicito con diferencia atrasada
	for j in range(-1,ny-1):			# Se parte desde el ultimo punto para facilitar la periodicidad
		for i in range(1,nx-1):			# Puntos 0 y nx-1 ya fueron calculados (u,v=0,0)
			u_mid[i,j] =  u[i,j] - dt*( u[i,j]*Dx_a(u,i,j,dx) + v[i,j]*Dy_a(u,i,j,dy) )
			u_mid[i,j] += dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) )

			v_mid[i,j] =  v[i,j] - dt*( u[i,j]*Dx_a(u,i,j,dx) + v[i,j]*Dy_a(u,i,j,dy) )
			v_mid[i,j] += dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) ) - dt*g

	# Lado derecho del sistema
	B = Pois_RHS(u_mid, v_mid, P, nx, ny, dx, dt, rho)

	# Resolver Sistema
	P_vec = spsolve(A_s,B)

	# Reescribir la presion como matriz
	k = 0
	for j in range(ny):
		for i in range(nx):
			P[i,j] = P_vec[k]
			k += 1

	# Actualizar Velocidades
	u_n = u.copy()
	v_n = v.copy()
	for j in range(-1,ny-1):
		for i in range(1,nx-1):
			u[i,j] = u_mid[i,j] + dt/(2*dx*rho)*(P[i+1,j]-P[i-1,j])
			v[i,j] = v_mid[i,j] + dt/(2*dy*rho)*(P[i,j+1]-P[i,j-1])

	error = np.max(abs(v-v_n))/np.max(abs(v))
	iteracion += 1
	print("It: ", iteracion, "Error: ", error)
	#if iteracion > 1000:
	#	break


print( "Numero de Iteraciones: ", iteracion)

print(v[:,(ny-1)/2])

u = np.transpose(u)
v = np.transpose(v)

X, Y = np.meshgrid(X,Y)
fig = plt.figure(figsize = (4,8), dpi=100 )
plt.quiver(X[::10, ::2], Y[::10, ::2], u[::10, ::2], v[::10, ::2], headaxislength=4);
#plt.quiver(X, Y, u, v);
#plt.streamplot(X, Y, u, v, color=u, linewidth=5*P/P.max())
plt.show();
