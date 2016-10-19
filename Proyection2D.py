import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# Derivadas como funciones:
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

# Solucion analitica del estado estacionario
def Vel_analitica(t,g,L,nu,nx):
	v = np.zeros(nx)
	x = np.linspace(0,L,nx)
	n_max = 100
	for n in range(1,n_max):
		C = 2*g*L/(nu*(n*np.pi)**2) * (1-np.exp(-nu*(n*np.pi/L)**2*t))
		print 'const: ', np.shape(C)
		v +=  C* np.sin(n*np.pi/L*x)
		print 'vel: ', np.shape(v) 
	return v

# Fucion lado Derecho
def Pois_RHS(u_mid, v_mid,P, nx, ny, delta, dt, rho):
	C = rho*delta/(2*dt)
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

# Definir Matriz para las Presiones
def Pois_Matrix(nx,ny):
	M = np.zeros((nx*ny,nx*ny))
	k = 0
	for j in range(ny):
		for i in range(nx):
			M[k,k] = -4

		# Neumann + Periodica 
		 # Diferencia centrada en P[0,j] => P[-1,j] = P[1,j] 
			if i==0 and j==0:

				M[k,k+1]  = 1 + 1
				M[k,k+nx] = 1
				M[k,-nx]  = 1
			elif i==0 and j==ny-1:

				M[k,k+1]  = 1 + 1
				M[k,k-nx] = 1
				M[k,0]    = 1
			elif i==nx-1 and j==0:

				M[k,k-1]  = 1 + 1
				M[k,k+nx] = 1
				M[k,-1]   = 1
			elif i==nx-1 and j==ny-1:

				M[k,k-1]  = 1 + 1
				M[k,nx-1] = 1
				M[k,k-nx] = 1

		# Solo Neumann
			elif i==0:

				M[k,k+1]  = 1 + 1

				M[k,k+nx] = 1
				M[k,k-nx] = 1
			elif i==nx-1:

				M[k,k-1]  = 1 + 1 
				M[k,k+nx] = 1
				M[k,k-nx] = 1

		# Solo Periodica
			elif j==0:
				M[k,k+1]   = 1
				M[k,k-1]   = 1
				M[k,k+nx]  = 1
				M[k,-nx+i] = 1
			elif j==ny-1:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,i]    = 1
				M[k,k-nx] = 1

		# Resto de la malla
			else:
				M[k,k+1]  = 1
				M[k,k-1]  = 1
				M[k,k+nx] = 1
				M[k,k-nx] = 1

		# Condicion para fijar Presion
			if j==(ny-1)/2 and i==(nx-1)/2:
				M[k,k] += -1
			k +=1

	M_s = sparse.csr_matrix(M)
	return M_s


nu = 0.1 # m^2/s
g = 9.81 # m/s**2
rho = 1  # kg/m^3

# Los limites y numero de puntos deben estar de tal forma que dx = dy !!!
nx = 41
ny = (nx-1)*2 + 1
xlim = [0,2]
ylim = [0,4]

X = np.linspace(xlim[0],xlim[1], nx)
Y = np.linspace(ylim[0],ylim[1], ny)
dx = X[1]-X[0]
dy = Y[1]-Y[0]

dt = 0.01

if nu*dt/dx**2>0.5:
	print 'Condicion de difusion inestable: ', nu*dt/dx**2

# Campo de Velocidades
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))

# Velocidades Intermedias
u_mid = np.zeros((nx,ny))
v_mid = np.zeros((nx,ny))

# Velocidades Iteracion Anterior
u_n = np.zeros((nx,ny))
v_n = np.zeros((nx,ny))

# Presion como matriz y como vector
P = np.zeros((nx,ny))
P_vec = np.zeros(nx*ny)


# Crear Matriz del Sistema
A_s = Pois_Matrix(nx,ny)

# Arreglos
cycle = np.arange(-1,ny-1)
inter = np.arange(1,nx-1)

# Solucion Analitica
V_a = Vel_est(g,xlim[1],nu,nx)

error = 10
iteracion = 0
while error>1:

	if np.max(v)*dt/dy>1:
		print 'Condicion de conveccion inestable: ', np.max(v)*dt/dy

	# Euler Explplicito con diferencia atrasada
	for j in range(-1,ny-1):			# Se parte desde el ultimo punto para facilitar la periodicidad
		for i in range(1,nx-1):			# Puntos 0 y nx-1 no se calculan (u,v=0,0)
			u_mid[i,j] =  u[i,j] - dt*( u[i,j]*Dx_a(u,i,j,dx) + v[i,j]*Dy_a(u,i,j,dy) )
			u_mid[i,j] += dt*nu*( DDx(u,i,j,dx) + DDy(u,i,j,dy) )

			v_mid[i,j] =  v[i,j] - dt*( u[i,j]*Dx_a(v,i,j,dx) + v[i,j]*Dy_a(v,i,j,dy) )
			v_mid[i,j] += dt*nu*( DDx(v,i,j,dx) + DDy(v,i,j,dy) ) - dt*g

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

	u[1:-1,:] = u_mid[1:-1,:] + dt/(2*dx*rho)*( P[2:,:] - P[:-2,:] )
	v[:,cycle] = v_mid[:,cycle] + dt/(2*dy*rho)*( P[:,cycle+1]  - P[:,cycle-1] )

	# Variacion con respecto a la iteracion anterior y error con respecto a la solucion analitica
	variacion = sum(np.abs(v_n[1:-1,(ny-1)/2]-v[1:-1,(ny-1)/2])/np.abs(v_n[1:-1,(ny-1)/2]))
	error = sum(np.abs(V_a[1:-1] - v[1:-1:,(ny-1)/2])/np.abs(V_a[1:-1]))
	print "Iteracion: ", iteracion, "Variacion: ", variacion, "Error: ", error
	iteracion += 1


fig1 = plt.figure(figsize = (4,4), dpi=100 )
plt.plot(X,V_a, c='k',ls='-');

u = np.transpose(u)
v = np.transpose(v)
P = np.transpose(P)

X, Y = np.meshgrid(X,Y)
fig2 = plt.figure(figsize = (4,8), dpi=100 )
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 8
plt.contourf(X,Y,P);
plt.colorbar();
plt.quiver(X[::20,::2], Y[::20,::2], u[::20,::2], v[::20,::2], headaxislength=5);
plt.xlabel('X')
plt.ylabel('Y');
plt.show();
