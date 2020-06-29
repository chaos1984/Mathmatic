import numpy as np
import sympy as sp
# from sympy.abc import x,y
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def residual_initial(A,b,x0):
	return b-A.dot(x0)

def residual(A,alpha,r0):
    return r0-alpha*A.dot(r0)

def stepsize(A,r0):
    return r0.dot(r0.T)/r0.dot(A).dot(r0.T)

def threeDplot(X,Y,Z,data):
	fig = plt.figure()
	ax = Axes3D(fig) 
	ax.plot3D(data[:,0],data[:,1],data[:,2],'r')
	ax.scatter(data[:,0],data[:,1],data[:,2],marker='*',color='r')
	ax.plot_surface(X, Y, Z, rstride = 10, cstride = 10, cmap = plt.get_cmap('ocean'),alpha=0.7)
#	ax.plot_wireframe(X, Y, Z, rstride = 10, cstride = 10, cmap = plt.get_cmap('rainbow'))
	ax.contour(X, Y, Z, zdir = 'z', offset = -2, cmap = plt.get_cmap('rainbow'))
	ax.set_alpha(0.9)
	plt.show()


def Quadratic(A,b,X,Y):
	return 0.5*(A[0,0]*X*X+(A[0,1]+A[1,0])*X*Y+A[1,1]*Y*Y)-b[0]*X-b[1]*Y


if __name__ == '__main__':
	
	A = np.array([[3,2],[2,6]])
	b = np.array([2.,-8])
	data = np.linspace(-4,4,100)
	X, Y = np.meshgrid(data,data)
	Z = Quadratic(A,b,X,Y)


	x0 = [-4,-2]
	r = residual_initial(A,b,x0)
	maxerr =1E9
	i = 1
	x_loc = [[x0[0],x0[1],Quadratic(A,b,x0[0],x0[1])]]
	while maxerr > 0.01:
		alpha = stepsize(A,r)
		print ('alpha',alpha)
		err = alpha*r   #err == stepsize * direction
		x0 = x0 + err
		z = Quadratic(A,b,x0[0],x0[1])
		maxerr =max(err)
		print ("Run %d\nErr %10.5f" %(i,max(err)))
		print ('X:',x0,'Z',z)
		r = r-alpha*A.dot(r.T) 
		i += 1
		temp = x0.tolist()
		temp.append(z)
		x_loc.append(temp)
threeDplot(X,Y,Z,np.array(x_loc))
