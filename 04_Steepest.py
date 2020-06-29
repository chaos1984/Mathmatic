import numpy as np
import sympy as sp
from sympy.abc import x,y
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def residual_initial(A,b,x0):
	return b-A.dot(x0)

def residual(A,alpha,r0):
    return r0-alpha*A.dot(r0)

def stepsize(A,r0):
    return r0.dot(r0.T)/r0.dot(A).dot(r0.T)

if __name__ == "__main__":
	A = np.array([[3.,2],[2,6]])
	b = np.array([2.,-8])
	x0 = [-2,-2]
	r = residual_initial(A,b,x0)
	maxerr =1E9
	i = 1
	x_loc = []
	while maxerr > 0.01:
		alpha = stepsize(A,r)
		err = alpha*r
		x0 = x0 + err
		maxerr =max(err)
		print ("Run %d\nErr %10.5f" %(i,max(err)))
		print ('X:',x0)
		r = r-alpha*A.dot(r.T)
		i += 1
		x_loc.append(x0.tolist())
	print (x0)

