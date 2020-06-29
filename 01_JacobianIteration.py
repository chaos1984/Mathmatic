# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:35:39 2020

@author: yujin.wang
"""
import numpy as np

def  matrixSplit(A):
    dia = A.diagonal()
    D = dia*(np.eye(len(dia)))
    # print (dia,np.eye(len(dia)))
    E = A - D
    return D,E

def matrixB_z(D,E,b):
    return -np.linalg.inv(D).dot(E),np.linalg.inv(D).dot(b)

def JacobianIteration(B,z,x,tol=1e-3):
    err =1e9
    xi = x
    i = 1
    while err > tol:
        
        xi = B.dot(x) + z
        err = max(xi - x)
        print ('Run',i)
        print ('x:',xi)
        print ('Err:%5.3f' %(err))
        x = xi
        i+=1
    return x
        
if __name__ == '__main__':
    A = np.array([[3,2],[2,6]])
    b = np.array([2,-8])
    D,E = matrixSplit(A)
    B,z = matrixB_z(D,E,b)
    x0 = [0,0]
    ####################################
    JacobianIteration(B,z,x0,tol=1e-3)
    