# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 12:41:40 2018

@author: yujin.wang
"""
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import time
import sys

def lanczos(A,b,nmax):
	m = np.size(A,1)
	alpha = []
	beta = [0]
	qprev = np.zeros(m)
	q = b / np.linalg.norm(b)
	for n in range(nmax):
		v = np.dot(q,A) #*sci.linalg.inv(q)
		temp = np.dot(v,q.T)
		alpha.append(temp[0][0,0])
		v = v - np.dot(beta[-1],qprev)-np.dot(alpha[-1],q)
		beta.append(sci.linalg.norm(v))
		qprev = q
		q = v/beta[-1]
	beta = beta[1:-1]
	T = np.diag(alpha) + np.diag(beta,1) +np.diag(beta,-1)
	return T


def qreigen(A,num=100):
	m = np.size(A,1)
	p = np.eye(m)
	for i in range(num):
		v = np.diag(A)
		q,r = sci.linalg.qr(A)
		A = np.dot(r,q)
		p = np.dot(p,q)
		tol = max(np.diag(A))-max(v)
		print 'TOL:',tol,'Max. Eig:',max(v)
		if np.abs(tol) <1e-6:
			break
#	s = np.diag(np.diag(r))
	return p,np.diag(A)

if __name__ == '__main__':
	A = np.matrix([[5,1,3],[1,2,0],[3,0,11]])
	print np.linalg.eig(A)[0]
	print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	p,s = qreigen(A)
	print s
	print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
	b = np.matrix([1,1,1])
	tridiag = lanczos(A,b,23)
	p,s = qreigen(tridiag,num=200)
	print max(s),min(s)