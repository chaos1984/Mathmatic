# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 12:41:40 2018

@author: yujin.wang
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.setrecursionlimit(10000) # 10000 is an example, try with different values

def get_line(xn, yn):
    def line(x):
        index = -1
         
        #找出x所在的区间
        for i in range(1, len(xn)):
            if x <= xn[i]:
                index = i-1
                break
            else:
                i += 1
         
        if index == -1:
            return -100
         
        #插值
        result = (x-xn[index+1])*yn[index]/float((xn[index]-xn[index+1])) + (x-xn[index])*yn[index+1]/float((xn[index+1]-xn[index]))
         
        return result
    return line
################Force###############
def CenterDiff(t,endtime,m,u,delta,du0):
	fext = flin(t)
	fint = klin(u)*u
	ddu = (fext - fint)/m
	du = du0 + ddu * delta
	u = u + du*delta
	t+=delta
	t_list.append(t);un_list.append(u)
	if t+delta > endtime:
		return t_list,un_list
	else:
		print 't:%f    u:%f' %(t,u)
		return CenterDiff(t,endtime,m,u,delta,du)

def NewmarkSolver(t,endT,deltaT,TOL,u,du,ddu,k,c,m,iternum):
	t += deltaT
	a1 = m/(beta*deltaT*deltaT) + gamma *c/(beta*deltaT)
	a2 = m/(beta*deltaT) + c*(gamma/beta-1)
	a3 = m*(1./(2*beta)-1) + deltaT*c*(gamma/(2*beta)-1)

	k_ = klin(u) + a1
	f_ = flin(t) + a1*u + a2*du + a3*ddu
	un = f_/k_
	Rn = f_ - klin(un)*un - a1*un
	print '##############TIME: %f####################' %(t)
	print time.ctime()
	print 'Iter.:%f Rn:%f' %(iternum,Rn)
#Newton-Raphson
#	print 'Rn',Rn,f_,f_sn,a1,un
	while abs(Rn) > TOL:
		iternum += 1
		k_= k_ + a1
		dun = Rn/k_
		un = un + dun
		Rn = f_ - klin(un)*un - a1*un
		print 'Iter.:%f Rn:%f' %(iternum,Rn)
	dun = gamma*(un-u)/(beta*deltaT) + du*(1-gamma/beta) + deltaT*ddu*(1-gamma/(2*beta))
	ddun = (un-u)/(beta*deltaT*deltaT) - du/(beta*deltaT) - ddu*(1./(2*beta)-1)
	t_list.append(t)
	un_list.append(un)
	if t+deltaT > endT:
		return  t_list,un_list
	else:
		print t,un
		return NewmarkSolver(t,endT,deltaT,TOL,un,dun,ddun,k,c,m,iternum)

if __name__ == '__main__':
	global iternum
	deltaT = 0.05;endT=10;
	T_n=1.
	gamma=0.5;
	beta=0.25;
	#Initial condition
	m = 1.
	k = [100.,200.]
	dl = [0,10.]
	f = [0,3,3,3]
	ft = [0,1,10,100]
	klin = get_line(dl,k)
	flin = get_line(ft,f)
	c= 0.0
	
	u = 0;du = 0;f = 0
	ddu = (flin(0) - c*du - f)/m
	t0 = 0
	f_s = 0
	nn = 0
	print gamma,beta
	for TOL in [0.1,0.01,0.001]:
		for deltaT in [0.01]:
			for gamma in [0.5,0.6,0.7,0.8,0.9,1.]:
				try:
					nn += 1
			#		TOL = 1e-3
					beta=0.25;
					t_list=[0];un_list=[0];f_sn_list=[0];k_T_list = [0];iternum = 0
					nt,nu = NewmarkSolver(t0,endT,deltaT,TOL,u,du,ddu,k,c,m,iternum)
					nmlabel = 'NM r=%3.2f b=%3.2f' %(gamma,beta)
					plt.plot(nt,nu,label=nmlabel)
					
					t_list=[0];un_list=[0];f_sn_list=[0];k_T_list = [0];iternum = 0
					u0 = 0;du0=0;ddu0=0
					ct,cu = CenterDiff(t0,endT,m,u0,deltaT,du0)
					plt.plot(ct,cu,label='cd')
					
					print u,du,ddu
					t_list=[0];un_list=[0];f_sn_list=[0];k_T_list = [0];iternum = 0
					beta=1./6;
					nt,nu = NewmarkSolver(t0,endT,deltaT,TOL,u,du,ddu,k,c,m,iternum)
					nmlabel = 'NM r=%3.2f b=%3.2f' %(gamma,beta)
					plt.plot(nt,nu,label=nmlabel)
					
					plt.grid()
					plt.legend()
					filename = '%d_TOL_%4.3f_deltaT_%4.3f_c_%3.2f.png' %(nn,TOL,deltaT,c)
					plt.title(filename)
					plt.xlabel('Time')
					plt.ylabel('Displacement')
					plt.savefig(filename,dpi=100)
					plt.close()
				except:
					pass