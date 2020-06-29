from sympy import *
import numpy as np

beta_val = 0.25
gamma_val = 0.5

var(['dt','wc','beta','gamma'])
A = Matrix([[1.,0.,-beta*dt*dt],[0.,1.,-gamma*dt],[wc*wc,0.,1.]])
B = Matrix([[1.,dt,(0.5-beta)*dt*dt],[0.,1.,(1.-gamma)*dt],[0.,0.,0.]])
C = A.inv() * B

A_eigen = C.subs([(beta,beta_val),(gamma,gamma_val)])
res = A_eigen.eigenvals()

# -dt**2*wc**2/2 - dt*wc*sqrt((dt*wc - 2)*(dt*wc + 2))/2 + 1
# -dt**2*wc**2/2 + dt*wc*sqrt((dt*wc - 2)*(dt*wc + 2))/2 + 1

#谱半径小于1 Ax=B 矩阵求解过程收敛
ct_list = list(res.keys())


ct1 = abs(ct_list[0]) <1
# print (ct1.subs(dt,2./wc))
ct2 = abs(ct_list[1]) <1
# print (ct2.subs(dt,2./wc))
plot_implicit(ct1)
plot_implicit(ct2)
