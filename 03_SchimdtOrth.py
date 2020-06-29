# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:10:09 2020

@author: yujin.wang
"""
import numpy as np

a = np.array([1.,0.])
b = np.array([1.,1.])
c = a - np.dot(np.dot(a,b)/np.dot(b,b),b)
print (c)

c = a - np.dot(a,b)
print (c)
