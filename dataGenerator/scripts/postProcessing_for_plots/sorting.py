#!/usr/bin/env python

import numpy as np
import sys

family = np.array(["1-1-4-1-1",
 "1-1-4-1-2",
 "1-1-4-2-1",
 "1-1-4-2-2",
 "1-2-4-1-1",
 "1-2-4-1-2",
 "1-2-4-2-1",
 "1-2-4-2-2",
 "2-1-4-1-1",
 "2-1-4-1-2",
 "2-1-4-2-1",
 "2-1-4-2-2",
 "2-2-4-1-1",
 "2-2-4-1-2",
 "2-2-4-2-1",
 "2-2-4-2-2",
 ])
file = sys.argv[1]
data = np.loadtxt(file)

Hong = data[:,0]
slope = data[:,1]
cg = data[:,2]
print Hong
print "hong mei"
index = np.argsort(Hong)
print index+1
print 
print family[index]
