#!/usr/bin/env python

import numpy as np
import sys

file = sys.argv[1]                # file with only data showTable -d
ref_avg = sys.argv[2]             # reference ratios
s = ref_avg.split('-')
ref_ratios = np.array([float(s[0]), float(s[1]), float(s[3]), float(s[4])])
data = np.loadtxt(file)
val = data[:,0]
ratio = data[:,1]

ref_mean = np.nansum(np.log(ref_ratios)) #ref_ratios[0]*ref_ratios[1]*ref_ratios[2]*ref_ratios[3]
act_mean = np.nansum(np.log(ratio))#ratio[0]*ratio[1]*ratio[2]*ratio[3]
scale = np.exp((ref_mean - act_mean)/4.)
newRatio = ratio*scale

print '1 head'
print 'family','ratio'
print val[0],newRatio[0]
print val[1],newRatio[1]
print val[2],newRatio[2]
print val[3],newRatio[3]
