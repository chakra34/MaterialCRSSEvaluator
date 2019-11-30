#!/usr/bin/env python

import numpy as np
import sys,os

name = sys.argv[1]
val = sys.argv[2]
target = val.split('-')
target = np.array(target,dtype=float)
target = np.delete(target,np.where(target == 4))
data1 = np.loadtxt('{}_0.01.txt'.format(name))
CRSS_1  = data1[:,1]
out = [abs(np.log(target/target[i]) - np.log(CRSS_1/CRSS_1[i])) for i in xrange(4)]
rms = np.nanmean([np.sqrt(np.nanmean(np.square(out[i]))) for i in xrange(4)])
#rms_error = np.sqrt(np.nanmean(np.square(np.log(target) - np.log(CRSS_1) ) ) )

print '1 head'
print 'rms_error'
#print rms_error
print rms

