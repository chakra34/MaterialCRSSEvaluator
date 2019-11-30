#!/usr/bin/env python

import numpy as np
import sys,os

name = sys.argv[1]
data1 = np.loadtxt('{}_0.01.txt'.format(name))
data2 = np.loadtxt('{}_0.1.txt'.format(name))
data3 = np.loadtxt('{}_0.2.txt'.format(name))
data4 = np.loadtxt('{}_0.3.txt'.format(name))
data5 = np.loadtxt('{}_0.4.txt'.format(name))

# CRSS_1  = np.log(data1[:,1])
# CRSS_2  = np.log(data2[:,1])
# CRSS_3  = np.log(data3[:,1])
# CRSS_4  = np.log(data4[:,1])
# CRSS_5  = np.log(data5[:,1])
# values = np.vstack((CRSS_1,CRSS_2,CRSS_3,CRSS_4,CRSS_5)).T
# print np.nanmean(np.nanstd(values,axis=1)) # old way

## new way of precision

data1_crss = data1[:,1]
data2_crss = data2[:,1]
data3_crss = data3[:,1]
data4_crss = data4[:,1]
data5_crss = data5[:,1]

allData = np.vstack(( data1_crss,data2_crss,data3_crss,data4_crss,data5_crss))
std = []

# finding SD of 1-2 for all lambda then 1-3 for all lambda and so on

pairs_1_2 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,0]) or np.isnan(allData[i,1]): pairs_1_2.append(np.nan)
  else: pairs_1_2.append(np.log(allData[i,0]/allData[i,1]))
pairs_1_2 = np.array(pairs_1_2)
pairs_1_2 = pairs_1_2[~np.isnan(pairs_1_2)]
if len(pairs_1_2) > 1: std.append(np.std(pairs_1_2))

pairs_1_3 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,0]) or np.isnan(allData[i,2]): pairs_1_3.append(np.nan)
  else: pairs_1_3.append(np.log(allData[i,0]/allData[i,2]))
pairs_1_3 = np.array(pairs_1_3)
pairs_1_3 = pairs_1_3[~np.isnan(pairs_1_3)]
if len(pairs_1_3) > 1: std.append(np.std(pairs_1_3))

pairs_1_4 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,0]) or np.isnan(allData[i,3]): pairs_1_4.append(np.nan)
  else: pairs_1_4.append(np.log(allData[i,0]/allData[i,3]))
pairs_1_4 = np.array(pairs_1_4)
pairs_1_4 = pairs_1_4[~np.isnan(pairs_1_4)]
if len(pairs_1_4) > 1: std.append(np.std(pairs_1_4))

pairs_2_3 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,1]) or np.isnan(allData[i,2]): pairs_2_3.append(np.nan)
  else: pairs_2_3.append(np.log(allData[i,1]/allData[i,2]))
pairs_2_3 = np.array(pairs_2_3)
pairs_2_3 = pairs_2_3[~np.isnan(pairs_2_3)]
if len(pairs_2_3) > 1: std.append(np.std(pairs_2_3))

pairs_2_4 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,1]) or np.isnan(allData[i,3]): pairs_2_4.append(np.nan)
  else: pairs_2_4.append(np.log(allData[i,1]/allData[i,3]))
pairs_2_4 = np.array(pairs_2_4)
pairs_2_4 = pairs_2_4[~np.isnan(pairs_2_4)]
if len(pairs_2_4) > 1: std.append(np.std(pairs_2_4))

pairs_3_4 = []
for i in xrange(allData.shape[1]):
  if np.isnan(allData[i,2]) or np.isnan(allData[i,3]): pairs_3_4.append(np.nan)
  else: pairs_3_4.append(np.log(allData[i,2]/allData[i,3]))
pairs_3_4 = np.array(pairs_3_4)
pairs_3_4 = pairs_3_4[~np.isnan(pairs_3_4)]
if len(pairs_3_4) > 1: std.append(np.std(pairs_3_4))


print '1 head'
print 'mean_std'
print np.mean(std) 

