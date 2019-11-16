#!/usr/bin/env python

###########################################################################
# CRSSextractor.py: Statistic machinary for extracting CRSS ratio         #
#=========================================================================#
# Authors                                                                 #
# Aritra Chakraborty                                                      #
# Chen Zhang                                                              #
# Philip Eisenlohr                                                        #
# DESCRIPTION                                                             #
# -----------  Based on Hongmei Li's technique of identification of       #
#  CRSS parameters based on observed and theoretical Schmid factors       #
###########################################################################

import os,re,sys
import math                                                                                         # noqa
import numpy as np
from scipy import interpolate
from scipy import stats
from optparse import OptionParser
import damask
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

#---------------------------------- functions ---------------------------------------------------- #

def DprobabilityDensity(steps, cumProb,axis=1,range=None):
  """
  calculates probability density at given number of steps based on cumulative probability
  stored as numpy array of [value,cumP] for N pairs.
  """
  if len(cumProb) == 1: return cumProb
  else:
    if range is None: range = [min(cumProb[:,axis]),max(cumProb[:,axis])]
    if steps > len(cumProb): steps = len(cumProb) - 1 
    
    probDensity = []
    index = np.searchsorted(cumProb[:,axis],np.linspace(min(range),max(range),steps+1,endpoint=True))     # array indices where to evaluate
    index = np.maximum(np.zeros_like(index),np.minimum( (len(cumProb)-1) * np.ones_like(index),index))

    removed_points = []
    for i in xrange(len(index)-1):
      if list(cumProb[index[i]]) not in removed_points:
        first = cumProb[index[i]]
      for point in cumProb[index[i+1]:]:
        if first[0] == point[0]:
          removed_points.append([list(same) for same in cumProb[np.where(cumProb[:,0] == point[0])]])
        else:
          second = point
          break
      if first[0] == second[0]:
        print "only a single point is unique; derivative is that point"
        probDensity.append(first[0])
        probDensity.append(first[1])
      else:
        x = 0.5*(first[0] + second[0])
        y = (second[1] - first[1])/-(second[0] - first[0])
        probDensity.append(x)
        probDensity.append(y)
    probDensity = np.array(probDensity).reshape(len(probDensity)/2,2)

    plt.plot(cumProb[:,0],cumProb[:,1],c='r',label='actual')
    plt.plot(probDensity[:,0],probDensity[:,1],c='g',marker='o',label='derivative')
  #  plt.xscale('log')
    plt.yscale('log')
  #  plt.ylim(-0.1,0.1)
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()
    return probDensity


def NewProbabilityDensity(cumProb_obs,cumProb_pop,family=None):
  """ new way of calculating the derivative based on a dynamic delta"""
  if len(cumProb_obs) <= 3:
    return np.array([ [np.nan,np.nan] ])             # not enough data point
  else:
    x_obs = cumProb_obs[:,0]
    lowerlimit = cumProb_obs[np.where(cumProb_obs[:,0] >= np.log(0.1))[-1][-1],0]
    x_min,x_max = max(np.amin(cumProb_obs[:,0]),lowerlimit),np.amax(cumProb_obs[:,0])
    index = np.where(cumProb_pop[:,0] >= x_min)[0]
    pop   = cumProb_pop[index]
  #   elif len(cumProb_obs) <= 3:
  #     x_axis = x_obs
  #     y_obs  = cumProb_obs[:,1]
    delta = 6 if len(cumProb_obs) > 6 else len(cumProb_obs) - 1
    x_axis = np.linspace(x_min,x_max,delta)
    y_obs = [cumProb_obs[np.where(cumProb_obs[:,0] == x_axis[0])[0],1][0] ]
    for i in xrange(1,len(x_axis)-1):
      x1 = x_axis[i]
      if len(np.where(cumProb_obs[:,0] == x1)[0]) == 1:
        y1 = cumProb_obs[np.where(cumProb_obs[:,0] == x1)[0],1]
      else:
        index_left  = np.where(cumProb_obs[:,0] < x1)[-1][0]
        index_right = np.where(cumProb_obs[:,0] > x1)[-1][-1]
        local_slope = (cumProb_obs[index_right,1] - cumProb_obs[index_left,1])/(cumProb_obs[index_right,0] - cumProb_obs[index_left,0])
        y1 = cumProb_obs[index_left,1] + (x1 - cumProb_obs[index_left,0])*local_slope
      y_obs.append(y1)
    y_obs.append(cumProb_obs[np.where(cumProb_obs[:,0] == x_axis[-1])[0],1][0])
    y_pop = []
    for i in xrange(len(x_axis)-1):
      x1 = x_axis[i]
      if len(np.where(cumProb_pop[:,0] == x1)[0]) == 1:
        y1 = cumProb_pop[np.where(cumProb_pop[:,0] == x1)[0],1][0]
      else:
        index_left  = np.where(cumProb_pop[:,0] < x1)[-1][0]
        index_right = np.where(cumProb_pop[:,0] > x1)[-1][-1]
        local_slope = (cumProb_pop[index_right,1] - cumProb_pop[index_left,1])/(cumProb_pop[index_right,0] - cumProb_pop[index_left,0])
        y1 = cumProb_pop[index_left,1] + (x1 - cumProb_pop[index_left,0])*local_slope
      y_pop.append(y1)
    y_pop.append(cumProb_pop[np.where(cumProb_obs[:,0] == x_axis[-1])[0],1][0])
    plt.plot(x_axis,y_obs,c='r',marker='o',label='inter_obs',linestyle='--')
    plt.plot(cumProb_obs[:,0],cumProb_obs[:,1],marker='^',c='g',label='raw_obs',linestyle='--')
    plt.plot(x_axis,y_pop,c='r',marker='o',label='inter_pop')
    plt.plot(cumProb_pop[:,0],cumProb_pop[:,1],c='g',label='raw_pop')
  #   plt.xscale('log')
  #   plt.yscale('log')
    plt.xlim(-2.0,-0.6)
    plt.ylim(0.0,0.6)
    plt.grid(True)
    plt.legend(loc='best')
  #  plt.savefig('{}_interpolated.pdf'.format(family))
    plt.close()

    first_slope_obs = - (y_obs[1] - y_obs[0])/(x_axis[1] - x_axis[0])
    first_slope_pop = - (y_pop[1] - y_pop[0])/(x_axis[1] - x_axis[0])
    x = (x_axis[1] + x_axis[0])*0.5
    ratio = first_slope_obs/first_slope_pop
    activity = [x]
    activity.append(ratio)
    for i in xrange(1,len(x_axis)-1):
      slope_obs = - (y_obs[i+1] - y_obs[i-1])/(x_axis[i+1] - x_axis[i-1])
      slope_pop = - (y_pop[i+1] - y_pop[i-1])/(x_axis[i+1] - x_axis[i-1])
      ratio = slope_obs/slope_pop
      activity.append(x_axis[i])
      activity.append(ratio)

    last_slope_obs = - (y_obs[-1] - y_obs[-2])/(x_axis[-1] - x_axis[-2])
    last_slope_pop = - (y_pop[-1] - y_pop[-2])/(x_axis[-1] - x_axis[-2])
    x = (x_axis[-1] + x_axis[-2])*0.5
    ratio = last_slope_obs/last_slope_pop
    activity.append(x)
    activity.append(ratio)
    activity = np.array(activity).reshape(len(activity)/2,2)
    return activity


def DataToSchmid(data):
  """
  takes the input data and gives a list of arrays corresponding to different families
  the first column in each array is SF, second column is weights for each family
  the input data consists of two columns 'Family' and 'SF'.
  """

  sorted_data = data[data[:,0].argsort()]
  all_families  = [1,2,3,4]
  families_seen = list(np.unique(sorted_data[:,0]))
  first_point = []
  for family in families_seen:                                                       # added the intrinsic first point of SF 0.5 (only to the observed families)
    first_point.append(family)
    first_point.append(0.5)
  first_point = np.array(first_point).reshape(len(first_point)/2,2)
  added_data  = np.vstack((data,first_point))
  sorted_data = added_data[added_data[:,0].argsort()]
  if len(families_seen) != 4:
    absent_families = list(set(all_families) - set(families_seen))
    artificial = []
    for family in absent_families:
      artificial.append(family)
      artificial.append(np.nan)
    artificial  = np.array(artificial).reshape(len(artificial)/2,2)
    new_data = np.vstack((added_data,artificial))
    sorted_data = new_data[new_data[:,0].argsort()]
    sorted_SF   = sorted_data[:,1]
    weights     = np.ones_like(sorted_SF)
    sorted_SF   = np.vstack((sorted_SF,weights,weights)).T
    family, members = np.unique(new_data[:,0], return_counts=True)
    SF_list = []
    count = 0
    for i in xrange(len(members)):
      SF_list.append(sorted_SF[count:count+members[i]])
      count += members[i]
  else:
    sorted_SF   = sorted_data[:,1]
    weights     = np.ones_like(sorted_SF)
    sorted_SF   = np.vstack((sorted_SF,weights,weights)).T
    family, members = np.unique(added_data[:,0], return_counts=True)
    SF_list = []
    count = 0
    for i in xrange(len(members)):
      SF_list.append(sorted_SF[count:count+members[i]])
      count += members[i]
  
  return SF_list


def cdf(SF_list,observations = None):                                                     # distribution function
  """
  takes the input data comprising of list of arrays with SF and weights and outputs 
  a list of arrays with SF, weights and cummulative probabilities
  for observation the weights are scaled based on frequency of observations
  e.g. for a family having total 10 and observed 3, the weights are 8,9,10
  """

  number = 0
#   if observations:
#     for i in xrange(len(SF_list)):
#       for j in xrange(len(SF_list[i])):        
#         SF_list[i][j,1] = observations[i] - len(SF_list[i]) + 1
  for i in xrange(len(SF_list)):
    SF_list[i] = SF_list[i][SF_list[i][:,0].argsort()][::-1]
    s = SF_list[i][0,1]
    for j in xrange(len(SF_list[i])):
#       if observations: SF_list[i][j,2] = s/np.sum(SF_list[i][:,1])  #np.sum(observations[i])
#       else :           SF_list[i][j,2] = s/np.sum(SF_list[i][:,1])
        SF_list[i][j,2] = s/np.sum(SF_list[i][:,1])  #np.sum(observations[i])
        s += 1
    SF_list[i][:,2] *= float(len(SF_list[i]))/observations[i]

  return SF_list


def activation_ratio(population_SF,observation_SF,plot=False):            # activation ratio based on observed Schmid Factor
  """
  activationRatio based on probability density of population and observation
  the population probability density are interpolated for observed ones to match SF's
  output is a list of arrays having the corresponding SF's and activationRatio per family
  """

  density_population  = []
  density_observation = []
  activationRatio           = []
  interpolatedDensity_population = []
  x = []
  c = 0

  for i in xrange(len(observation_SF)):
    if i == 3 or i == 2: delta = 4
    else: delta = 4
    density_observation.append(DprobabilityDensity(delta,cumProb_obs,axis=1,range=None,))
    density_observation[i][:,0] = np.exp(density_observation[i][:,0])
    density_population.append(DprobabilityDensity(delta,cumProb_pop,axis=1,range=None,))
    density_population[i][:,0] = np.exp(density_population[i][:,0])  
    x.append(density_observation[i][:,0])

    interpolation_population_density = interpolate.UnivariateSpline(density_population[i][:,0],
                                                                    density_population[i][:,1],
                                                                    k=1,ext='extrapolate') 
    interpolatedDensity_population.append(interpolation_population_density(x[-1]))
      
    activationRatio.append(np.vstack((x[-1],
                                     (density_observation[i][:,1]/interpolatedDensity_population[i]))).T)
    if plot == True:
      legend = ['basal','prism','pyr<a>','pyr<c+a>']
      color  = ['b'    ,'r'    ,'g'     ,'y',      ]
      plt.title('Activation ratio vs SF')
      plt.xlim(0.5,0.1)
      plt.ylim(1e-3,1e3)
      plt.grid(True)
      plt.loglog(activationRatio[i][:,0],activationRatio[i][:,1],c=color[c],linewidth=1.2,marker='+',label=legend[c])
      c += 1
  if plot == True:
    plt.legend(loc='best')
    plt.show()
  return activationRatio


def NewActivation_Ratio(population_SF,observation_SF,plot=False):            # activation ratio based on observed Schmid Factor
  """
  New method
  """
  
  activationRatio           = []
  legend = ['basal','prism','pyr<a>','pyr<c+a>']
  color  = ['b'    ,'r'    ,'g'     ,'y',      ]
  for i in xrange(len(observation_SF)):
    cumProb_obs = np.vstack((np.log(observation_SF[i][:,0]),observation_SF[i][:,2])).T
    cumProb_pop = np.vstack((np.log(population_SF[i][:,0]),population_SF[i][:,2])).T
    activity = NewProbabilityDensity(cumProb_obs,cumProb_pop,family=legend[i])
    activity[:,0] = np.exp(activity[:,0])
    activationRatio.append(activity)

  return activationRatio

def ScaledDistribution(SF_distributions,shift=None):                                                     # distribution function
  """
  takes the input data comprising of list of arrays with SF and weights and outputs
  scales the cum probability based on respective observations
  """
  total_number = sum([len(i) for i in SF_distributions])
  present = 0.0
  for dist in SF_distributions:
    dist[:,2] = dist[:,2]*len(dist)/total_number
    if shift == True:
      present += len(dist)
      dist[:,2] = dist[:,2] + float( (total_number - present)/total_number)
  return SF_distributions

def function(x,m,c):
  return m*x + c

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
CRSS extractor
""", version = scriptID)

parser.add_option('-o','--obs',
                  dest = 'obs_file',
                  help = 'observation file name')

parser.add_option('-l','--label',
                  dest = 'labels',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing slip family [%default]')

parser.add_option('-S','--SF',
                  dest = 'SF',
                  type = 'string', metavar = 'string',
                  help = 'heading of columns containing Schmid Factors [%default]')

parser.add_option('-p', '--plot',
                  dest = 'plots',
                  action = 'store_true',
                  help = 'whether to output the plots [%default]')

parser.add_option('-e', '--exp',
                  dest = 'expectation',
                  action = 'store_true',
                  help = 'whether to scale population distribution [%default]')

parser.add_option('--slope',
                  dest = 'slope',
                  action = 'store_true',
                  help = 'whether output crss is ratio of activity slope [%default]')

parser.add_option('--cg',
                  dest = 'cg',
                  action = 'store_true',
                  help = 'whether output crss is ratio of mid point of activity curves [%default]')

parser.add_option('--hcg',
                  dest = 'horicg',
                  action = 'store_true',
                  help = 'whether output crss is ratio of mid point of activity curves [%default]')

parser.add_option('--out',
                  dest = 'output',
                  action = 'store_true',
                  help = 'whether output obs, pop, and activity data [%default]')

parser.add_option( '--lambda',
                  dest = 'l',
                  type = 'string',
                  help = 'lambda [%default]')


parser.set_defaults(labels = "Family",
                    SF  = "SF",
                    expectation = False,
                    plots = False,
                    l     = None,
                    slope = False,
                    cg    = False,
                    corrcoef = False,
                    output = False,
                    horicg = False,
                   )

(options,filenames) = parser.parse_args()
legend = ['basal','prism','pyr<a>','pyr<c+a>']
color  = ['b'    ,'r'    ,'g'     ,'y',      ]

n = 3.0                               # variable
out_dir = os.path.dirname(os.path.realpath(options.obs_file))                      # files are generated in the directory of the observed file
if options.obs_file:
#------------------- Reading the simulated tip from file -----------------

  table_obs = damask.ASCIItable(name = options.obs_file, buffered = False,readonly=True)
  table_obs.head_read()
  table_obs.data_readArray()
  table_obs.data = np.vstack((table_obs.data[:,table_obs.label_index(options.labels)],
                              table_obs.data[:,table_obs.label_index(options.SF)]) ).T
  SF_observation = DataToSchmid(table_obs.data)

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              readonly = True,
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()
  table.data_readArray()
  table.data = np.vstack((table.data[:,table.label_index(options.labels)].T,
                          table.data[:,table.label_index(options.SF)].T ) ).T


# --------------------------- Finding number of members in a family (population count) -------------- #


  SF_population = DataToSchmid(table.data)
  observations  = [len(i) for i in SF_population]
  theoretical_SF_distribution = cdf(SF_population,observations)  
  observed_SF_distribution = cdf(SF_observation,observations)
  print [len(j) for j in SF_observation]
  inter = []

# --------------------------- sanity check for interpolation method --------------------------------- #

  if options.output:
    for i in xrange(len(observed_SF_distribution)):
      with open('{}/SF_observationDistribution_sys{}.txt'.format(out_dir,i+1),'w') as observed_SF:
        observed_SF.write('1 head\n')
        observed_SF.write('SF cumProb\n')
        for j in xrange(len(observed_SF_distribution[i])):
            observed_SF.write('{} {}\n'.format(observed_SF_distribution[i][j,0],
                                               observed_SF_distribution[i][j,2]))      
    for i in xrange(len(theoretical_SF_distribution)):
      with open('{}/SF_theoreticalDistribution_sys{}.txt'.format(out_dir,i+1),'w') as observed_SF:
        observed_SF.write('1 head\n')
        observed_SF.write('SF cumProb\n')
        for j in xrange(len(theoretical_SF_distribution[i])):
            observed_SF.write('{} {}\n'.format(theoretical_SF_distribution[i][j,0],
                                               theoretical_SF_distribution[i][j,2]))
  delta = 6
# ----------------- slope of the distribution  Population  ---------------------------------------- #

  activationRatio = NewActivation_Ratio(theoretical_SF_distribution, observed_SF_distribution,plot=options.plots)

  if options.output:
    for i in xrange(len(activationRatio)):
      with open('{}/ActivationRatioDistribution_sys{}.txt'.format(out_dir,i+1),'w') as observed_SF:
        observed_SF.write('1 head\n')
        observed_SF.write('SF Activity\n')
        for j in xrange(len(activationRatio[i])):
            observed_SF.write('{} {}\n'.format(activationRatio[i][j,0],
                                               activationRatio[i][j,1]))

    for i in xrange(len(activationRatio)):
      if len(activationRatio[i][~np.isnan(activationRatio[i][:,1])]) > 1:
        plt.loglog(activationRatio[i][:,0],activationRatio[i][:,1],c=color[i],linewidth=1.2,marker='o',label=legend[i])
    
    plt.title('Activation ratio vs SF')
    plt.xlim(0.1,0.5)
    plt.ylim(1e-3,1e1)
    plt.grid(True)
    plt.legend(loc=2)
#    plt.savefig('{}_{}_Activity.pdf'.format(options.obs_file.split('*.txt')[0],options.l))
  #  plt.show()
    plt.close()


#----------------------------------------------------------------------------------------------------- #
  
  ratio = []
  valid_ratio = np.nan
  if options.slope:
    slopes = []
    valid_slope = np.nan
    for i in xrange(len(activationRatio)):
      valid_activationRatio = activationRatio[i][~np.isnan(activationRatio[i][:,1])]  
      if len(valid_activationRatio) > 0:
        x = np.log(valid_activationRatio[:,0])
        y = np.log(valid_activationRatio[:,1])
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        valid_slope = slope
        print 'R2',r_value**2
        x_a = np.log(np.linspace(0.1,0.5,5))
        y_fitted = np.exp(slope*x_a + intercept)
        with open('{}/Slope_ActivationRatioDistribution_sys{}.txt'.format(out_dir,i+1),'w') as activity_slope:
          activity_slope.write('1 head\n')
          activity_slope.write('SF slopeFitted\n')
          for j in xrange(len(y_fitted)):
              activity_slope.write('{} {}\n'.format(np.exp(x_a[j]),y_fitted[j]))
#         plt.scatter(np.exp(x),np.exp(y),marker='^',c=color[i])
#         plt.plot(np.exp(x),np.exp(slope*x + intercept),c=color[i])
#         plt.grid(True)
#         plt.xscale('log')
#         plt.yscale('log')
#         plt.show()
        slopes.append(slope)
      else: 
        slopes.append(np.nan)
        with open('{}/Slope_ActivationRatioDistribution_sys{}.txt'.format(out_dir,i+1),'w') as activity_slope:
          activity_slope.write('1 head\n')
          activity_slope.write('SF slopeFitted\n')
          for j in xrange(1):
              activity_slope.write('{} {}\n'.format(np.nan,np.nan))

    if np.isnan(slopes[0]) == True or slopes[0] <= 0.0: slopes = (slopes/valid_slope)
    else: slopes = slopes/slopes[0]
    ratio = slopes

  elif options.cg:
    slopes = []
    valid_slope = np.nan
    min_SF = np.nanmax([np.amin(i[:,0]) for i in activationRatio])
    max_SF = np.nanmin([np.amax(i[:,0]) for i in activationRatio])
    spread = max_SF - min_SF
    print min_SF,max_SF
    for i in xrange(len(activationRatio)):
      valid_activationRatio = activationRatio[i][~np.isnan(activationRatio[i][:,1])]  
      if len(valid_activationRatio) > 0:
        if spread <= 0.1:
          slopes.append(1/np.average(valid_activationRatio[:,1]))
        else:
          index = np.where((valid_activationRatio[:,0] >= min_SF) & (valid_activationRatio[:,0] <= max_SF) )[0]
          region = valid_activationRatio[index]
          slopes.append(1/np.average(region[:,1]))
        valid_slope = slopes[-1]
      else: slopes.append(np.nan)
    ratio   = slopes/valid_slope

  elif options.horicg:
    slopes = []
    valid_slope = np.nan
    for i in xrange(len(activationRatio)):
      valid_activationRatio = activationRatio[i][~np.isnan(activationRatio[i][:,0])]  
      if len(valid_activationRatio) > 0:
        slopes.append(np.average(activationRatio[i][:,0]))
        valid_slope = slopes[-1]
      else: slopes.append(np.nan)
    ratio   = slopes/valid_slope

  else:
    SF = np.nanmin(np.array([np.amax(activationRatio[0][:,0]),np.amax(activationRatio[1][:,0]),
             np.amax(activationRatio[2][:,0]),np.amax(activationRatio[3][:,0]) ]) )
    for i in xrange(len(activationRatio)):
      valid_activationRatio = activationRatio[i][~np.isnan(activationRatio[i][:,1])]
      if (len(valid_activationRatio) > 0) :                                                     # In case no observation "nan"
        index = np.where(activationRatio[i][:,0] == SF)[0]
        if len(index) == 1:
          ratio.append(activationRatio[i][index,1][0])
        else:
          if len(np.where(activationRatio[i][:,0] < SF)[0]) >= 1 and len(np.where(activationRatio[i][:,0] > SF)[0]) >= 1:
            index_left  = np.where(activationRatio[i][:,0] < SF)[0][-1]          
            index_right = np.where(activationRatio[i][:,0] > SF)[0][0]
            local_slope = (activationRatio[i][index_right,1] - activationRatio[i][index_left,1])/(activationRatio[i][index_right,0] - activationRatio[i][index_left,0])
            slope = activationRatio[i][index_left,1] + local_slope * (SF - activationRatio[i][index_left,0])
          else: slope = activationRatio[i][-1,1]
          ratio.append(slope)
        diff = ratio[-1]
        valid_ratio = diff
      else: 
        diff = np.nan
        ratio.append(diff)                                                         # appending for last point
    if np.isnan(ratio[0]) == True or ratio[0] <= 0.0: ratio = 1/(ratio/valid_ratio)
    else: ratio = 1/(ratio/ratio[0])# **(1/float(n))


  if options.l:
    with open('{}/CRSS_results_{}.txt'.format(out_dir,options.l),'w') as crss:
      crss.write("1 head\n")
      crss.write("family ratio\n")
      left = np.ones_like(ratio) * -1
      for i in xrange(len(ratio)):
        if ratio[i] == np.nan: ratio[i] = -1
        left[i] = ratio[i]
        crss.write("{} {} {} {} {} {}\n".format(i+1,ratio[i],left[0],left[1],left[2],left[3]))

  damask.util.croak("ratio {}".format(ratio))
  table.close()