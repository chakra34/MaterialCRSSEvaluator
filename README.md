# MaterialCRSSEvaluator
Evaluates Critical Resolved Shear Stresses (CRSS) from slip trace data

###########################################################################
# Machinery to evaluate material's CRSS values from slip trace data       #
#                                                                         #
#=========================================================================#
# Authors                                                                 #
# Aritra Chakraborty                                                      #
# Chen Zhang                                                              #
# Philip Eisenlohr                                                        #
# DESCRIPTION                                                             #
#  Reading manual                                                         #
###########################################################################



#--- Prerequisites --#

Python 2.7.*          # The code is written in Python version 2.7

#------------------------#

Usage

"CRSSextractor.py"

This code takes a slip trace data consisting of "theoretical population of traces" (i.e., total possible slip traces of all the grains in the observed patch) as well as the "actual observed traces" (slip traces observed from the slip trace analysis after small deformation) and calculates the critical resolved shear stress (CRSS) values based on three different methodologies.

By default it uses the methodology described in 

Li, H., D. E. Mason, T. R. Bieler, C. J. Boehlert, and M. A. Crimp. "Methodology for estimating the critical resolved shear stress ratios of α-phase Ti using EBSD-based trace analysis." Acta Materialia 61, no. 20 (2013): 7555-7567.

Additionally, it has two other means to get the CRSS, one using the slope of the observed activity ; and another using a numerial average of the activity distributions.

In this code an exemplary data set with a pre-set CRSS of 1-1-1-1 ( basal, prism, pyr<a>, and pyr<c+a>) is included in the "data" directory that was generated artifically using Crystal Plasticity Fast Fourier Transform (CP-FFT) simulations, using the open source package DAMASK (https://damask.mpie.de/) from which some of the utilities of file parsing are used in the code as well.

syntax:

go into the code directory and execute the following as a test:

./CRSSextractor.py ../data/Population_1-1-4-1-1.txt  -o ../data/Observation.txt

The resulting ratio is:

ratio 1.         1.04075518 1.22855641 1.05571215

which is very close to actual 1.,1.,1.,1.

There are other flags that may or may not be used, and the code may be modified as well based on the user's preference.

Note: There are some interpolations (Univariate spline) that are done which might be tuned based on the density of a data set as well as the way the numerical differentiation is done to calculate the probability density function from the cumulative distribution function. 
The way the code is now has been optimized for most of the hexagonal data set on which it was tested for cp-Titanium.

Fuirther details about the methodology and citation for the present code is available in the following journal:

Chakraborty, Aritra, Chen Zhang, Shanoob Balachandran, Thomas R. Bieler, and Philip Eisenlohr. "Assessment of surface and bulk-dominated methodologies to measure critical resolved shear stresses in hexagonal materials." Acta Materialia (2019).





