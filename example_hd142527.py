# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 10:57:15 2013

@author: anthony
Run indep_clps calibration and fitting on data
"""

import matplotlib.pyplot as plt
import numpy as np
import pymask
import pickle
import scipy


#oifits directory
wdir1='/home/anthony/data/gpi_data/'
wdir=wdir1+'gpi_May14/Analysis/hd142527_J/oifits/'

# OIFITS file name
oif='HD142527_oneblock_C1_python.oifits'

#Scale the errors to make up for the actual number of indep wavs and indep clps
err_scale=np.sqrt(37./17)*np.sqrt(10./3.)

F=False

detec_lims=F
chi2_grid=F
get_params=F

#################################
#Load the oifits
icpo=pymask.cpo(wdir+oif)

#################################
#Now use the pymask routines for fitting and detection limits
#Set the initial parameters and priors
ip=[77.0,118.0,70.]
sep_prior=[50,100]
pa_prior=[90,150]
crat_prior=[20,150]
extra_error=0. #Error to add in quadrature to the closure phase errors (deg)

print 'Using',extra_error,'degrees of extra error'

## FITTING ROUTINES

plt.close('all')

# A coarse chi2 grid over seps, pas, crats
if chi2_grid:
    data2=pymask.coarse_grid(icpo,cmax=150)

# Use MCMC to get the parameters ( you may need to do "sudo pip install emcee")
if get_params:
    data=pymask.hammer(icpo,ivar=ip,plot=True,niters=500,model='constant',nwalcps=10,
                       sep_prior=sep_prior,pa_prior=pa_prior,crat_prior=crat_prior,
                       err_scale=err_scale)
# DETECTION LIMITS
if detec_lims:   
    lims_data=pymask.detec_limits(icpo,threads=3,nsim=500,cmax=250,nsep=15,ncon=30,smax=75,
                                  nth=30,save=False,errscale=err_scale)
