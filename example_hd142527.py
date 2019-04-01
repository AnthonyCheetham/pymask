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
import os


# OIFITS file name
wdir = os.path.dirname(os.path.realpath(__file__))+'/'
oif='HD142527_SPHERE_IRDIS_July2017.oifits'

#Scale the errors to make up for the actual number of indep wavs and indep clps
# The SPHERE mask has 7 holes, and there are n_holes/3 times more closure phases than independent closure phases
err_scale=np.sqrt(7./3.) 

chi2_grid=False # Run a coarse chi2 grid over separation, positition angle and contrast ratio
detec_lims=False # Calculate the detection limits in contrast ratio as a function of separation
get_params=False # Run emcee (Affine-Invariant MCMC) on data
nest=False # Run nested sampling on data
estimate_extra_error=False # Calculate how much extra uncertainty needs to be added in quadrature to give chi2_reduced = 1

#################################
#Load the oifits
icpo=pymask.cpo(wdir+oif)

#################################
#Now use the pymask routines for fitting and detection limits
#Set the initial parameters and priors
ip=[77.0,118.0,70.]
sep_prior=[20,250] # mas
pa_prior=[0,360] # degrees
crat_prior=[5,500] # contrast ratio
extra_error=0. #Error to add in quadrature to the closure phase errors (deg)

print 'Using',extra_error,'degrees of extra error'

## FITTING ROUTINES

# A coarse chi2 grid over seps, pas, crats
if chi2_grid:
    data2=pymask.coarse_grid(cpo,nsep=40,nth=40,ncon=40,thmin=pa_prior[0],thmax=pa_prior[1],
                  smin=sep_prior[0],smax=sep_prior[1],cmin=crat_prior[0],cmax=crat_prior[1])

# Use MCMC to get the parameters ( you may need to do "sudo pip install emcee")
if get_params:
    hammer_data=pymask.hammer(icpo,ivar=ip,nwalcps=50,niters=500,model='constant',
                       sep_prior=sep_prior,pa_prior=pa_prior,crat_prior=crat_prior,
                       err_scale=err_scale,extra_error=extra_error,plot=True)

    # Example plot to show the chains in separation
    chain = hammer_data['chain']
    plt.figure('chain')
    plt.clf()
    plt.plot(chain[:,:,0].T)
    plt.xlabel('Iteration')
    plt.ylabel('Sep (mas)')

# DETECTION LIMITS
if detec_lims:   
    lims_data=pymask.detec_limits(icpo,threads=3,nsim=500,cmax=250,nsep=15,ncon=30,smax=75,
                                  nth=30,save=False,err_scale=err_scale,extra_error=extra_error)

if nest:
    paramlimits = [sep_prior[0],sep_prior[1],pa_prior[0],pa_prior[1],crat_prior[0],crat_prior[1]]
    nest_data =pymask.nest(cpo,paramlimits=paramlimits,eff=0.3,multi=True,err_scale=err_scale,
                extra_error=extra_error,plot=True,npoints=1000)

if estimate_extra_error:
    plt.figure(1)
    pymask.find_extra_error(ip,cpo,err_scale=err_scale)