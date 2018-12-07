import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfits as pf
import copy
import pickle
import os
import sys
import pdb
import glob
import gzip
#import pymultinest
import os, threading, subprocess
import matplotlib.pyplot as plt
import json
import oifits
#import oifits_2015 as oifits
#import oifits_old as oifits
import time
from scipy.io.idl import readsav

from cp_tools import *

'''------------------------------------------------------------------------
cpo.py - Python class for manipulating oifits format closure phase data.
------------------------------------------------------------------------'''

# =========================================================================

class cpo():
    ''' Class used to manipulate multiple closure phase datasets'''

    def __init__(self, oifits):
        # Default instantiation.

        # if the file is a complete (kpi + kpd) structure
        # additional data can be loaded.
        try:
            self.extract_from_oifits(oifits)
        except:
            print('Invalid file.')


    def extract_from_oifits(self,filename):
        '''Extract closure phase data from an oifits file.'''

        data = oifits.open(filename)
        self.name = ''

        self.ndata = len(data.t3)

        for j in data.wavelength:
            wavel = data.wavelength[j].eff_wave
            self.wavel = wavel
            
        self.nwavs = len(self.wavel)

        self.target = data.target[0].target

        t3data = []
        t3err = []
        self.u = np.zeros((self.ndata,3))
        self.v = np.zeros((self.ndata,3))

        for j, t3 in enumerate(data.t3):
            t3data.append(t3.t3phi)
            t3err.append(t3.t3phierr)
            self.u[j,:] = [t3.u1coord,t3.u2coord,-(t3.u1coord+t3.u2coord)]
            self.v[j,:] = [t3.v1coord,t3.v2coord,-(t3.v1coord+t3.v2coord)]
            
        self.t3data = np.array(t3data)
        self.t3err = np.array(t3err)
        
        # Also load the v2
        nv2=len(data.vis2)
        vis2data=[]
        vis2err=[]
        self.v2_u=np.zeros((nv2))
        self.v2_v=np.zeros((nv2))
        for j,v2 in enumerate(data.vis2):
            vis2data.append(v2.vis2data)
            vis2err.append(v2.vis2err)
            self.v2_u[j]=v2.ucoord
            self.v2_v[j]=v2.vcoord
        self.vis2data=np.array(vis2data)
        self.vis2err=np.array(vis2err)

class icpo():
    ''' Class used to manipulate multiple closure phase datasets.
    This one structures the data differently, so it can be used with a covariance
    matrix, or with projected closure phases'''

    def __init__(self, oifits):
        # Default instantiation.

        # if the file is a complete (kpi + kpd) structure
        # additional data can be loaded.
        try:
           self.extract_from_oifits(oifits)
        except:
            print('Invalid file.')


    def extract_from_oifits(self,filename):
        '''Extract closure phase data from an oifits file.'''

        data = oifits.open(filename)
        self.name = ''

        self.ndata = len(data.t3)

        for j in data.wavelength:
            wavel = data.wavelength[j].eff_wave
            self.wavel = wavel
            break
        self.nwavs = len(self.wavel)

        self.target = data.target[0].target

        t3data = []
        t3err = []
        self.u = np.zeros((self.ndata,3))
        self.v = np.zeros((self.ndata,3))

        for j, t3 in enumerate(data.t3):
            t3data.append(t3.t3phi)
            t3err.append(t3.t3phierr)
            self.u[j,:] = [t3.u1coord,t3.u2coord,-(t3.u1coord+t3.u2coord)]
            self.v[j,:] = [t3.v1coord,t3.v2coord,-(t3.v1coord+t3.v2coord)]
            
        self.t3data = np.array(t3data)
        self.t3err = np.array(t3err)
        
        # Also load the v2
        nv2=len(data.vis2)
        vis2data=[]
        vis2err=[]
        self.v2_u=np.zeros((nv2))
        self.v2_v=np.zeros((nv2))
        for j,v2 in enumerate(data.vis2):
            vis2data.append(v2.vis2data)
            vis2err.append(v2.vis2err)
            self.v2_u[j]=v2.ucoord
            self.v2_v[j]=v2.vcoord
        self.vis2data=np.array(vis2data)
        self.vis2err=np.array(vis2err)

# =========================================================================

class PolDiffObs():
    ''' Class used to manipulate polarimetric differential observables'''

    def __init__(self, filename):
        # Default instantiation.

        try:
           self.extract_from_idlvar(filename)
        except:
            print('Invalid files.')


    def extract_from_idlvar(self,filename):
        '''Extract the differential visibility and closure phase data from 
        the files output by the IDL masking pipeline.'''
            
        # Use scipy's idl reading function to get a dictionary of all the saved variables
        data = readsav(filename)
        
        # The differential visibility quantities
        self.u=data['u']
        self.v=data['v']
        self.vis=data['vis']
        self.vis_err=data['vis_err']

        # Closure phase quantities        
        self.bs_u=data['bs_u']
        self.bs_v=data['bs_v']
        self.cp=data['cp']
        self.cp_err=data['cp_err']
        
        # Differential phase quantities
        self.ph=data['ph']
        self.ph_err=data['ph_err']
        
        # Generic quantities (that arent implemented yet)
#        self.wavel=data['wavel']
        