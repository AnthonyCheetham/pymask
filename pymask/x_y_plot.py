# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 13:17:03 2014

@author: anthony
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool
from pymask import cp_model
from matplotlib.cm import register_cmap, cmap_d
 # =========================================================================
    
def chi2_grid(everything):
    '''Function for multiprocessing, does 2d chi2 grid for xy_grid'''
    cpo=everything['cpo']
    data_cp=cpo.t3data
    chi2=np.zeros((len(everything['ys']),len(everything['cons'])))
    x=everything['x']
    ys=everything['ys']
    seps=np.sqrt(x**2+ys**2)
    pas=np.angle(np.complex(0,1)*ys + np.complex(1,0)*x,True)
    for ix in range(ys.size):
        for k,con in enumerate(everything['cons']):
            mod_cps = cp_model([seps[ix],pas[ix],con],cpo.u,cpo.v,cpo.wavel)
            chi2[ix,k]=np.sum(((data_cp.ravel()-mod_cps.ravel())/cpo.t3err.ravel())**2)       
    return chi2

# =========================================================================

 # =========================================================================

def xy_grid(cpo,nxy=30,ncon=32,xymax='Default', cmin=10.,cmax=500.,
             threads=0,err_scale=1.,extra_error=0.,plot_chi2=False,fix_crat=False,
             cmap='ds9cool'):

    '''An attempt to copy Sylvestre's chi2 grid plots, using x and y instead
    of separation and position angle.

    Written by A Cheetham, with some parts stolen from other pysco/pymask routines.'''
    
    ds9cool = {'red': lambda v : np.interp(v, [0, 0.29, 0.76, 1], [0, 0, .1, 1]),
            'green': lambda v : np.interp(v, [0, 0.22, 0.96, 1], [0, 0, 1, 1]),
            'blue': lambda v : np.interp(v, [0, 0.53, 1], [0, 1,1])}
    
    register_cmap('ds9cool', data=ds9cool)
    #------------------------
    # first, load your data!
    #------------------------

    ndata=cpo.ndata
    
    u,v = cpo.u,cpo.v
    
    cpo.t3err=np.sqrt(cpo.t3err**2+extra_error**2)
    cpo.t3err*=err_scale

    wavel = cpo.wavel

    w = np.array(np.sqrt(u**2 + v**2))

    if xymax == 'Default':
#        xymax = cpt.rad2mas(1./np.min(w/np.max(wavel)))
        xymax=(1./np.min(w/np.max(wavel)))/np.pi*(180*3600*1000)

    #------------------------
    # initialise grid params
    #------------------------

    xys= np.linspace(-xymax,xymax,nxy)
    cons = cmin  + (cmax-cmin)  * np.linspace(0,1,ncon)
    
    if fix_crat != False:
        cons=np.array([fix_crat])
        ncon=1
    
    #------------------------
    # Calculate chi2 at each point
    #------------------------

    tic = time.time() # start the clock
    chi2=np.zeros((nxy,nxy,ncon))
    if threads ==0:
        toc=time.time()
        for ix,x in enumerate(xys):
            everything={'x':x,'cons':cons,'ys':xys,'cpo':cpo,'ix':ix}
            chi2[ix,:,:]=chi2_grid(everything)
            if (ix % 50) ==0:
                tc=time.time()
                print('Done '+str(ix)+'. Time taken: '+str(tc-toc)+'seconds')
                toc=tc
    else:
        all_vars=[]
        for ix in range(nxy):
            everything={'x':xys[ix],'cons':cons,'ys':xys,'cpo':cpo,'ix':ix}
            all_vars.append(everything)
        pool = Pool(processes=threads)
        chi2=pool.map(chi2_grid,all_vars)
    tf = time.time()
    if tf-tic > 60:
        print('Total time elapsed: '+str((tf-tic)/60.)+'mins')
    elif tf-tic <= 60:
        print('Total time elapsed: '+str(tf-tic)+' seconds')
    chi2=np.array(chi2)
    best_ix=np.where(chi2 == np.amin(chi2))
    
    #hack: if the best chi2 is at more than one location, take the first.
    bestx=xys[best_ix[0][0]]
    besty=xys[best_ix[1][0]]
    sep=np.sqrt(bestx**2+besty**2)
    pa=np.angle(np.complex(bestx,besty),True)
    best_params=[sep,pa,cons[best_ix[2][0]]]
    best_params=np.array(np.array(best_params).ravel())
    print('Separation '+str(best_params[0])+' mas')
    print('Position angle '+str(best_params[1])+' deg')
    print('Contrast Ratio '+str(best_params[2]))
    # ---------------------------------------------------------------
    #                        sum over each variable so we can visualise it all
    # ---------------------------------------------------------------
    temp_chi2=ndata*chi2/np.amin(chi2)
    like=np.exp(-(temp_chi2-ndata)/2)
    x_y=np.sum(like,axis=2)

    # ---------------------------------------------------------------
    #                        contour plot!
    # ---------------------------------------------------------------
    plt.figure(0)
    plt.clf()
    if plot_chi2:
        im=np.min(chi2,axis=2)
        # convert to reduced chi2
#        im/= (cpo.t3data.size-3)
#        print 'reduced chi2'
    else:
        im=x_y
#    plt.imshow(im,extent=[np.amin(xys),np.amax(xys),np.amin(xys),np.amax(xys)],aspect='auto',cmap=cmap)
    # Plot it with RA on the X axis
    plt.imshow(im,extent=[np.amin(xys),np.amax(xys),np.amin(xys),np.amax(xys)],aspect='auto',cmap=cmap)
    plt.colorbar()
    plt.ylabel('Dec (mas)')
    plt.xlabel('RA (mas)')
#    plt.title('X vs Y')
    plt.title('Chi2 of single companion fit')
    
    plt.plot([0],[0],'wo')
    plt.xlim(xys[-1],xys[0])
    plt.ylim(xys[0],xys[-1])
    
    # ---------------------------------------------------------------
    #               And the detection limits that come for free!
    # ---------------------------------------------------------------
    chi2_null=np.sum((cpo.t3data/cpo.t3err)**2)
    # Define the detec limits to be the contrast at which chi2_binary - chi2_null < 25
    detecs=(chi2-chi2_null) < 25
    detec_lim=np.zeros((nxy,nxy))
    for x_ix in range(nxy):
        for y_ix in range(nxy):
            detectable_cons=cons[detecs[x_ix,y_ix,:]]
            if len(detectable_cons) ==0:
                detec_lim[x_ix,y_ix]=cons[-1]
            else:
                detec_lim[x_ix,y_ix]=np.min(detectable_cons)
                            
    plt.figure(1)
    plt.clf()
#    plt.imshow(detec_lim,extent=(xys[0],xys[-1],xys[0],xys[-1]),cmap=cmap)
    # Plot it with RA on the X axis
    plt.imshow(detec_lim,extent=(xys[0],xys[-1],xys[0],xys[-1]),cmap=cmap)
    plt.colorbar()
    plt.title('Detection limits')
    plt.xlabel('RA (mas)')
    plt.ylabel('Dec (mas)')
    plt.xlim(xys[-1],xys[0])
    plt.ylim(xys[0],xys[-1])
    
    # And we should also print whether the likelihood peak is a detection 
    #  according to the limits we just calculated
    limit_at_pos=detec_lim[best_ix[0][0],best_ix[1][0]]
    print('Contrast limit at best fit position: '+str(limit_at_pos))
    if limit_at_pos > best_params[2]:
        print('Detection!')
    else:
        print('No significant detection found')
    
 
    data = {   'chi2'  : chi2,
            'like': like,
            'xys': xys,
            'cons'  : cons,
            'best_params': best_params} 
    return data
 # =========================================================================