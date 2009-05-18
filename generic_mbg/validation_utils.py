from __future__ import division
import pymc as pm
import pylab as pl
import tables as tb
import numpy as np
from scipy import integrate
from histogram_utils import *
import time
import os

__all__ = ['validation_fns']

def roc(p_samps, n_samps, pos, neg):
    """
    Plots the receiver-operator characteristic.
    """
    t = np.linspace(0,1,500)
    tp = [1]
    fp = [1]
    tot_pos = np.float(np.sum(pos))
    tot_neg = np.float(np.sum(neg))
    
    marginal_p = np.mean(p_samps,axis=0)
    
    for i in xrange(len(t)):
        where_yes = np.where(marginal_p>t[i])
        where_no = np.where(marginal_p <= t[i])
        fp_here = np.sum(neg[where_yes])/tot_neg
        if fp_here < 1:
            tp.append(np.sum(pos[where_yes])/tot_pos)
            fp.append(fp_here)
        if fp_here == 0:
            break
            
    fp = np.array(fp)
    tp = np.array(tp)
        
    pl.fill(np.concatenate((fp,[1,1])),np.concatenate((tp,[0,1])), facecolor='.8',edgecolor='k',linewidth=1)
    pl.plot([0,1],[0,1],'k-.')
    
    pl.xlabel('False positive rate')
    pl.ylabel('True positive rate')
    
    AUC = -np.sum(np.diff(fp)*(tp[1:] + tp[:-1]))/2.
    if np.isnan(AUC):
        raise ValueError, 'AUC is NaN.'
    pl.title('AUC: %.3f'%AUC)
    pl.axis('image')
    
def scatter(p_samps, n_samps, pos, neg):
    """
    Plots the expected fraction positive against the observed fraction 
    positive.
    """
    p_pred = np.mean(p_samps, axis=0)
    p_obs = pos/(pos+neg).astype('float')
    
    pl.plot(p_obs, p_pred,'k.', markersize=2)
    urc = max(p_obs.max(),p_pred.max())*1.1
    
    pl.plot([0,urc],[0,urc],'k-.')
    pl.xlabel('Observed fraction positive')
    pl.ylabel('Expected fraction positive')
    pl.axis('image')
    
def coverage(p_samps, n_samps, pos, neg):
    """
    Plots the coverage plot in the lower right panel of mbgw.
    """
    obs_quantiles = np.array([np.sum(n_samps[:,i] > pos[i]) for i in xrange(len(pos))], dtype='float')/n_samps.shape[0]
    pt = np.linspace(0,1,500)
    cover = np.array([np.sum(obs_quantiles<pti) for pti in pt], dtype='float')/len(obs_quantiles)
    pl.plot([0,1],[0,1],'k-.')
    pl.plot(pt,cover,'k-')
    pl.xlabel('Predictive quantile')
    pl.ylabel('Fraction of observations below quantile')
    pl.axis('image')
            
validation_fns = [roc,scatter,coverage]