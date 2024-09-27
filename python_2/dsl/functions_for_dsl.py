#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:03:43 2023

@author: lorenzosangelantoni
"""
import numpy as np
import pandas as pd
from collections import OrderedDict
import xarray as xr
import scipy
from scipy import special
from scipy.stats import norm
from scipy.stats import ks_2samp as ks_2samp
import os
import matplotlib as plt


class structtype():
    pass


def count_dsl(array,thresh,Count,i,cval,pacchetti,dummy_index,index_event_start,index_event_end,index_event,index_events) : 

    while i < len(array)-1 : 
        
        
        if array[i] < thresh : 
            
            i = i + 1;
                
            if i >= len(array) : 
                    
                break
                                    
            c = 1;
                
            while array[i] < thresh : # Start counting (c) time steps
                    
                c = c + 1
                i = i + 1
                    
                if i >= len(array) :
                        
                    break
                        
                    
            if c >= cval : # Start counting (Count) when "c" is above 2 consecutive days
    
                    Count = Count + 1
    
                    pacchetti[Count] = c
                    
                    dummy_index[Count] = i
                    index_event_end[Count]   = int(dummy_index[Count]) #-1 
                    index_event_start[Count] = int(dummy_index[Count] - pacchetti[Count]) 
                    index_event[Count] = np.array(np.arange(int(index_event_start[Count]), int(index_event_end[Count]),1))
                    
                    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
                    dummy_len=len(np.array(np.arange(int(index_event_start[Count]), int(index_event_end[Count]),1)))
                    index_events[Count,range(0,dummy_len)] = np.array(np.arange(int(index_event_start[Count]), int(index_event_end[Count]),1))
                    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
                    
        else :
            i = i + 1

    return pacchetti, index_event, index_events, index_event_start, index_event_end




def adjust_colorbar(data_array,colormap,levels):
    '''
    this function adjust the colorbar regulating upper and lower bins depending if dataset values exceed or not 
    the upper and lower limits of the colorbar. (to be consistent between different CMIP6 models).
    
    Parameters
    ----------
    data_array : xarray dataarray
        dataarray of values (2D lat,lon)
    colormap : string
        string identifying the chosen buit-in matplotloib colormap (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
    levels : numpy 1D array (1D)
        array containing  the bins edges of the colorbar
    
    Returns
    ----------
    my_cmap : matplotlib colormap
        customized colormaps
    
    '''

    #base colormaps 
    base = plt.cm.get_cmap(colormap)
    #max number of bins for the colorbar (N)
    N = len(levels)+1
    #define auxiliar variables useful for the extremes of the colorbar
    step_begin = 0
    step_end = 0
    if data_array.min() > levels[0]: #if minimum dataarray value is higher than lowermost threshold 
        step_begin = 1 #flag to start at bin +1
    if data_array.max() < levels[-1]: #if maximum dataarray value is lower than the uppermost threshold 
        step_end = 1   #flag to stop at bin N-1
    #attributes to create the customized colormap are:
    #number of bins fof the colorbar 
    M = N - step_begin - step_end
    #color range, from 0 to 1 if both extremes are saturated. If not, bins at the at the edges are appropriately removed
    color_list = base(np.linspace(0+step_begin/N, 1-step_end/N, M))
    #colormap name
    cmap_name = base.name + str(M)
    #create customized colorbar
    my_cmap = base.from_list(cmap_name, color_list, M)
    return my_cmap

