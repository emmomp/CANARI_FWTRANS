#!/usr/bin/env python3
# coding: utf-8
"""
montecarlo_function.py

Calculate bootstrapped p-values for correlation and explained variance calculations

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jan 2026

@author: emmomp@bas.ac.uk Emma J D Boland, based on original matlab code by Chris Hughes, Liverpool
"""
import numpy as np

def MonteCarlo_function(series1, series2, R_actual, times_replicate=1000):
    series1M=montecarlo(series1,times_replicate)
    series2M=montecarlo(series2,times_replicate)
    RR=np.zeros((times_replicate,times_replicate))
    for i in range (0,times_replicate):
        for j in range(0,times_replicate):
            corr=np.corrcoef(series1M[:,i], series2M[:,j])
            RR[i,j]=corr[0,1]
    RRR=RR.ravel()
    RRR_above=RRR[np.abs(RRR)>=np.abs(R_actual)]

    fracLevP=len(RRR_above)/len(RRR)

    return fracLevP

def MonteCarlo_EV(series1, series2, EV_actual, times_replicate=1000):
    series1M=montecarlo(series1,times_replicate)
    series2M=montecarlo(series2,times_replicate)
    RR=np.zeros((times_replicate,times_replicate))
    for i in range (0,times_replicate):
        var_sol=np.var(series1M[:,i])
        for j in range(0,times_replicate):
            RR[i,j]=1-np.var(series1M[:,i]-series2M[:,j])/var_sol
    RRR=RR.ravel()
    RRR_above=RRR[RRR>=EV_actual]

    fracLevP=len(RRR_above)/len(RRR)

    return fracLevP

def montecarlo( tin,numt ):
    """
    generate numt surrogate timeseries
    with same characteristics as tin

    input: 
    tin = time series (assumes a 1D time series, converts internally to 1D if not)
    numt = number of surrogate time series to return

    output: 
    tsurr = array dimensions (tlen,numt), where tlen is the
            number of elements in tin. Each column of tsurr
            is a single surrogate timeseries with the same
            power as tin at each frequency, but random phase.
    """

    # double length of timeseries by adding its reflection,
    # and take fourier transform
    tft=np.fft.fft( np.concatenate([tin.ravel(),np.flip(tin.ravel())]))
    tlen=tin.size
    # make numt sets of complex numbers with amplitude 1 and random phase

    trand=np.exp(1j*2.0*np.pi*np.random.rand(tlen-1,numt))
    # make an initial value for the numt transforms, the same as for the
    # transform of tin
    tr1=trand[0,:]*0.0+1.0
    #
    # make a middle value for the numt transforms which is purely real
    #
    trmid=np.sin(np.random.rand(numt))
    #
    # make numt transformed time series, each with same mean as original
    # 2nd half is reversed complex conjugate of 1st half to ensure 
    # it represents the fft of a real timeseries.
    #
    trrev=np.flip(np.conjugate(trand),axis=0)
    trand=np.concatenate([tr1[np.newaxis],trand,trmid[np.newaxis],trrev])
    tft_new=np.multiply(np.tile(tft,(numt,1)).T,trand)
    #
    # transform back to timeseries, and extract first half only
    #
    tsurr=np.fft.ifft(tft_new,axis=0)

    return np.real(tsurr[:tlen])

