#!/usr/bin/env python3
# coding: utf-8
"""
adj_stats.py

Calculate mean and abs mean time series of adjoint sensitivities

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import numpy as np
import ecco_v4_py as ecco
sys.path.insert(0,'/users/emmomp/Python')
import xadjoint as xad
import utils as ut
from inputs import *

mths=['Mar','Jun','Sep','Dec']
adj_freq=604800
nt=260
adj_vars=['adxx_qnet','adxx_empmr','adxx_tauu','adxx_tauv']

for mth in mths:
    for year in eyears:
        expt=f'ad_5y_denstr_horflux_fw_{mth}_noparam_7d_{year}/'
        startdate=f'{int(year)-4}-01-01'
        lag0=f'{year}-{mthi[mth]:02.0f}-{calendar.monthrange(int(year),mthi[mth])[1]}'
        print(expt,startdate,lag0)
        myexp = xad.Experiment(GRIDDIR,f'{EXPDIR}/{expt}',start_date=startdate,lag0=lag0,nt=nt,adj_freq=adj_freq)
        myexp.load_vars(adj_vars)
    
        myexp.data['adxx_tauu']=-myexp.data['adxx_tauu'].rename({'i_g':'i'})
        myexp.data['adxx_tauv']=-myexp.data['adxx_tauv'].rename({'j_g':'j'})
    
        myexp.data=myexp.data.assign_coords({'eyear':year,'month':mth,'fc':myexp.fc}).swap_dims({'time':'lag_years'})
        data_stats=ut.calc_tseries(myexp.data)
        data_stats.to_netcdf({f'{EXPDIR}/{expt}/{expt}_stats.nc'})