#!/usr/bin/env python
# coding: utf-8
"""
calc_eccov4_varbybasin.py

Code to calculate variance of ECCOv4r4 forcing variables by basin.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Apr 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr 
import utils as ut
from inputs import EXPDIR, DATA_DIR, ecco_grid

ds_climanom=xr.open_dataset(f'{EXPDIR}/fwd_26y/exf_climanoms.nc').drop_vars(['oceTAUX','oceTAUY'])
ds_ca_monvar=ds_climanom.groupby(ds_climanom.time.dt.month).var('time')

my_masks=ut.load_canari_masks()
ds_ca_mv_bybasin=ut.calc_tseries(ds_ca_monvar,my_masks,weight=ds_ca_monvar.rA)
rA_bybasin=[]
for basin in my_masks:
    rA_bybasin.append(ecco_grid['rA'].where(my_masks[basin]).sum().assign_coords({"mask": basin}))
rA_bybasin=xr.concat(rA_bybasin,'mask',coords='minimal')
ds_basinmean=ds_ca_mv_bybasin/rA_bybasin

ds_basinmean.to_netcdf(f'{DATA_DIR}/eccov4_varbybasin.nc')