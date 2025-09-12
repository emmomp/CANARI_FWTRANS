#!/usr/bin/env python3
# coding: utf-8
"""
calc_pert_fwanoms.py

Calculate fresh water anomalies from ECCOv4r4 perturbation experiments

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xmitgcm
import glob
import os
import sys
import numpy as np
import xarray as xr
sys.path.insert(0,'/users/emmomp/Python/ECCOv4-py')
import ecco_v4_py as ecco
from inputs import GRIDDIR, ECCO_GRID, SOLN_DIR, DATA_DIR, EXPDIR

pert_vars=['SALT',]
rhoconst = 1029
year_start = 1992 # Start of simulation
year_pert = 1996 # Start of perturbation
pert_ny=year_pert - year_start
calc_years = 4 # Number of years to calculate over
freq='mon'
Sref=35

def load_pert(pert,freq):
    print(f'Loading {freq} {pert} data')
    exf_ds=[]
    for evar in pert_vars:
        print(evar)
        try:
            var_ds = xmitgcm.open_mdsdataset(
                        data_dir=f'{pert}/diags/{evar}_{freq}_mean',
                        grid_dir=GRIDDIR,
                        prefix=[f'{evar}_{freq}_mean'],
                        geometry="llc",
                        delta_t=3600,
                        ref_date=startdate,
                        read_grid=False,
                    )
            var_ds = var_ds.rename({"face": "tile"})
            exf_ds.append(var_ds)
        except:
            print(f'No {freq} {evar}')
    if not exf_ds:
        print(f'No {freq} at all for {pert}')
        return None
    exf_ds=xr.merge(exf_ds)
    return exf_ds

all_pert_plus=glob.glob('../experiments/pert_*pulseplus*')
for pert_plus in all_pert_plus:
    pert_minus=pert_plus.replace('plus','minus')
    pert_lab=pert_plus.split('/')[-1].replace('plus','')
    
    startdate=f'{year_start}-01-01'

    fout=f'{DATA_DIR}/perts/{pert_lab}_{freq}_pertFW.nc'
    if os.path.isfile(fout):
        print(f'Found {fout}, skipping')
    else:
    
        ds_plus=load_pert(pert_plus,freq)
        ds_plus['FWC']=((1-ds_plus.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC).sum('k')
        ds_minus=load_pert(pert_minus,freq)
        ds_minus['FWC']=((1-ds_minus.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC).sum('k')
        if (not ds_plus) and (not ds_minus) :
            continue
        else:
            year_min=ds_plus.time.dt.year.min()
            year_max=ds_plus.time.dt.year.max()
            years=[int(year) for year in np.arange(year_min,year_max+1)]
            if freq=='week':
                ds_ctrl=load_pert(f'{EXPDIR}/fwd_26y',freq)
            else:
                ds_ctrl= ecco.recursive_load_ecco_var_from_years_nc(SOLN_DIR, \
                                           vars_to_load=pert_vars,years_to_load=years)       
            ds_ctrl['FWC']=((1-ds_ctrl.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC).sum('k')
        
        ds_all=[]
        ds_pert_lin=((ds_plus-ds_minus)*0.5).assign_coords({'pert':'lin'})
        ds_pert_nonlin=((ds_plus+ds_minus)*0.5-ds_ctrl).assign_coords({'pert':'nonlin'})
        ds_all.append(ds_pert_lin)
        ds_all.append(ds_pert_nonlin)

        ds_all=xr.concat(ds_all,'pert',coords='minimal',compat='override')
        if ds_all:
            ds_out=ds_all.isel(time=slice(pert_ny*12+1,pert_ny*12+1+12*calc_years))
            print(f'Writing {fout}')
            ds_out.drop_vars('time_bnds').to_netcdf(fout)

