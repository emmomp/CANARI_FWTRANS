#!/usr/bin/env python
# coding: utf-8
"""
calc_eccov4_meanssh.py

Code to calculate solution mean ECCOv4r4 Sea Surface Height.

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Apr 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import sys
sys.path.insert(0,'/users/emmomp/Python/ECCOv4-py')
import ecco_v4_py as ecco   
from inputs import SOLN_DIR, DATA_DIR

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

ds_ctrl=ecco.recursive_load_ecco_var_from_years_nc(SOLN_DIR,vars_to_load=['SSH'],years_to_load='all')    
ssh_mean=ds_ctrl.SSH.mean('time')
ssh_mean.attrs.update(attrs)
ssh_mean.to_netcdf(f'{DATA_DIR}/eccov4_meanssh.nc')