#!/usr/bin/env python3
# coding: utf-8
"""
denmark_strait_ctrlvars.py

Extract variables from Denmark Strait in perturbations to the ECCOv4r4 solution

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jun 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import sys
import xarray as xr
import numpy as np

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import SOLN_DIR, ecco_grid, EXPDIR, GRIDDIR

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "Data from perturbation experiments based on ECCOv4r4 from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

SECTION = "Denmark Strait"
pert='pert_10y_tauu_NGlandJANpulse'
startdate=np.datetime64('1992-01-01T12:00:00')
Sref=35

[section_pt1, section_pt2] = ecco.get_section_endpoints(SECTION)
line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(    
    section_pt1, section_pt2, ecco_grid
)

extra_variables = dict( SSH = dict(dims=['j','i'], attrs=dict(standard_name='sea_surface_height_above_geoid', long_name='Dynamic sea surface height anomaly', units='m')))

ds_plus=[]
ds_minus=[]
for var in ['MXLDEPTH','SIarea','SALT','THETA','SSH']:
    ds_plus.append(ecco.load_ecco_vars_from_mds(f'{EXPDIR}/{pert}plus/diags/{var}_mon_mean',mds_grid_dir=GRIDDIR,mds_files=f'{var}_mon_mean',output_freq_code='AVG_MON',
                                                    model_start_datetime=startdate,read_grid=False,extra_variables=extra_variables
                                                   ))
    ds_minus.append(ecco.load_ecco_vars_from_mds(f'{EXPDIR}/{pert}minus/diags/{var}_mon_mean',mds_grid_dir=GRIDDIR,mds_files=f'{var}_mon_mean',output_freq_code='AVG_MON',
                                                    model_start_datetime=startdate,read_grid=False,extra_variables=extra_variables
                                                    ))
ds_plus=xr.merge(ds_plus).assign_coords({'type':'plus'})
ds_minus=xr.merge(ds_minus).assign_coords({'type':'minus'})
ds_perts=xr.concat([ds_plus,ds_minus],'type')
ds_perts['FWC']=((1-ds_perts.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC*ecco_grid.rA).sum('k')
ds_perts['FWC_z']=((1-ds_perts.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC*ecco_grid.rA)
ds_section=ds_perts.chunk({'time':121,'k':50,'tile':1,'j':20,'i':20}).where(line_maskC,drop=True).squeeze().stack({'ji':['j','i']}).dropna('ji').sortby('XC')
ds_section.attrs.update(attrs)
ds_section.reset_index('ji').to_netcdf(f'../data_out/perts/{SECTION.replace(" ","_")}_{pert}.nc')