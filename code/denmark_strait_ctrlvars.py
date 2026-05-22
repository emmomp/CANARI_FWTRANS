#!/usr/bin/env python3
# coding: utf-8
"""
denmark_strait_ctrlvars.py

Extract variables from Denmark Strait in the ECCOv4r4 solution

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jun 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import sys
import numpy as np
import xarray as xr

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import EXPDIR, ecco_grid, GRIDDIR

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

SECTION = "Denmark Strait"
startdate=np.datetime64('1992-01-01T12:00:00')
Sref=35

[section_pt1, section_pt2] = ecco.get_section_endpoints(SECTION)
line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(    
    section_pt1, section_pt2, ecco_grid
)

extra_variables = dict( SSH = dict(dims=['j','i'], attrs=dict(standard_name='sea_surface_height_above_geoid', long_name='Dynamic sea surface height anomaly', units='m')))

ds_ctrl=[]
for var in ['MXLDEPTH','SIarea','SALT','THETA','SSH']:
    ds_ctrl.append(ecco.load_ecco_vars_from_mds(f'{EXPDIR}/fwd_26y/diags/{var}_mon_mean',mds_grid_dir=GRIDDIR,mds_files=f'{var}_mon_mean',output_freq_code='AVG_MON',
                                                    model_start_datetime=startdate,read_grid=False,extra_variables=extra_variables
                                                   )) 
ds_ctrl=xr.merge(ds_ctrl+[ecco_grid,])

ds_ctrl['FWC']=((1-ds_ctrl.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC*ecco_grid.rA).sum('k')
ds_ctrl['FWC_z']=((1-ds_ctrl.SALT/Sref)*ecco_grid.drF*ecco_grid.hFacC*ecco_grid.rA)
ctrl_section=ds_ctrl.chunk({'time':121,'k':50,'tile':1,'j':20,'i':20}).where(line_maskC,drop=True).squeeze().stack({'ji':['j','i']}).dropna('ji').sortby('XC')
ctrl_section.attrs.update(attrs)
ctrl_section.time.encoding.clear()
ctrl_section.reset_index('ji').to_netcdf(f'../data_out/{SECTION.replace(" ","_")}_vars_mon.nc')

# Define halocline and find base
base_mask=(ctrl_section['SALT']<=34.8)&((-ctrl_section['Z'])<ctrl_section['Depth'])
base_depth=ecco_grid.Z[base_mask.argmin('k').load()-1]
base_mask=xr.where(ecco_grid.Z>base_depth,1,0)
base_mask.name='halocline_mask'
base_mask=base_mask.assign_coords({'Z':ecco_grid['Z']})
base_depth.name='halocline_depth'
ds_out=xr.merge([base_mask,base_depth.drop_vars(['Z','PHrefC','drF'])])
ds_out.attrs.update(attrs)
ds_out.reset_index('ji').to_netcdf(f'../data_out/{SECTION.replace(" ","_")}_halocline_mon.nc')
