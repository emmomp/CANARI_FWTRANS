#!/usr/bin/env python3
# coding: utf-8
"""
denmark_strait_ctrlvars.py

Extract variables from Denmark Strait in the ECCOv4r4 solution

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jun 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import sys
import xarray as xr

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import SOLN_DIR, ecco_grid

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

SECTION = "Denmark Strait"

[section_pt1, section_pt2] = ecco.get_section_endpoints(SECTION)
line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(    
    section_pt1, section_pt2, ecco_grid
)

ds_ctrl=[]
for var in ['MXLDEPTH','SIarea','SALT','THETA','SSH']:
    ds_ctrl.append(ecco.recursive_load_ecco_var_from_years_nc(SOLN_DIR, \
                                               vars_to_load=[var,    ],years_to_load='all')    )  
ds_ctrl=xr.merge(ds_ctrl+[ecco_grid,])

ctrl_section=ds_ctrl.chunk({'time':121,'k':50,'tile':1,'j':20,'i':20}).where(line_maskC,drop=True).squeeze().stack({'ji':['j','i']}).dropna('ji').sortby('XC')
ctrl_section.time.encoding.clear()
ctrl_section.reset_index('ji').to_netcdf(f'../data_out/{SECTION.replace(" ","_")}_vars.nc')

# Define halocline and find base
base_mask=(ctrl_section['SALT']<=34.8)&((-ctrl_section['Z'])<ctrl_section['Depth'])
base_depth=ecco_grid.Z[base_mask.argmin('k').load()-1]
base_mask=xr.where(ecco_grid.Z>depth,1,0)
base_mask.name='halocline_mask'
base_mask=base_mask.assign_coords({'Z':ecco_grid['Z']})
base_depth.name='halocline_depth'
ds_out=xr.merge([base_mask,base_depth.drop_vars(['Z','PHrefC','drF'])])
ds_out.attrs.update(attrs)
ds_out.reset_index('ji').to_netcdf(f'../data_out/{SECTION.replace(" ","_")}_halocline.nc')
