#!/usr/bin/env python3
# coding: utf-8
"""
denmark_strait_transports.py

Calculate Denmark Strait Transports in the ECCOv4r4 solution

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Feb 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import sys
import xarray as xr
from inputs import SOLN_DIR, ecco_grid

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco

SECTION = "Denmark Strait"
SREF = 35

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

grid = ecco.get_llc_grid(ecco_grid)

ds = ecco.recursive_load_ecco_var_from_years_nc(
    SOLN_DIR,
    vars_to_load=[
        "THETA",
        "SALT",
        "ADVx_SLT",
        "ADVy_SLT",
        "DFxE_SLT",
        "DFyE_SLT",
        "ADVx_TH",
        "ADVy_TH",
        "DFxE_TH",
        "DFyE_TH",
        "UVELMASS",
        "VVELMASS",
        "GM_PsiX",
        "GM_PsiY",
    ],
    years_to_load="all",
)

ds = xr.merge([ecco_grid, ds])

# paramaters for observation fw fluxes, taken from Tesdal et al (2020)
obs_params = {
    "34p8": {
        "sref": 34.8,
        "masktype": "ws",
        "maskW": grid.interp(ds.SALT, "X", boundary="extend") < 34.8,
        "maskS": grid.interp(ds.SALT, "Y", boundary="extend") < 34.8,
    },
    "34p9": {
        "sref": 34.9,
        "masktype": "ws",
        "maskW": grid.interp(ds.SALT, "X", boundary="extend") < 34.9,
        "maskS": grid.interp(ds.SALT, "Y", boundary="extend") < 34.9,
    },
}

print("Data loaded, calculating transports")

fw_transports = []
print(SECTION)
[section_pt1, section_pt2] = ecco.get_section_endpoints(SECTION)
line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(
    section_pt1, section_pt2, ds
)
full_heat_transport = ecco.calc_section_heat_trsp(
    ds, maskS=line_maskS, maskW=line_maskW
)
full_vol_transport = ecco.calc_section_vol_trsp(ds, maskS=line_maskS, maskW=line_maskW)
fw_transports.append(
    ecco.calc_section_fw_trsp(
        ds, Sref=SREF, maskS=line_maskS, maskW=line_maskW, grid=grid
    ).assign_coords({"Sref": SREF})
)
for obs,obs_dict in obs_params.items():
    if obs_dict["masktype"] == "z":
        maskS = line_maskS * obs_dict["mask"]
        maskW = line_maskW * obs_dict["mask"]
    elif obs_dict["masktype"] == "ws":
        maskS = line_maskS * obs_dict["maskS"]
        maskW = line_maskW * obs_dict["maskW"]
    SREF = obs_dict["sref"]
    fw_transports.append(
        ecco.calc_section_fw_trsp(
            ds, Sref=SREF, maskS=maskS, maskW=maskW, grid=grid
        ).assign_coords({"Sref": SREF})
    )
full_fw_transport = xr.concat(fw_transports, "Sref")
full_transport = xr.merge(
    [full_heat_transport, full_vol_transport, full_fw_transport], compat="override"
)
full_transport = full_transport.assign_coords({"section": SECTION})

print("All calculations done, writing to file")
full_transport.attrs.update(attrs)
full_transport.to_netcdf(
    f'../data_out/{SECTION.replace(" ","_")}_full_depth_transports.nc'
)

print("All done")
