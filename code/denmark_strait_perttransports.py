#!/usr/bin/env python3
# coding: utf-8
"""
denmark_strait_perttransports.py

Calculate Denmark Strait FW & Vol 2-D Transports in perturbation experiments to ECCOv4r4

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import sys
import xarray as xr
from inputs import SOLN_DIR, ecco_grid, DATA_DIR

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco

SECTION = "Denmark Strait"
SECTION_LABEL='Den'
SREF = 35

[section_pt1, section_pt2] = ecco.get_section_endpoints(SECTION)
line_maskC, line_maskW, line_maskS = ecco.get_section_line_masks(
    section_pt1, section_pt2, ecco_grid
)

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

grid = ecco.get_llc_grid(ecco_grid)

exps = [
    "pert_10y_tauu_NGlandJANpulseplus",
    "pert_10y_tauu_NGlandJANpulseminus",
    "pert_10y_tauu_NAlaskaJANpulseplus",
    "pert_10y_tauu_NAlaskaJANpulseminus",
]

for exp in exps:
    print(f"Calculating {exp} {SECTION} transports")
    ds = ecco.recursive_load_ecco_var_from_years_nc(
        SOLN_DIR,
        vars_to_load=[
            "SALT",
            "UVELMASS",
            "VVELMASS",
            "GM_PsiX",
            "GM_PsiY",
        ],
        years_to_load="all",
    )

    ds = xr.merge([ecco_grid, ds])

    print("Data loaded, calculating transports")

    full_vol_transport = ecco.calc_section_vol_trsp(
        ds, maskC=line_maskC, grid=grid, along_section=True
    )
    full_fw_transport = ecco.calc_section_fw_trsp(
        ds, Sref=SREF, maskC=line_maskC, along_section=True, grid=grid
    ).assign_coords({"Sref": SREF})

    full_transport = xr.merge(
        [full_vol_transport, full_fw_transport], compat="override"
    )
    full_transport = full_transport.assign_coords({"exp": exp})

    print("All calculations done, writing to file")
    full_transport.attrs.update(attrs)
    full_transport.reset_index("ij").to_netcdf(
        f"{DATA_DIR}/perts/{exp}_{SECTION_LABEL}_2dtransports.nc"
    )

print("All done")
