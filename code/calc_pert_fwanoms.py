#!/usr/bin/env python3
# coding: utf-8
"""
calc_pert_fwanoms.py

Calculate fresh water anomalies from ECCOv4r4 perturbation experiments

Requires:
- Diagnostic from files from perturbation experiments to ECCOv4r4, using the xx_tauu.0000000129.data.{pert}_{pert_sign} files 
- xx_tauu.0000000129.data.{pert}_{pert_sign}, in ../other_data, were created using make_perts.py

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import glob
import os
import sys
import numpy as np
import xarray as xr
import xmitgcm

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import GRIDDIR, ecco_grid, SOLN_DIR, DATA_DIR, EXPDIR

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "Data from perturbation experiments based on ECCOv4r4 from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

pert_vars = [
    "SALT",
]

YEAR_START = 1992  # Start of simulation
YEAR_PERT = 1996  # Start of perturbation
PERT_NY = YEAR_PERT - YEAR_START
CALC_YEARS = 4  # Number of years to calculate over
FREQ = "mon"
SREF = 35

def load_pert(pert, freq):
    print(f"Loading {freq} {pert} data")
    exf_ds = []
    for evar in pert_vars:
        print(evar)
        try:
            var_ds = xmitgcm.open_mdsdataset(
                data_dir=f"{pert}/diags/{evar}_{freq}_mean",
                grid_dir=GRIDDIR,
                prefix=[f"{evar}_{freq}_mean"],
                geometry="llc",
                delta_t=3600,
                ref_date=STARTDATE,
                read_grid=False,
            )
            var_ds = var_ds.rename({"face": "tile"})
            exf_ds.append(var_ds)
        except:
            print(f"No {freq} {evar}")
    if not exf_ds:
        print(f"No {freq} at all for {pert}")
        return None
    exf_ds = xr.merge(exf_ds)
    return exf_ds


all_pert_plus = glob.glob("../experiments/pert_*pulseplus*")
for pert_plus in all_pert_plus:
    pert_minus = pert_plus.replace("plus", "minus")
    pert_lab = pert_plus.split("/")[-1].replace("plus", "")

    STARTDATE = f"{YEAR_START}-01-01"

    FOUT = f"{DATA_DIR}/perts/{pert_lab}_{FREQ}_pertFW.nc"
    if os.path.isfile(FOUT):
        print(f"Found {FOUT}, skipping")
    else:

        ds_plus = load_pert(pert_plus, FREQ)
        ds_plus["FWC"] = (
            (1 - ds_plus.SALT / SREF) * ecco_grid.drF * ecco_grid.hFacC
        ).sum("k")
        ds_minus = load_pert(pert_minus, FREQ)
        ds_minus["FWC"] = (
            (1 - ds_minus.SALT / SREF) * ecco_grid.drF * ecco_grid.hFacC
        ).sum("k")
        year_min = ds_plus.time.dt.year.min()
        year_max = ds_plus.time.dt.year.max()
        years = [int(year) for year in np.arange(year_min, year_max + 1)]
        if FREQ == "week":
            ds_ctrl = load_pert(f"{EXPDIR}/fwd_26y", FREQ)
        else:
            ds_ctrl = ecco.recursive_load_ecco_var_from_years_nc(
                SOLN_DIR, vars_to_load=pert_vars, years_to_load=years
            )
        ds_ctrl["FWC"] = (
            (1 - ds_ctrl.SALT / SREF) * ecco_grid.drF * ecco_grid.hFacC
        ).sum("k")

        ds_all = []
        ds_pert_lin = ((ds_plus - ds_minus) * 0.5).assign_coords({"pert": "lin"})
        ds_pert_nonlin = ((ds_plus + ds_minus) * 0.5 - ds_ctrl).assign_coords(
            {"pert": "nonlin"}
        )
        ds_all.append(ds_pert_lin)
        ds_all.append(ds_pert_nonlin)

        ds_all = xr.concat(ds_all, "pert", coords="minimal", compat="override")
        if ds_all:
            ds_out = ds_all.isel(
                time=slice(PERT_NY * 12 + 1, PERT_NY * 12 + 1 + 12 * CALC_YEARS)
            )
            print(f"Writing {FOUT}")
            ds_out.attrs.update(attrs)
            ds_out.drop_vars("time_bnds").to_netcdf(FOUT)
