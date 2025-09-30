#!/usr/bin/env python3
# coding: utf-8
"""
generate_masks.py

Generate perturbation xx_ files based on adjoint sensitivities

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xmitgcm
import numpy as np
import xarray as xr
from inputs import ecco_grid, DATA_DIR

# Adjoint sensitivity experiment to base perturbation on
EXPT = "ad_5y_denstr_horflux_fw_Dec_noparam_7d_2000/"

# Perturbation experiment parameters
NY_PERT = 10
PERT_DIR = "../other_data/"

pert_dict = {"NGland": "tauu", "NAlaska": "tauu"}
pert_amp = {"tauu": 0.1}


def write_pert_pulse(pvar, label, pertds, pulse_time="1996-01-01"):
    extra_metadata = xmitgcm.utils.get_extra_metadata(domain="llc", nx=90)
    xxtowrite = pertds.rename({"tile": "face"})
    facets = xmitgcm.utils.rebuild_llc_facets(xxtowrite, extra_metadata)
    compact = xmitgcm.utils.llc_facets_2d_to_compact(facets, extra_metadata)
    time = np.arange(
        np.datetime64("1992-01-01T12:00:00"),
        np.datetime64("2002-01-01T12:00:00"),
        np.timedelta64(7, "D"),
    )
    pert_time = (time < np.datetime64(pulse_time)).sum()

    compact_expand = np.array([])
    for tt in range(0, pert_time):
        compact_expand = np.concatenate([compact_expand, np.zeros_like(compact)])
    for tt in range(0, 4):
        compact_expand = np.concatenate([compact_expand, compact])
    for tt in range(pert_time, len(time)):
        compact_expand = np.concatenate([compact_expand, np.zeros_like(compact)])

    xmitgcm.utils.write_to_binary(
        compact_expand, f"{PERT_DIR}/xx_{pvar}.0000000129.data.{label}_plus"
    )
    xmitgcm.utils.write_to_binary(
        -compact_expand, f"{PERT_DIR}/xx_{pvar}.0000000129.data.{label}_minus"
    )

myexp_data = xr.open_dataset(f"{DATA_DIR}/expts/{EXPT}/tauu.nc")
ds_4y = myexp_data.swap_dims({"time": "lag_years"}).sel(lag_years=-4, method="nearest")

pert_mask = {}

NG_mask = (
    (ecco_grid.XC < -14)
    & (ecco_grid.XC > -100)
    & (ecco_grid.YC > 82)
    & (ecco_grid.YC < 85)
)
pert_mask["NGland"] = NG_mask
NAlaska_mask = (
    (ecco_grid.XC < -125)
    & (ecco_grid.XC > -158)
    & (ecco_grid.YC > 70)
    & (ecco_grid.YC < 74)
)
Bering_mask = (ecco_grid.XC < -154) & (ecco_grid.YC > 73)
mask_comb = NAlaska_mask & ~Bering_mask
pert_mask["NAlaska"] = mask_comb

for pert,var in pert_dict.items():
    pert_ds = ds_4y[f"adxx_{var}"].where(pert_mask[pert])
    if pert == "NAlaska":
        pert_ds = pert_ds.where(pert_ds < 0)
    else:
        pert_ds = pert_ds.where(pert_ds > 0)
    pert_ds = pert_ds * pert_amp[var] / np.abs(pert_ds).max()
    print(f"Writing {pert}")
    write_pert_pulse(var, pert + "JANpulse", pert_ds.fillna(0))
print("Writing finished")
