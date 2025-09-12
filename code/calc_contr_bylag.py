#!/usr/bin/env python3
# coding: utf-8
"""
calc_contr_bylag.py

Calculate 2-D contributions to Denmark Strait FW transport in ECCOv4r4 from convolutions

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import utils as ut
from inputs import CONV_DIR, CONTR_DIR, ecco_convs_2d, eyears

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 Denmark Strait FW flux reconstructions from Boland\
    et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state estimate,\
    see ecco-group.org",
}

lags = [0, -0.25, -0.5, -1.5, -4]
lag_labels = [
    "0 to -3m lag",
    "-3m to -6m lag",
    "-6m to -18m lag",
    "-18m to -4 y lag",
    "0 to -4 y lag",
]
lag_flabels = ["0to3m", "3to6m", "6to18m", "18to4y", "0to4y"]

for eyear in eyears:
    print(f"Calculating contributions for {eyear}")

    conv_ecco, cexps_mdict, cexps_edict = ut.load_ecco_convs(CONV_DIR, eyear)
    conv_ecco = conv_ecco[ecco_convs_2d].chunk("auto")

    for ilag in range(0, len(lags) - 1):
        print(f"Calculating {lag_labels[ilag]}")
        FOUT = f"{CONTR_DIR}/{eyear}/ecco_2dconvs_{eyear}_{lag_flabels[ilag]}.nc"
        ds_exp = (
            conv_ecco[ecco_convs_2d]
            .sel(lag_years=slice(lags[ilag], lags[ilag + 1]))
            .sum("lag_years")
            .squeeze()
        )
        ds_exp = ds_exp.assign_coords(
            {
                "lag_range": lag_labels[ilag],
                "month": ("exp", [cexps_mdict[exp] for exp in ds_exp.exp.data]),
            }
        )
        ds_exp.to_netcdf(FOUT)
    ilag += 1
    FOUT = f"{CONTR_DIR}/{eyear}/ecco_convs_{eyear}_{lag_flabels[ilag]}.nc"
    ds_exp = (
        conv_ecco[ecco_convs_2d]
        .sel(lag_years=slice(0, lags[-1]))
        .sum("lag_years")
        .squeeze()
    )
    ds_exp = ds_exp.assign_coords(
        {
            "lag_range": lag_labels[-1],
            "month": [cexps_mdict[exp] for exp in ds_exp.exp.data],
        }
    )
    ds_exp.attrs.update(attrs)
    ds_exp.to_netcdf(FOUT)
