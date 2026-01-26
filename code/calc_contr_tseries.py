#!/usr/bin/env python
# coding: utf-8
"""
calc_contr_tseries.py

Code to calculate basin integrals of contributions to ECCOv4r4 Denmark Strait Freshwater Flux.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import os.path
import glob
from datetime import date
import utils as ut
from inputs import CONV_DIR, CONTR_DIR, FCNAME, ecco_convs_2d, eyears

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 Denmark Strait FW flux basin contributions \
       from Boland et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state \
       estimate, see ecco-group.org",
}

my_masks = ut.load_canari_masks()

for eyear in eyears:
    print(f"Calculating contributions for {eyear}")
    cexps_full = glob.glob(f"{CONV_DIR}/transfw_*_7d_{eyear}")
    cexps = [exp.split("/")[-1] for exp in cexps_full]

    for exp in cexps:
        FOUT = f"{CONTR_DIR}/{eyear}/{FCNAME}_contr_tseries_{exp}.nc"
        if os.path.isfile(FOUT):
            print(f"Found {FOUT}, skipping")
        else:
            print(f"Doing {exp}")

        conv_ecco, cexps_mdict, _ = ut.load_ecco_convs(CONV_DIR, eyear, exp=exp)

        dJ_exp = (
            conv_ecco[ecco_convs_2d]
            .sel(exp=exp)
            .chunk({"year": 1, "lag_years": 130, "tile": 13, "j": 90, "i": 90})
        )

        dJpred_ecco_cumsum = dJ_exp.cumsum("lag_years").assign_coords(
            dates=(dJ_exp.dates)
        )
        dJpred_ecco_cumsum["wind_OCE"] = (
            dJpred_ecco_cumsum["adxx_tauuXoceTAUU"]
            + dJpred_ecco_cumsum["adxx_tauvXoceTAUV"]
        )
        dJpred_ecco_cumsum["all_OCE"] = (
            dJpred_ecco_cumsum[ecco_convs_2d].to_array().sum("variable")
        )
        dJpred_ecco_cumsum = dJpred_ecco_cumsum.assign_coords(conv_ecco.coords)
        dJpred_ecco_cumsum = dJpred_ecco_cumsum.assign_coords(
            {"month": cexps_mdict[exp]}
        )

        dJ_tseries = ut.calc_tseries(dJpred_ecco_cumsum, my_masks)
        dJ_tseries.attrs.update(attrs)
        dJ_tseries.to_netcdf(FOUT)
