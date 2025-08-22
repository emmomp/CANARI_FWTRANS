#!/usr/bin/env python
# coding: utf-8
"""
calc_eccov4_varbybasin.py

Code to calculate time series of ECCOv4r4 adjoint sensitivities by basin vs lag.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Aug 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import sys
import calendar
from datetime import date
import xarray as xr
import numpy as np

sys.path.insert(0, "/users/emmomp/Python")
import xadjoint as xad
import utils as ut
from inputs import EXPDIR, GRIDDIR, eyears, mthi, FCNAME

mths = ["Mar", "Jun", "Sep", "Dec"]
ADJ_FREQ = 604800
NT = 260

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 Denmark Strait FW adjoint sensitivities from Boland\
       et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state estimate,\
       see ecco-group.org",
}

my_masks = ut.load_canari_masks()

ds_all = []
for mth in mths:
    ds_all_mth = []
    for year in eyears:
        EXPT = f"ad_5y_denstr_horflux_fw_{mth}_noparam_7d_{year}/"
        STARTDATE = f"{year-4}-01-01"
        lag0 = f"{year}-{mthi[mth]:02.0f}-{calendar.monthrange(int(year),mthi[mth])[1]}"
        print(EXPT, STARTDATE, lag0)
        myexp = xad.Experiment(
            GRIDDIR,
            EXPDIR + EXPT,
            start_date=STARTDATE,
            lag0=lag0,
            nt=NT,
            adj_freq=ADJ_FREQ,
        )
        myexp.load_vars(["adxx_qnet", "adxx_tauu", "adxx_tauv", "adxx_empmr"])

        myexp.data["adxx_tauu"] = -myexp.data["adxx_tauu"].rename({"i_g": "i"})
        myexp.data["adxx_tauv"] = -myexp.data["adxx_tauv"].rename({"j_g": "j"})

        myexp.data = (
            myexp.data.assign_coords({"eyear": year, "month": mth, "fc": myexp.fc})
            .swap_dims({"time": "lag_years"})
            .chunk({"lag_years": 260, "tile": 1, "j": 20, "i": 20})
        )
        ds_tseries = ut.calc_tseries(myexp.data, my_masks).interp(
            lag_years=np.arange(-5, 0, 0.025)
        )
        ds_all_mth.append(ds_tseries)
    ds_all.append(xr.concat(ds_all_mth, "eyear"))
ds_all = xr.concat(ds_all, "month")
ds_all.attrs.update(attrs)
ds_all.to_netcdf(f"../data_out/{FCNAME}_adjsens_basin_tseries_vlag.nc")
