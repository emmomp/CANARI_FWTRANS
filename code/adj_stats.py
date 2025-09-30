#!/usr/bin/env python3
# coding: utf-8
"""
adj_stats.py

Calculate mean and abs mean time series of adjoint sensitivities, write adjoint sensitivities to file

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import calendar
import sys
from datetime import date

sys.path.insert(0,'/users/emmomp/Python/ECCOv4-py')
import ecco_v4_py as ecco
sys.path.insert(0, "/users/emmomp/Python")
import xadjoint as xad
import utils as ut
from inputs import GRIDDIR, EXPDIR, eyears, mthi, DATA_DIR

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'ECCOv4r4 Denmark Strait FW flux reconstructions from Boland et al. 2025 (inprep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the ECCOv4r4 global ocean state estimate, see ecco-group.org'}

mths = ["Mar", "Jun", "Sep", "Dec"]
ADJ_FREQ = 604800
NT = 260
adj_vars = ["adxx_qnet", "adxx_empmr", "adxx_tauu", "adxx_tauv"]

for mth in mths:
    for year in eyears:
        EXPT = f"ad_5y_denstr_horflux_fw_{mth}_noparam_7d_{year}/"
        STARTDATE = f"{int(year)-4}-01-01"
        lag0 = f"{year}-{mthi[mth]:02.0f}-{calendar.monthrange(int(year),mthi[mth])[1]}"
        print(EXPT, STARTDATE, lag0)
        myexp = xad.Experiment(
            GRIDDIR,
            f"{EXPDIR}/{EXPT}",
            start_date=STARTDATE,
            lag0=lag0,
            nt=NT,
            adj_freq=ADJ_FREQ,
        )
        myexp.load_vars(adj_vars)

        myexp.data["adxx_tauu"] = -myexp.data["adxx_tauu"].rename({"i_g": "i"})
        myexp.data["adxx_tauv"] = -myexp.data["adxx_tauv"].rename({"j_g": "j"})

        myexp.data = myexp.data.assign_coords(
            {"eyear": year, "month": mth}).swap_dims({"time": "lag_years"})
        if hasattr(myexp,'fc'):
            myexp.data = myexp.data.assign_coords({ "fc": myexp.fc})
        myexp.to_nctiles(split_timesteps=False)
        
        data_stats = ut.calc_tseries(myexp.data.chunk({'tile':-1,'j':-1,'i':-1}))
        data_stats.attrs.update(attrs)
        data_stats.to_netcdf(f"{EXPDIR}/{EXPT}/{EXPT[:-1]}_stats.nc")
