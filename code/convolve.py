#!/usr/bin/env python
# coding: utf-8
"""
convolve.py

Code to reconstruct ECCOv4r4 Denmark Strait Freshwater Flux from adjoint sensitivities.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Feb 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import os
import calendar
import sys
import xarray as xr
import numpy as np

sys.path.insert(0, "/users/emmomp/Python")
import xadjoint as xad
import pandas as pd
from datetime import date
import convolve_fns

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'ECCOv4r4 Denmark Strait FW flux reconstructions from Boland et al. 2025 (inprep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the ECCOv4r4 global ocean state estimate, see ecco-group.org'}

DATA_DIR = "/data/expose/ECCOv4-r4/Version4/Release4/input_forcing/"
ROOTDIR = "/users/emmomp/data/canari/experiments/"
GRIDDIR = "/users/emmomp/data/orchestra/grid2/"
CONV_DIR = "../data_out/denm_X_ECCOclimanom_hfreq"
TRANSP = "fw"

ecco_grid = xr.open_dataset("~/data/orchestra/other_data/ECCO_r3_alt/ECCOv4r3_grid.nc")
eyear = sys.argv[1]  # Should be 2000, 2006 or 2014
print(eyear)
if eyear not in ["2000", "2006", "2014"]:
    raise ValueError("Expecting 2000, 2006 or 2014")

startdates = {"2000": "1996-01-01", "2006": "2002-01-01", "2014": "2010-01-01"}
startdate1 = startdates[eyear]
NT_ADJ = 260
ADJFREQ = 604800

all_mths = dict(
    zip(
        range(0, 12),
        [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ],
    )
)

mth_num = [3 - 1, 6 - 1, 9 - 1, 12 - 1]
# adj_diag_map={'adxx_qnet':['EXFqnet','oceQnet'],'adxx_tauu':['EXFtauu','oceTAUU'],'adxx_tauv':['EXFtauv','oceTAUV'],'adxx_empmr':['EXFempmr','oceFWflx']}
adj_diag_map = {
    "adxx_qnet": [
        "oceQnet",
    ],
    "adxx_empmr": [
        "oceFWflx",
    ],
}
# adj_diag_map={'adxx_tauu':['EXFtauu','oceTAUU'],'adxx_tauv':['EXFtauv','oceTAUV']}
# adj_diag_map={'adxx_qnet':['EXFqnet',],'adxx_empmr':['EXFempmr']}
# adj_diag_map={'adxx_tauu':['EXFtauu',],'adxx_tauv':['EXFtauv',]}

DATA_DIR = "/users/emmomp/data/canari/experiments/fwd_26y/"

ds_climanom = xr.open_dataset(f"{DATA_DIR}/exf_climanoms.nc")

for iem in mth_num:
    # for iem in [10,]:
    expt = f"ad_5y_denstr_horflux_fw_{all_mths[iem]}_noparam_7d_{eyear}/"
    print(expt)
    for lag_mth in [-1, 0, 1]:
        #   for lag_mth in [-1,0]:
        print(all_mths[np.mod(iem + lag_mth, 12)])
        STARTDATE = (
            str(np.datetime64(startdate1, "M") + np.timedelta64(lag_mth, "M")) + "-01"
        )
        lmth = np.mod(iem + lag_mth, 12) + 1
        if iem == 11 and lag_mth == 1:
            lag0 = np.datetime64(
                f"{int(eyear)+1}-{lmth:02.0f}-{calendar.monthrange(int(eyear),lmth)[1]}"
            )
        else:
            lag0 = np.datetime64(
                f"{eyear}-{lmth:02.0f}-{calendar.monthrange(int(eyear),lmth)[1]}"
            )
        print(f"start date {STARTDATE}, lag 0 {lag0}")

        myexp = xad.Experiment(
            GRIDDIR,
            ROOTDIR + expt,
            start_date=STARTDATE,
            lag0=lag0,
            nt=NT_ADJ,
            adj_freq=ADJFREQ,
        )
        myexp.load_vars(["adxx_qnet", "adxx_tauu", "adxx_tauv", "adxx_empmr"])
        myexp.data["adxx_tauu"] = -myexp.data["adxx_tauu"].rename({"i_g": "i"})
        myexp.data["adxx_tauv"] = -myexp.data["adxx_tauv"].rename({"j_g": "j"})

        myexp.data = myexp.data.assign_coords({"TRANSP": TRANSP, "year": eyear})

        if lag_mth == 0:
            of = f"trans{TRANSP}_{all_mths[iem]}_noparam_7d_{eyear}"
        else:
            of = f"trans{TRANSP}_{all_mths[np.mod(iem + lag_mth, 12)]}_from{all_mths[iem]}_7d_{eyear}"

        if not os.path.isdir(f"{CONV_DIR}/{of}"):
            os.mkdir(f"{CONV_DIR}/{of}")
            for iyear in range(1992, 2018):
                os.mkdir(f"{CONV_DIR}/{of}/{iyear}")

        for iyear in range(1992, 2018):
            tempdata = myexp.data.copy()
            if iem == 11 and lag_mth == 1:
                tempdata["time"] = (
                    pd.DatetimeIndex(tempdata["time"])
                    + pd.DateOffset(years=iyear - int(eyear) - 1)
                ).to_numpy()
            else:
                tempdata["time"] = (
                    pd.DatetimeIndex(tempdata["time"])
                    + pd.DateOffset(years=iyear - int(eyear))
                ).to_numpy()
            tempdata = tempdata.assign_coords({"year": iyear})
            data_convolve=convolve_fns.convolve_ecco(
                tempdata,
                ds_climanom,
                adj_diag_map,
                regrid=False,
                ecco_grid=ecco_grid,
                dir_out=f"{CONV_DIR}/{of}/{iyear}",
                smooth=True,
                overwrite=False,
                attrs=attrs
            )

print("Done!")
