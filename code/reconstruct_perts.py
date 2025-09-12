#!/usr/bin/env python
# coding: utf-8
"""
reconst_denm_highfreq_all.py

Code to reconstruct ECCOv4r4 Denmark Strait Freshwater Flux perturbations from
adjoint sensitivities and applied perturbation fields.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import os
import calendar
import sys
import glob
from datetime import date
import xarray as xr
import numpy as np

sys.path.insert(0, "/users/emmomp/Python")
import xadjoint as xad
import pandas as pd
import convolve_fns
from inputs import (
    adj_diag_map,
    imth,
    ecco_grid,
    EXPDIR,
    GRIDDIR,
    CONV_DIR,
    DATA_DIR,
    eyears,
)


attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 Denmark Strait FW flux reconstructions from Boland et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state estimate, see ecco-group.org",
}

mth_num = [3, 6, 9, 12]
startdates = {"2000": "1996-01-01", "2006": "2002-01-01", "2014": "2010-01-01"}
NT_ADJ = 260
ADJFREQ = 604800
NYEARS_ADJ = 5

perts = glob.glob(f"{DATA_DIR}/*pertfields.nc")
perts = [x.split("/")[-1] for x in perts]
for pert in perts:
    if "pulse" in pert:
        print(pert)
        pert_lab = pert.split("/")[-1].removesuffix("_pertfields.nc")
        ds_pert = xr.open_dataset(f"{DATA_DIR}/{pert}")
        pert_year_max = ds_pert.time.dt.year.max()
        pert_year_min = ds_pert.time.dt.year.min()

        for eyear in eyears:
            print(eyear)
            STARTDATE1 = startdates[eyear]

            for iem in mth_num:
                for lag_mth in [-1, 0, 1]:
                    lmth = np.mod(iem - 1 + lag_mth, 12) + 1
                    print(imth[lmth])
                    STARTDATE = (
                        str(
                            np.datetime64(STARTDATE1, "M")
                            + np.timedelta64(lag_mth, "M")
                        )
                        + "-01"
                    )

                    if iem == 12 and lag_mth == 1:
                        lag0 = np.datetime64(
                            f"{int(eyear)+1}-{lmth:02.0f}-{calendar.monthrange(int(eyear),lmth)[1]}"
                        )
                    else:
                        lag0 = np.datetime64(
                            f"{eyear}-{lmth:02.0f}-{calendar.monthrange(int(eyear),lmth)[1]}"
                        )
                    print(f"start date {STARTDATE}, lag 0 {lag0}")

                    FOUT = f"{CONV_DIR}/{pert_lab}_{eyear}_{imth[lmth]}_recon.nc"
                    if os.path.isfile(FOUT):
                        print(f"Found {FOUT}, skipping")
                        continue
                    expt = f"ad_5y_denstr_horflux_fw_{imth[iem]}_noparam_7d_{eyear}/"
                    print(expt)

                    myexp = xad.Experiment(
                        GRIDDIR,
                        f"{EXPDIR}/{expt}",
                        start_date=STARTDATE,
                        lag0=lag0,
                        nt=NT_ADJ,
                        adj_freq=ADJFREQ,
                    )
                    myexp.load_vars(
                        ["adxx_qnet", "adxx_tauu", "adxx_tauv", "adxx_empmr"]
                    )
                    myexp.data["adxx_tauu"] = myexp.data["adxx_tauu"].rename(
                        {"i_g": "i"}
                    )
                    myexp.data["adxx_tauv"] = myexp.data["adxx_tauv"].rename(
                        {"j_g": "j"}
                    )

                    for var in myexp.data:
                        myexp.data[var].load()

                    ds_out = []
                    for iyear in range(pert_year_min.data, pert_year_max.data + 1):
                        print(iyear)
                        tempdata = myexp.data.copy()
                        if iem == 12 and lag_mth == 1:
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
                        data_convolve = convolve_fns.convolve_ecco(
                            tempdata,
                            ds_pert,
                            adj_diag_map,
                            regrid=False,
                            ecco_grid=ecco_grid,
                            smooth=False,
                            attrs=attrs,
                        )
                        ds_out.append(
                            data_convolve[
                                [x for x in data_convolve if x[-4:] == "_sum"]
                            ]
                        )

                    ds_out = xr.concat(ds_out, "year")
                    ds_out.attrs.update(attrs)
                    ds_out.to_netcdf(FOUT)


print("Done!")
