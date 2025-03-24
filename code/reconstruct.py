#!/usr/bin/env python
# coding: utf-8
"""reconstruct.py 

Code to reconstruct ECCOv4r4 Denmark Strait Freshwater Flux from adjoint sensitivities.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""

import glob
from datetime import date
import dask
import xarray as xr
import numpy as np
import inputs
import utils as ut

EXPT = "fwd_26y"

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 Denmark Strait FW flux reconstructions from Boland\
    et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state estimate,\
    see ecco-group.org",
}

print("Loading solution")
fc = ut.get_soln(FCNAME, f"{EXPDIR}/{EXPT}")
fc_climanom, fc_mth = ut.soln_anoms(fc)

for eyear in eyears[1:]:
    print(f"Calculating reconstructions for {eyear}")

    conv_ecco,cexps_mdict,cexps_edict=ut.load_ecco_convs(CONV_DIR, eyear)
    conv_ecco=conv_ecco[ecco_convs['all']].chunk('auto')

    dJpred_ecco = conv_ecco.sum(dim="lag_years")
    dJpred_ecco["wind_EXF"] = (
        dJpred_ecco["adxx_tauuXEXFtauu_sum"] + dJpred_ecco["adxx_tauvXEXFtauv_sum"]
    )
    dJpred_ecco["wind_OCE"] = (
        dJpred_ecco["adxx_tauuXoceTAUU_sum"] + dJpred_ecco["adxx_tauvXoceTAUV_sum"]
    )
    dJpred_ecco["all_OCE"] = dJpred_ecco[ecco_convs["OCE"]].to_array().sum("variable")
    dJpred_ecco["all_EXF"] = dJpred_ecco[ecco_convs["EXF"]].to_array().sum("variable")

    dJpred_ecco_cumsum = conv_ecco.cumsum("lag_years").assign_coords(
        dates=dJpred_ecco.dates
    )
    dJpred_ecco_cumsum["wind_EXF"] = (
        dJpred_ecco_cumsum["adxx_tauuXEXFtauu_sum"]
        + dJpred_ecco_cumsum["adxx_tauvXEXFtauv_sum"]
    )
    dJpred_ecco_cumsum["wind_OCE"] = (
        dJpred_ecco_cumsum["adxx_tauuXoceTAUU_sum"]
        + dJpred_ecco_cumsum["adxx_tauvXoceTAUV_sum"]
    )
    dJpred_ecco_cumsum["all_OCE"] = (
        dJpred_ecco_cumsum[ecco_convs["OCE"]].to_array().sum("variable")
    )
    dJpred_ecco_cumsum["all_EXF"] = (
        dJpred_ecco_cumsum[ecco_convs["EXF"]].to_array().sum("variable")
    )
    dJpred_ecco_cumsum = dJpred_ecco_cumsum.assign_coords(conv_ecco.coords)

    dJ_vars = dJpred_ecco.squeeze().stack(yearexp=["exp", "year"]).sortby("dates")

    cum_ev = xr.open_dataset(
        f"{EV_DIR}/{eyear}/horflux_fw_denm_cumev_bylag_byvar_{eyear}.nc"
    )
    cum_ev_bym = xr.open_dataset(
        f"{EV_DIR}/{eyear}/horflux_fw_denm_cumev_bylag_byvar_bymonth_{eyear}.nc"
    )
    cum_ev_bym["lag_years"] = cum_ev["lag_years"]
    lagmax = cum_ev_bym.idxmax("lag_years").squeeze().load()
    lagmax_ds = xr.concat([lagmax[var] for var in lagmax], "var").assign_coords(
        var=list(lagmax)
    )
    lagmax_ds.name = "lag_max"

    for var in ecco_convs["all"] + ["wind_OCE", "wind_EXF", "all_OCE", "all_EXF"]:
        print("Writing reconstructions")
        print(f"Full reconstruction, {var}")
        YY = (
            dJ_vars[var]
            .swap_dims({"yearexp": "dates"})
            .rename({"dates": "time"})
            .assign_coords({"eyear": eyear})
        )
        YY.attrs.update(attrs)
        YY.drop_vars(["yearexp", "exp"]).to_netcdf(
            f"{RECON_DIR}/denstr_fwflux_4yrecon_{eyear}_{var}.nc"
        )

        # Peak year reconstruction
        lmax = cum_ev[var].idxmax("lag_years")
        dJpred_maxev = (
            dJpred_ecco_cumsum[var]
            .sel(lag_years=lmax, method="nearest")
            .stack(yearexp=["exp", "year"])
            .sortby("dates")
        )
        dJpred_maxev = (
            dJpred_maxev.swap_dims({"yearexp": "dates"})
            .drop_vars("time")
            .rename({"dates": "time"})
        )
        print(f"{var} recon max {lmax:2.1f}y")

        YY = dJpred_maxev.assign_coords({"eyear": eyear})
        print(f"Peak reconstruction, {var}")
        YY.attrs.update(attrs)
        YY.drop_vars(["yearexp", "exp"]).to_netcdf(
            f"{RECON_DIR}/denstr_fwflux_peakEVrecon_{eyear}_{var}.nc"
        )

        # Peak seasonal reconstruction
        lmax = lagmax_ds.sortby("month").sel(var=var)
        dJpred_maxev = (
            dJpred_ecco_cumsum[var]
            .squeeze()
            .swap_dims({"exp": "month"})
            .sortby("month")
            .sel(lag_years=lmax, method="nearest")
            .sortby("month")
        )
        dJpred_maxev = dJpred_maxev.stack(yearmonth=["year", "month"]).drop_vars("time")

        YY = (
            dJpred_maxev.swap_dims({"yearmonth": "dates"})
            .rename({"dates": "time"})
            .assign_coords({"eyear": eyear})
        )
        YY.attrs.update(attrs)
        print(f"Peak monthly reconstruction, {var}")
        YY.drop_vars(["yearmonth", "exp"]).to_netcdf(
            f"{RECON_DIR}/denstr_fwflux_mthEVrecon_{eyear}_{var}.nc"
        )

        print(f"Done {var}")
