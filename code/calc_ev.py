#!/usr/bin/env python
# coding: utf-8
"""
calc_ev.py

Code to calculate the explained variance of a reconstruction.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Feb 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import os.path
from datetime import date
import xarray as xr
import utils as ut
from inputs import *

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "Explained variances of ECCOv4r4 Denmark Strait FW flux \
    reconstructions from Boland et al. 2025 (inprep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 global ocean state \
    estimate,see ecco-group.org",
}

# Lags for 2D EV fields
lags = [None, -0.25, -0.5, -1.5, -4]
lag_labels = [
    "0 to -3m lag",
    "-3m to -6m lag",
    "-6m to -18m lag",
    "-18m to -4 y lag",
    "0 to -4 y lag",
]

EXPT = "fwd_26y"
fc = ut.get_soln(FCNAME, f"{EXPDIR}/{EXPT}")
fc_climanom, fc_mth = ut.soln_anoms(fc)
    
def main():
    out={}
    
    for eyear in eyears:
        print(eyear)
        out[eyear]={}

        conv_ecco,_,_ = ut.load_ecco_convs(CONV_DIR, eyear)

        print("Calculating global EV time series")

        dJpred_var_bylag = conv_ecco[ecco_convs["all"]].squeeze().chunk('auto')

        fout = f"{EV_DIR}/{eyear}/{FCNAME}_ev_bylag_byvar_{eyear}.nc"
        out[eyear]['ev_total'] = calc_ev_1d(fout, dJpred_var_bylag)
        fout = f"{EV_DIR}/{eyear}/{FCNAME}_ev_bylag_byvar_bymonth_{eyear}.nc"
        out[eyear]['ev_total_bym'] = calc_ev_1d(fout, dJpred_var_bylag, bymonth=True)

        dJpred_ecco_cumsum = (
             dJpred_var_bylag.cumsum("lag_years").assign_coords(conv_ecco.coords)
        )

        fout = f"{EV_DIR}/{eyear}/{FCNAME}_cumev_bylag_byvar_{eyear}.nc"
        out[eyear]['cum_ev_total'] = calc_ev_1d(fout, dJpred_ecco_cumsum)
        fout = f"{EV_DIR}/{eyear}/{FCNAME}_cumev_bylag_byvar_bymonth_{eyear}.nc"
        out[eyear]['cum_ev_total_bym'] = calc_ev_1d(fout, dJpred_ecco_cumsum, bymonth=True)

        print("calculating 2-D EV")
        out[eyear]['ev_2d_full']={}
        out[eyear]['ev_2d_bym']={}
        for var in ecco_convs_2d:
            print(var)
            dJ_2d = conv_ecco[var]

            fout = f"{EV_DIR}/{FCNAME}_fullEV2d_{var}_{eyear}.nc"
            out[eyear]['ev_2d_full'][var] = calc_ev_2d(fout, dJ_2d, lags, lag_labels)
            fout = f"{EV_DIR}/{FCNAME}_monthlyEV2d_{var}_{eyear}.nc"
            out[eyear]['ev_2d_bym'][var] = calc_ev_2d(fout, dJ_2d, lags, lag_labels, bymonth=True)
        
    return out

            
def calc_ev_1d(file, dJ, bymonth=False):
    """
    Calculates expected variance time series of provided reconstructions and writes to file
    
    Parameters
    ----------
    file : str
        Location of file to save output. If file already exists, nothing happens
    dJ : xarray dataset
        Contains reconstructions of objective function
    by_month : logical, optional
        Default False. If true, calculates ev of monthly time series, otherwise full time series

    Returns
    -------
    ev : xarray dataset
        Contains explained variance of dJ tested against fc_climanom or fc_mth

    """
    if os.path.isfile(file):
        print(f"Found {file}, skipping")
        return None
    dJ = dJ.chunk({"year": -1}).sel(year=slice(1996, None))
    if bymonth:
        dJ = dJ.swap_dims({"exp": "month"})
        fc = fc_mth.sel(year=slice(1996, None))
        vardim = "year"
    else:
        dJ = dJ.stack(yearexp=["exp", "year"]).swap_dims({"yearexp": "dates"})
        fc = fc_climanom.sel(time=slice("1996-01-01", None)).rename({"time": "dates"})
        vardim = "dates"

    dJ["wind_EXF"] = dJ["adxx_tauuXEXFtauu_sum"] + dJ["adxx_tauvXEXFtauv_sum"]
    dJ["wind_OCE"] = dJ["adxx_tauuXoceTAUU_sum"] + dJ["adxx_tauvXoceTAUV_sum"]
    dJ["all_EXF"] = dJ[ecco_convs["EXF"]].to_array().sum("variable")
    dJ["all_OCE"] = dJ[ecco_convs["OCE"]].to_array().sum("variable")

    ev = 1 - (fc - dJ).var(vardim) / fc.var(vardim)
    print(f"Writing to {file}")
    ev.attrs.update(attrs)
    ev.to_netcdf(file)
    return ev

def calc_ev_2d(file, dJ, lag_list, laglabels, bymonth=False):
    """
    Calculates 2D expected variance of provided contribution fields, integrated over given lags, writes to file
    
    Parameters
    ----------
    file : str
        Location of file to save output. If file already exists, nothing happens
    dJ : xarray dataset
        Contains 2D reconstructiosn of objective function
    lag_list : list of lists
        Contains pairs of lags to denote limits of lag integrals
    laglabels : list of str
        Contains strings to provide coordinate labels for lags, should be same length as lag_list
    by_month : logical, optional
        Default False. If true, calculates ev of monthly time series, otherwise full time series

    Returns
    -------
    ev_out : xarray dataset
        Contains explained variance of dJ contributions tested against fc_climanom or fc_mth

    """
    if os.path.isfile(file):
        print(f"Found {file}, skipping")
        return None
    ev_out = []
    dJ = dJ.chunk({"year": -1}).sel(year=slice(1996, None))
    if bymonth:
        dJ = dJ.swap_dims({"exp": "month"})
        fc = fc_mth.sel(year=slice(1996, None))
        vardim = "year"
    else:
        dJ = dJ.stack(yearexp=["exp", "year"]).swap_dims({"yearexp": "dates"})
        fc = fc_climanom.sel(time=slice("1996-01-01", None)).rename({"time": "dates"})
        vardim = "dates"

    for ilag in range(0, len(lag_list) - 1):
        dJ_lag = dJ.sel(lag_years=slice(lag_list[ilag], lag_list[ilag + 1])).sum(
            "lag_years"
        )
        ev_2d = (1 - (fc - dJ_lag).var(vardim) / fc.var(vardim)).assign_coords(
            lag=laglabels[ilag]
        )
        ev_out.append(ev_2d / ecco_grid.rA)
    ilag += 1
    dJ_lag = dJ.sel(lag_years=slice(None, lag_list[-1])).sum("lag_years")
    ev_2d = (1 - (fc - dJ_lag).var(vardim) / fc.var(vardim)).assign_coords(
        lag=laglabels[ilag]
    )
    ev_out.append(ev_2d / ecco_grid.rA)

    ev_out = xr.concat(ev_out, "lag")
    ev_out.attrs.update(attrs)
    ev_out.to_netcdf(file)
    return ev_out

if __name__=='__main__':
    main()
    