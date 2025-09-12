#!/usr/bin/env python
# coding: utf-8
"""
calc_eccov4_climanoms.py

Code to calculate climatological anomalies of ECCOv4r4 weekly mean ocean forcing fields, 
generated from diagnostic fields found in {EXPDIR}/fwd_26y/diags/.

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Sep 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import sys
import xmitgcm
import xarray as xr

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import EXPDIR, GRIDDIR, ecco_grid, oce_vars

attrs = {
    "contact": "emmomp@bas.ac.uk",
    "references": "ECCOv4r4 data from Boland et al (in prep)",
    "date": "Created on " + date.today().strftime("%d/%m/%Y"),
    "notes": "Data produced by analysis of the ECCOv4r4 state estimate, see ecco-group.org",
}

exf_ds = []
STARTDATE = "1992-01-01"
RHOCONST = 1029

print("Loading data")
for evar in oce_vars[:1]:
    print(evar)
    var_ds = xmitgcm.open_mdsdataset(
        data_dir=EXPDIR + "/fwd_26y/diags/",
        grid_dir=GRIDDIR,
        prefix=[f"{evar}_week_mean"],
        geometry="llc",
        delta_t=3600,
        ref_date=STARTDATE,
        read_grid=False,
    )
    var_ds = var_ds.rename({"face": "tile"})
    exf_ds.append(var_ds)
exf_ds = xr.merge(exf_ds)

if "oceTAUX" in oce_vars:
    exf_ds["oceTAUX"] = exf_ds["oceTAUX"].load()
    exf_ds["oceTAUY"] = exf_ds["oceTAUY"].load()
    [exf_ds["oceTAUU"], exf_ds["oceTAUV"]] = ecco.vector_calc.UEVNfromUXVY(
        exf_ds["oceTAUX"], exf_ds["oceTAUY"], ecco_grid
    )
    exf_ds.drop_vars(["oceTAUX", "oceTAUY"])

if "oceFWflx" in oce_vars:
    exf_ds["oceFWflx"] = (
        -exf_ds["oceFWflx"] / RHOCONST
    )  # Convert to m/s, same sign as EMPMR

if "oceQnet" in oce_vars:
    exf_ds["oceQnet"] = -exf_ds["oceQnet"]  # same sign as EXF

print("Calculating climatology and anomalies, saving")
ds_clim = exf_ds.groupby(exf_ds.time.dt.month).mean(dim="time").load()
ds_climanom = exf_ds.groupby(exf_ds.time.dt.month) - ds_clim
ds_clim.attrs.update(attrs)
ds_climanom.attrs.update(attrs)
ds_climanom.to_netcdf(f"{EXPDIR}/fwd_26y/exf_climanoms.nc")
ds_clim.to_netcdf(f"{EXPDIR}/fwd_26y/exf_clim.nc")
print("Done!")
