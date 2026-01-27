#!/usr/bin/env python3
# coding: utf-8
"""
calc_pert_fwanoms.py

Calculate surface forcing anomalies from ECCOv4r4 perturbation experiments

Requires:
- Diagnostic from files from perturbation experiments to ECCOv4r4, 
  using the xx_tauu.0000000129.data.{pert}_{pert_sign} files
- xx_tauu.0000000129.data.{pert}_{pert_sign}, in ../other_data, were created using make_perts.py

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jan 2026

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import glob
import os
import sys
import xmitgcm
import xarray as xr

sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import GRIDDIR, ecco_grid

exf_vars = [
    "EXFqnet",
    "EXFempmr",
    "EXFtaux",
    "EXFtauy",
    "oceTAUX",
    "oceTAUY",
    "oceFWflx",
    "oceQnet",
]

RHOCONST = 1029


def load_pert(pert, dfreq):
    print(f"Loading {dfreq} {pert} data")
    exf_ds = []
    for evar in exf_vars:
        print(evar)
        try:
            var_ds = xmitgcm.open_mdsdataset(
                data_dir=f"{pert}/diags/{evar}_{dfreq}_mean",
                grid_dir=GRIDDIR,
                prefix=[f"{evar}_{dfreq}_mean"],
                geometry="llc",
                delta_t=3600,
                ref_date=startdate,
                read_grid=False,
            )
            var_ds = var_ds.rename({"face": "tile"})
            exf_ds.append(var_ds)
        except:
            print(f"No {dfreq} {evar}")
    if not exf_ds:
        print(f"No {dfreq} at all for {pert}")
        return None
    exf_ds = xr.merge(exf_ds)
    return exf_ds


all_pert_plus = glob.glob("../experiments/pert_*plus*")
for pert_plus in all_pert_plus:
    pert_minus = pert_plus.replace("plus", "minus")
    pert_lab = pert_plus.split("/")[-1].replace("plus", "")

    if "2002" in pert_plus:
        year_start = "2002"
    else:
        year_start = "1996"
    startdate = f"{year_start}-01-01"

    for freq in ["week", "mon"]:
        fout = f"../data_out/perts/{pert_lab}_{freq}_pertfields.nc"
        if os.path.isfile(fout):
            print(f"Found {fout}, skipping")
        else:

            ds_plus = load_pert(pert_plus, freq)
            ds_minus = load_pert(pert_minus, freq)
            if (not ds_plus) and (not ds_minus):
                continue
            ds_ctrl = load_pert(f"../experiments/pert_10y_ctrl_{year_start}", freq)

            ds_all = []
            if ds_plus and ds_ctrl:
                ds_pert_plus = (ds_plus - ds_ctrl).assign_coords({"pert": "plus"})
                ds_all.append(ds_pert_plus)
            if ds_minus and ds_ctrl:
                ds_pert_minus = (ds_minus - ds_ctrl).assign_coords({"pert": "minus"})
                ds_all.append(ds_pert_minus)
            if ds_plus and ds_minus and ds_ctrl:
                ds_pert_lin = ((ds_plus - ds_minus) * 0.5).assign_coords(
                    {"pert": "lin"}
                )
                ds_pert_nonlin = ((ds_plus + ds_minus) * 0.5 - ds_ctrl).assign_coords(
                    {"pert": "nonlin"}
                )
                ds_all.append(ds_pert_lin)
                ds_all.append(ds_pert_nonlin)

            ds_all = xr.concat(ds_all, "pert")
            if ds_all:
                if "EXFtaux" in ds_all:
                    ds_all["EXFtaux"] = ds_all["EXFtaux"].rename({"i": "i_g"}).load()
                    ds_all["EXFtauy"] = ds_all["EXFtauy"].rename({"j": "j_g"}).load()
                    [ds_all["EXFtauu"], ds_all["EXFtauv"]] = (
                        ecco.vector_calc.UEVNfromUXVY(
                            ds_all["EXFtaux"], ds_all["EXFtauy"], ecco_grid
                        )
                    )
                    ds_all.drop_vars(["EXFtaux", "EXFtauy"])

                if "oceTAUX" in ds_all:
                    ds_all["oceTAUX"] = ds_all["oceTAUX"].load()
                    ds_all["oceTAUY"] = ds_all["oceTAUY"].load()
                    [ds_all["oceTAUU"], ds_all["oceTAUV"]] = (
                        ecco.vector_calc.UEVNfromUXVY(
                            ds_all["oceTAUX"], ds_all["oceTAUY"], ecco_grid
                        )
                    )
                    ds_all.drop_vars(["oceTAUX", "oceTAUY"])

                if "oceFWflx" in ds_all:
                    ds_all["oceFWflx"] = (
                        -ds_all["oceFWflx"] / RHOCONST
                    )  # Convert to m/s, same sign as EMPMR

                if "oceQnet" in ds_all:
                    ds_all["oceQnet"] = -ds_all["oceQnet"]  # same sign as EXF

                print(f"Writing {fout}")
                ds_all.to_netcdf(fout)
