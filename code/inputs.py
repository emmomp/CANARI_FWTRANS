#!/usr/bin/env python3
# coding: utf-8
"""
Inputs for other scripts, mostly directory locations and metadata

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr

CONV_DIR = "../data_out/denm_X_ECCOclimanom_hfreq"
PERTCONV_DIR = "../data_out/denm_X_perts_hfreq"
EV_DIR = "../data_out/ev"
RECON_DIR = "../data_out/reconstructions"
EXPDIR = "/users/emmomp/data/canari/experiments"
GRIDDIR = "/users/emmomp/data/orchestra/grid2/"
SOLN_DIR='/data/expose/ECCOv4-r4/Version4/Release4/nctiles_monthly'
CONTR_DIR = "../data_out/contrs"
DATA_DIR= "../data_out"
MASK_DIR = "/data/smurphs/emmomp/canari/masks"

TRANSP = "fw"
FCNAME = "horflux_fw_denm"

ecco_grid=xr.open_dataset('/data/expose/ECCOv4-r4/Version4/Release4/nctiles_grid/ECCO-GRID.nc') 

eyears = ["2000", "2006", "2014"]
adj_diag_map = {
    "adxx_qnet": "oceQnet",
    "adxx_tauu": "oceTAUU",
    "adxx_tauv": "oceTAUV",
    "adxx_empmr":  "oceFWflx",
}

adj_units=dict(zip(['adxx_qnet','adxx_empmr','adxx_tauu','adxx_tauv'],['W/m$^2$','m/s','m$^2$/s','m$^2$/s']))
adj_labels=dict(zip(['adxx_qnet','adxx_empmr','adxx_tauu','adxx_tauv'],['Net Heat Flux','Net Freshwater Flux','Zonal Wind Stress','Meridional Wind Stress']))

oce_vars=["oceQnet","oceTAUU", "oceTAUV", "oceFWflx","wind_OCE"]
exf_units=dict(zip(oce_vars+['EXFtauu',],['W/m$^2$','m$^2$/s','m$^2$/s','m/s','m$^2$/s','m$^2$/s']))
exf_labels=dict(zip(oce_vars,['Net Heat Flux','Zonal Wind Stress','Meridional Wind Stress','Net Freshwater Flux','Wind Stress']))

masks_labels=dict(zip(['global','egland','natl','arct','gin','hudson', 'north', 'baffin','barents','gland','norw'],['Global','E Gland','N Atlantic','Arctic','GIN','Hudson','North','Baffin','Barents','Greenland','Norwegian']))

ecco_convs =  []
ecco_convs_2d =  []
for adj, diag in adj_diag_map.items():
    ecco_convs = ecco_convs + [adj + "X" + diag + "_sum",]
    ecco_convs_2d = ecco_convs_2d + [adj + "X" + diag,]

oexps = [
    "transfw_Mar_noparam_7d",
    "transfw_Jun_noparam_7d",
    "transfw_Sep_noparam_7d",
    "transfw_Dec_noparam_7d",
]
mth = [
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
]
mthi = dict(zip(mth, list(range(1, 13))))

imth = dict(zip(range(1, 13), mth))
