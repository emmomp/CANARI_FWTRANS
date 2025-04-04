#!/usr/bin/env python3
# coding: utf-8
"""
Inputs for other scripts, mostly directory locations and metadata

Required to reproduce data for Boland et al. 2025 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr

FCNAME = "horflux_fw_denm"
CONV_DIR = "../data_out/denm_X_ECCOclimanom_hfreq"
EV_DIR = "../data_out/ev"
RECON_DIR = "../data_out/reconstructions"
EXPDIR = "/users/emmomp/data/canari/experiments"
GRIDDIR = "/users/emmomp/data/orchestra/grid2/"
CONTR_DIR = "../data_out/contrs"

TRANSP = "fw"
FCNAME = "horflux_fw_denm"

ecco_grid = xr.open_dataset("~/data/orchestra/other_data/ECCO_r3_alt/ECCOv4r3_grid.nc")

eyears = ["2006", "2014", "2000"]
adj_diag_map = {
    "adxx_qnet": ["oceQnet","EXFqnet" ],
    "adxx_tauu": ["oceTAUU", "EXFtauu"],
    "adxx_tauv": ["oceTAUV", "EXFtauv"],
    "adxx_empmr": [ "oceFWflx","EXFempmr"],
}

adj_units=dict(zip(['adxx_qnet','adxx_empmr','adxx_tauu','adxx_tauv'],['W/m$^2$','m/s','m$^2$/s','m$^2$/s']))
adj_labels=dict(zip(['adxx_qnet','adxx_empmr','adxx_tauu','adxx_tauv'],['Net Heat Flux','Net Freshwater Flux','Zonal Wind Stress','Meridional Wind Stress']))

masks_labels=dict(zip(['global','egland','natl','arct','gin'],['Global','E Gland','N Atl','Arctic','GIN']))

ecco_convs = {}
ecco_convs["all"] = []
ecco_convs["OCE"] = []
ecco_convs["EXF"] = []
ecco_convs_2d = {}
ecco_convs_2d["all"] = []
ecco_convs_2d["OCE"] = []
ecco_convs_2d["EXF"] = []
for adj, diag in adj_diag_map.items():
    ecco_convs["all"] = ecco_convs["all"] + [adj + "X" + v + "_sum" for v in diag]
    ecco_convs_2d["all"] = ecco_convs_2d["all"] + [adj + "X" + v for v in diag]
    ecco_convs["OCE"] = ecco_convs["OCE"] + [
        adj + "X" + v + "_sum" for v in diag if v[:3] == "oce"
    ]
    ecco_convs["EXF"] = ecco_convs["EXF"] + [
        adj + "X" + v + "_sum" for v in diag if v[:3] == "EXF"
    ]
    ecco_convs_2d["OCE"] = ecco_convs_2d["OCE"] + [
        adj + "X" + v  for v in diag if v[:3] == "oce"
    ]
    ecco_convs_2d["EXF"] = ecco_convs_2d["EXF"] + [
        adj + "X" + v  for v in diag if v[:3] == "EXF"
    ]

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
