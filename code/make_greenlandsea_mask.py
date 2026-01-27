#!/usr/bin/env python3
# coding: utf-8
"""
make_greenlandsea_mask.py

Generate custom mask of the Greenland Sea

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Jan 2026

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import numpy as np
from cartopy.io.shapereader import Reader
from shapely import Point

ecco_grid = xr.open_dataset("~/data/orchestra/other_data/ECCO_r3_alt/ECCOv4r3_grid.nc")

shp=Reader('../other_data/iho_greenlandsea/iho.shp')
geoms=shp.geometries()
geom=next(geoms)

gland_mask=xr.full_like(ecco_grid.rA,False,dtype=np.bool)

for tile in [2,6]:
    hfac=ecco_grid.hFacC.isel(k=0,tile=tile)
    for i in range(0,90):
        for j in range(0,90):
            point=hfac.sel(i=i,j=j)
            if point>0:
                this_point = Point(point.XC.data,point.YC.data)
                res = geom.contains(this_point)
                gland_mask.loc[dict(tile=tile,i=i,j=j)]=res
            else:
                gland_mask.loc[dict(tile=tile,i=i,j=j)]=False

gland_mask.to_netcdf('../other_data/greenlandsea_mask.nc')              