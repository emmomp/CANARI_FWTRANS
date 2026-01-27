#!/usr/bin/env python3
# coding: utf-8
"""
generate_masks.py

Generate space and time masks for adjoint sensitivity experiments described in Boland et al. 2026, using code found in https://doi.org/10.5281/zenodo.17225253

Required to reproduce data for Boland et al. 2026 (in prep)
See https://github.com/emmomp/CANARI_FWTRANS for details

Updated Mar 2025

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import sys
import xmitgcm
import numpy as np
sys.path.insert(0, "/users/emmomp/Python/ECCOv4-py")
import ecco_v4_py as ecco
from inputs import MASK_DIR, ecco_grid


sections = [
    "Denmark Strait",
]
extra_metadata = xmitgcm.utils.get_extra_metadata(domain="llc", nx=90)

# Spatial masks
for section in sections:
    [section_pt1, section_pt2] = ecco.get_section_endpoints(section)
    maskC, maskW, maskS = ecco.get_section_line_masks(
        section_pt1, section_pt2, ecco_grid
    )
    masks_out = {"C": maskC, "W": maskW, "S": maskS}
    for mask in ["C", "W", "S"]:
        masktowrite = masks_out[mask].rename({"tile": "face"})
        facets = xmitgcm.utils.rebuild_llc_facets(masktowrite, extra_metadata)
        compact = xmitgcm.utils.llc_facets_2d_to_compact(facets, extra_metadata)
        xmitgcm.utils.write_to_binary(
            compact, f"{MASK_DIR}/{section.replace(' ','')}_mask{mask}"
        )

# Time masks
NY = 5
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
for mth in range(0, 12):
    maskT = np.zeros(NY * 12 + 6)
    maskT[(NY - 1) * 12 + mth] = 1
    xmitgcm.utils.write_to_binary(maskT, f"{MASK_DIR}/{NY}y_{all_mths[mth]}_maskT")
