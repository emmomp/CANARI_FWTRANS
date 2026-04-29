# CANARI_FWTRANS
This repositary contains python code and notebooks to accompany the manuscript "Local and Remote Drivers of Liquid Freshwater Transport through Denmark Strait" Boland et al 2026. The contents will allow for the reproduction of all tables and figures in the paper, as well the reproduction of the data files necessary for the figures. See below for more details.

Feel free to use or reproduce the code and figures but please attribute as outlined in the license.

For more details on the ECCOv4r4 ocean state estimate, used to produce the data in this study, see https://ecco-group.org/products-ECCO-V4r4.htm

Emma Boland Jan 2026 [emmomp@bas.ac.uk](mailto:emmomp@bas.ac.uk)

## Requirements

Third-party requirements are found notebooks/requirements.txt and code/requirements.txt. Additionally, two further requirements are:
- My forked version of [ECCOv4-py](https://github.com/emmomp/ECCOv4-py/tree/fw_transports) which includes the ability to calculate freshwater transports from model diagnostics as well as some updates to plotting tools.
- [xadjoint](https://github.com/emmomp/xadjoint), a tool for reading and analysing adjoint sensitivity experiments from MITgcm, which automatically generates time step and lagged time metadata, based on xmitgcm.

## Steps to reproduce the paper's figures and tables

To reproduce the paper's figures and tables, follow these steps:
- Download the data required for the figures from the Figshare repository [https://doi.org/10.6084/m9.figshare.31169446]. Alternatively this data can be re-generated from the original model output using the python files in the [code](code/) directory - see "Steps to reproduce the paper's analysis from model output".
- Install necessary libraries (see requirements above or [figure_notebooks/requirements.txt](figure_notebooks/requirements.txt)).
- Clone the [figure_notebooks](figure_notebooks/) directory into the same directory that contains 'data_in'.
- Run the notebooks.

## Steps to reproduce the paper's analysis from model output

To reproduce the data files required to produce the figures, follow these steps:
- Rerun the paper's experiments (see "Steps to rerun the model experiments")
- Install the necessary python libraries (see [code/requirements.txt](code/requirements.txt))
- Run the python scripts in [code](code/)

## Steps to rerun the model experiments

Take my published version of [ECCOv4r4](https://github.com/emmomp/ECCO-v4-Configurations/releases/tag/v1.1), 
using code_noparam to compile and namelist_adjsen when running, with masks generated using code/generate_masks.py 
