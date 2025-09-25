# CANARI_FWTRANS
This repositary contains python code and notebooks to accompany the manuscript "Local and Remote Drivers of Liquid Freshwater Export through Denmark Strait" Boland et al 2025. The contents will allow for the reproduction of all tables and figures in the paper, as well the reproduction of the data files necessary for the figures. See below for more details.

Feel free to use or reproduce the code and figures but please attribute as outlined in the license.

For more details on the ECCOv4r4 ocean state estimate, used to produce the data in this study, see https://ecco-group.org/products-ECCO-V4r4.htm

Emma Boland Sep 2025 [emmomp@bas.ac.uk](mailto:emmomp@bas.ac.uk)

## Requirements

Third-part requirements are found notebooks/requirements.txt and code/requirements.txt, which includes a requirement for my forked version of [ECCOv4-py](https://github.com/emmomp/ECCOv4-py.git@24b62eafc276d962ab2f78affc3fabdf773b09ca)

## Steps to reproduce the paper's figures and tables

To reproduce the paper's figures and tables, follow these steps:
- Download the data required for the figures from the Figshare repository XXX and place in a directory named 'data_out'. Alternatively this data can be re-generated from the original model output using the python files in the [code](code/) directory - see "Steps to reproduce the paper's analysis from model output".
- Install necessary libraries (see requirements above or [figure_notebooks/requirements.txt](figure_notebooks/requirements.txt)).
- Clone the [figure_notebooks](figure_notebooks/) directory into the same directory that contains 'data_in'.
- Run the notebooks.

## Steps to reproduce the paper's analysis from model output

To reproduce the data files required to produce the figures, follow these steps:

## Steps to rerun the model experiments
