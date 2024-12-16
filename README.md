# DFTB3-CH-CC-RPs_For_HydrideTransfers

This repository contains the data needed to replicate/evaluate our research work titled: "Adapted DFTB3 repulsive potentials reach DFT accuracy for hydride transfer reactions in enzymes"

If you need more technical information about our research, don't hesitate to contact us via email at josevlibera2010@gmail.com.

## Repository Structure

### DFTB3_Fitted_Parameters
Folder containing the C-C and C-H fitted parameters to describe the hydride transfer reaction in Ccr.

### IRC
Folder containing the output of the IRC calculations describing Ccr's hydride transfer reaction in vacuo.

### String_Calculations
Folder containing the inputs and results from the Adaptive String Method's (ASM) calculations on Ccr and DHFR enzymes' hydride transfer reactions.

#### Ccr/
Contains all calculations related to the Crotonyl-CoA Carboxylase/Reductase enzyme:
- **Ccr_DFT/**: DFT-level calculations for the hydride transfer reaction
- **Ccr_DFTB3_CH_CC_mod3/**: Calculations using modified DFTB3 parameters (CH_4_CC_3)
- **Ccr_DFTB3_CH_mod5/**: Calculations using modified DFTB3 parameters (CH_5)
- **Ccr_DFTB3_Original_parameters/**: Calculations using original DFTB3 parameters
- **INPUTS/**: Input files for ASM calculations on Ccr hydride transfer reaction step

#### DHFR/
Contains all calculations related to the Dihydrofolate Reductase enzyme:
- **DHFR_DFTB3_CH_CC_mod3/**: Calculations using modified DFTB3 parameters (CH_4_CC_3)
- **DHFR_DFTB3_Original_parameters/**: Calculations using original DFTB3 parameters
- **INPUTS/**: Input files for ASM calculations on DHFR hydride transfer reaction

### Parametrization
Contains the scripts and data used for DFTB3 parameter optimization:
- **LIB/**: Library folder containing the python clases in utils.py module
- **PARAM/**: Folder with DFTB3's original parameter files and for storing the modified parameters 
- **RUNS/**: Calculation outputs using sqm program from AmberTools suite
- **STRUCT/**: Ordered IRC trajectory in xyz format 
- **AnalysisRepPot-DistDistrib-IRC.ipynb**: Jupyter notebook for analyzing repulsive potentials and distance distributions from IRC calculations
- **FitHrm_CH_EXP1_deltaE.ipynb**: Jupyter notebook for fitting C-H using linear combinations of harmonic functions
- **FitHrm_CH_CC_EXP1_deltaE.ipynb**: Jupyter notebook for fitting C-H and C-C RPs using linear combinations of harmonic functions
- **PlotGaussAndHarmonics.ipynb**: Jupyter notebook for visualizing Gaussian and harmonic functions


## Citation
If you use these parameters or data in your research, please cite our paper: (https://doi.org/10.26434/chemrxiv-2024-kkwwg-v2)
