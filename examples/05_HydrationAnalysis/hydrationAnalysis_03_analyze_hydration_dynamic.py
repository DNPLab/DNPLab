# %%
"""
.. _hydrationAnalysis_02_analyze_hydration_dynamic:

==============================================================================
Analyze Hydration Dynamic from ODNP-enhanced NMR and T1 Inversion-Recovery Data
==============================================================================

This example demonstrates how to loaded and analyze the hydration dynamic of data from h5 file

The sample is 1 mM TEMPO in water

Loaded from h5 file, the dictionary must have: 

    1. hydration_info['odnp_enhancements']: the enhancements of NMR spectra at different microwave power
    2. hydration_info['enhancement_fit_coefficients']: the coefficients and the errors of enhancements vs. power fitting curve
    3. hydration_info['odnp_ir_coefficients']: the coefficients and the errors of ODNP_enhanced T1 inversion recovery integral
    4. hydration_info['ir_coefficients']: the coefficients and the errors of T1 inversion recovery integrals curve (no radical)

The hydration analysis results are DNPData objects and stored as a dictionary saved in a h5 file:
    1. hydration_info['t1_coefficients']: the coefficients of t1 vs. microwave power fitting curve
    2. hydration_info['f']: the leakage factor and its error
    3. hydration_info['ksigma_sp_coefficients']: the coefficients of building ksigma*s(p)
    4. hydration_info['ksigma_sp_errors']: the errors of building ksigma*s(p)
    5. hydration_info['ksigma_smax']: the ksigma*smax and its error
    6. hydration_info['krho']: the krho and its error
    7. hydration_info['xi_smax']: the coupling factor*smax and its error
    
    If smax is defined in hydration_info['sample_information']:

    8. hydration_info['xi']: the coupling factor and its error
    9. hydration_info['ksigma']: the ksigma and its error
"""

# %%
# First you need to prepare the Python environment by importing DNPLab. In this example we also have to imp, also numpy is imported. However, this is only necessary if you would like to access numpy functions and procedures directly.
import dnplab as dnp
import numpy as np

# The sample is defined for plotting
sampleTag = '1 mM TEMPO in Water'

# %%
# Hydration Analysis
# After you have all required DNPData objects in h5 file, the following function will analyze the data and return the hydration results
# The updating to h5 file is not on by default
# If you want to save the hydration analysis results to h5 file, make save_file True

hydration_results = dnp.analysis.hydration.hydration_analysis('../../data/h5/1mM_TEMPO_water_hydration_info.h5', save_file = True) # The function can read the h5 file and output the hydration analysis results

# The hydration results will be shown automatically.