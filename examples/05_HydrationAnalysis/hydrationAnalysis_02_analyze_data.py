# %%
"""
.. _hydrationAnalysis_02_analyze_data:

==============================================================================
Analyze Processed ODNP-enhanced NMR and T1 Inversion-Recovery Experiments Data
==============================================================================

This example demonstrates how to loaded the processed data from h5 file, analyze the data and save the analysis results to h5 file for hydration analysis.

The sample is 1 mM TEMPO in water

Loaded from h5 file, the dictionary must have: 

    1. hydration_info['odnp_proc']: processed ODNP_enhanced NMR spectra
    2. hydration_info['odnp_ir_proc']: processed ODNP_enhanced T1 Inversion Recovery spectra
    3. hydration_info['odnp_ir_proc']: processed T1 Inversion Recovery spectra (no radical)
    4. hydration_info['sample_information']: the sample information and must include 'radical_concentration'

The analysis results are DNPData objects and stored as a dictionary saved in a h5 file for hydration analysis:
    1. hydration_info['odnp_enhancements']: the enhancements of NMR spectra at different microwave power
    2. hydration_info['t1_integrals']: the integrals of T1 Inversion Recovery spectra at different microwave power
    3. hydration_info['t10_integrals']: the integrals of T1 Inversion Recovery spectrum (no radical)
    4. hydration_info['enhancement_fit_coefficients']: the coefficients and the errors of enhancements vs. power fitting curve
    5. hydration_info['odnp_ir_coefficients']: the coefficients and the errors of ODNP_enhanced T1 inversion recovery integrals curves
    6. hydration_info['ir_coefficients']: the coefficients and the errors of T1 inversion recovery integrals curve (no radical)
"""

# %%
# First you need to prepare the Python environment by importing DNPLab. In this example we also have to imp, also numpy is imported. However, this is only necessary if you would like to access numpy functions and procedures directly.
import dnplab as dnp
import numpy as np

# The sample is defined for plotting
sampleTag = '1 mM TEMPO in Water'

# %%
# Load h5 file and assign the dictionary to 'hydration_info'

hydration_info = dnp.load('../../data/prospa/1mM_TEMPO_Water/h5/hydration_info.h5')

# %%
# Analyze ODNP_enhanced NMR Spectra

enh_data = hydration_info['odnp_proc'] # access the processed ODNP_enhanced NMR spectra
enhancements = dnp.integrate(enh_data, dim = 'f2', regions = [(4.2, 5.0)]) # integrate spectra
enhancements.values /= enhancements.values[0] # Normalize the first enhancement to 1
hydration_info['odnp_enhancements'] = enhancements # store enhancements in dictionary
enh_out = dnp.fit(dnp.relaxation.buildup_function, enhancements, dim = 'power', p0 = (-40, 0.5)) # calculate the coefficients for enhancements fitting curve
enh_fit = enh_out['fit'] # fitting curve
enh_coefficients = enh_out['popt'] # fitting coefficients, [E_max, p_1/2]
enh_errors = enh_out['err'] # fitting errors, [E_max error, p_1/2 error]
hydration_info['odnp_fit_coefficients'] = enh_coefficients # store coefficients in a dictionary
hydration_info['odnp_fit_errors'] = enh_errors # store errors in a dictionary

# %%
# Plot Enhancements vs. Microwave Power

dnp.plt.figure()
# dnp.fancy_plot(enh_fit, title = sampleTag)
dnp.plot(enhancements, ls = 'none', marker = 'o', color = '#F37021')
dnp.plot(enh_fit, ls = '-', color = '#F37021', label = '$E_{max}$ = %0.03f +/- %0.03f' %(enh_coefficients.values[0], enh_errors.values[0]))
dnp.plt.xlabel('Power (W)')
dnp.plt.ylabel('Enhancements (a.u.)')
dnp.plt.title('Enhancements vs. Power')
dnp.plt.grid(':')
dnp.plt.legend()
dnp.plt.show()

# %%
# Analyze ODNP_enhanced T1 Inversion Recovery Spectra

t1_data = hydration_info['odnp_ir_proc'] # access the processed T1 Inversion Recovery spectra
t1_integrals = dnp.integrate(t1_data, dim = 'f2', regions = [(4.2, 5.0)]) # integrate spectra
t1_integrals.values /= t1_integrals.values[-1] # Normalize the last integrals to 1
hydration_info['t1_integrals'] = t1_integrals # store t1 integrals in dictionary
t1_out = dnp.fit(dnp.relaxation.t1, t1_integrals, dim = 't1', p0 = (2, -1.5, 1)) # calculate the coefficients for building inversion recovery curve
t1_fit = t1_out['fit'] # fitting curve
t1_coefficients = t1_out['popt'] # t1 inversion recovery coefficients, [T1, M_0, M_inf]
t1_errors = t1_out['err'] # t1 inversion recovery errors, [T1 error, M_0 error, M_inf error]
hydration_info['odnp_ir_coefficients'] = t1_coefficients # store coefficients in a dictionary
hydration_info['odnp_ir_errors'] = t1_errors # store errors in a dictionary
# %%
# Plot T1 Inversion Recovery at Different Microwave Power

for m, p in enumerate(t1_data.coords['power']):
    dnp.plt.figure('Inversion Recovery, MW = %0.01f W' %p)
    dnp.plot(t1_integrals['power', m], ls = 'none', marker = 'o', color = '#F37021')
    dnp.plot(t1_fit['power', m], ls = '-', color = '#F37021', label = '$T_1$ = %0.03f +/- %0.03f s' %(t1_coefficients.values[0,m], t1_errors.values[0,m]))
    dnp.plt.xlabel('Time (s)')
    dnp.plt.ylabel('Signal (a.u.)')
    dnp.plt.title('Inversion Recovery\nMW = %0.01f W' %p)
    dnp.plt.grid(':')
    dnp.plt.legend()
dnp.plt.show()

# %%
# Analyze T1 Inversion Recovery Spectra (No Radical)

t10_data = hydration_info['ir_proc'] # access the processed T1 Inversion Recovery spectra
t10_integrals = dnp.integrate(t10_data, dim = 'f2', regions = [(4.2, 5.0)]) # integrate spectra
t10_integrals.values /= t10_integrals.values[-1] # Normalize the last integrals to 1
hydration_info['t10_integrals'] = t10_integrals # store t10 integrals in dictionary
t10_out = dnp.fit(dnp.relaxation.t1, t10_integrals, dim = 't1', p0 = (2, -1.5, 1)) # calculate the coefficients for building inversion recovery curve
t10_fit = t10_out['fit'] # fitting curve
t10_coefficients = t10_out['popt'] # t10 inversion recovery coefficients, [T10, M_0, M_inf]
t10_errors = t10_out['err'] # t10 inversion recovery errors, [T10 error, M_0 error, M_inf error]
hydration_info['ir_coefficients'] = t10_coefficients # store coefficients in a dictionary
hydration_info['ir_errors'] = t10_errors # store errors in a dictionary

# %%
# Plot T1 Inversion Recovery (No Radical)

dnp.plt.figure()
dnp.plot(t10_integrals, ls = 'none', marker = 'o', color = '#F37021')
dnp.plot(t10_fit, ls = '-', color = '#F37021', label = '$T_1$ = %0.03f +/- %0.03f s' %(t10_coefficients.values[0], t10_errors.values[0]))
dnp.plt.xlabel('Time (s)')
dnp.plt.ylabel('Signal (a.u.)')
dnp.plt.title('Inversion Recovery (No Radical)')
dnp.plt.grid(':')
dnp.plt.legend()
dnp.plt.show()

# %%
# Saving h5 file
dnp.save(hydration_info, '../../data/prospa/1mM_TEMPO_Water/h5/hydration_info.h5', overwrite = True)