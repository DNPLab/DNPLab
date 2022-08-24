# %%
"""
.. _hydrationAnalysis_01_load_and_process_data:

========================================================================
Load and Process ODNP-enhanced NMR and T1 Inversion-Recovery Experiments
========================================================================

This example demonstrates how to load, process and save the dataset that is used to perform hydration analysis of Prospa data from ODNP-enhanced NMR and T1 Inversion-Recovery Experiments through a function

The sample is 1 mM TEMPO in water

To perform the hydration analysis, the following datasets are required:

    1. ODNP_enhanced NMR Experiment (NMR Spectra vs. Microwave Power)
    2. ODNP_enhanced T1 Inversion Recovery Experiment (T1 Inversion Recovery NMR Spectra vs. Microwave Power)
    3. T1 Inversion Recovery Experiment (T1 Inversion Recovery NMR Spectra without radical)

They are processed to DNPData objects and stored as a dictionary saved in a h5 file for hydration analysis.
"""

# %%
# First you need to prepare the Python environment by importing DNPLab. In this example we also have to imp, also numpy is imported. However, this is only necessary if you would like to access numpy functions and procedures directly.
import dnplab as dnp
import numpy as np

# Second you need to create dictionary to store DNPData objects
hydration_info = {}

# The sample is defined for plotting
sampleTag = '1 mM TEMPO in Water'

# %%
# Load And Process ODNP_enhanced NMR Experiment Data
# --------------------------------------------------
# In this example we use 1D NMR spectra acquired using Prospa. Example data is located in the data folder. 

enh_file = [
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/1/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/2/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/3/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/4/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/5/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/6/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/7/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/8/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/9/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/10/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/11/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/12/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/13/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/14/data.2d',
    '../../data/prospa/1mM_TEMPO_Water/1Pulse_20220516/15/data.2d',
]

# Microwave Power for this sub_dataset
enh_power = np.r_[-99., 0., 5., 10., 15., 20., 25., 26., 27., 28., 29., 30., 31., 32., 33.] # dBm
enh_power_w = 10.**(enh_power/10.)/1000. # W

enh_data = dnp.load(enh_file, dim = 'power', coord = enh_power_w) # load raw data and create 'power' dimension
hydration_info['odnp_raw'] = enh_data # store raw data in dictionary
enh_data = dnp.apodize(enh_data, dim = 't2', lw = 4.) # you can change apodize method and linewidth if necessary
enh_data = dnp.fourier_transform(enh_data, dim = 't2') # fourier_transform
enh_data = dnp.phase_cycle(enh_data, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
enh_data = dnp.ndalign(enh_data) # aligning spectra using correlation method
enh_data = dnp.autophase(enh_data) # auto phase the spectra so that the signal is negative maximum

enh_data = enh_data.sum('Average')/16. # averaging the 'Average' dimension
enh_data.coords['f2'] -= enh_data.abs.argmax('f2').values[0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['odnp_proc'] = enh_data # store processed dataset in dictionary

# %%
# Plot processed ODNP_enhanced NMR Spectra

enh_data.attrs['experiment_type'] = 'nmr_spectrum'
dnp.plt.figure()
dnp.fancy_plot(enh_data, xlim = [20, -10], title = sampleTag)
dnp.plt.show()

# %%
# Load And Process ODNP_enhanced T1 Inversion Recovery Experiment
# --------------------------------------------------
# The processing steps are similar with that of ODNP-enhanced NMR Experiment Data
t1_file = [
    '../../data/prospa/1mM_TEMPO_Water/T1-IR-FID_20220516/1/data.3d',
    '../../data/prospa/1mM_TEMPO_Water/T1-IR-FID_20220516/2/data.3d',
    '../../data/prospa/1mM_TEMPO_Water/T1-IR-FID_20220516/3/data.3d',
]

# Microwave Power for T1 Inversion Recovery Dataset
t1_power = np.r_[27., 30., 31.8] # dBm
t1_power_w = 10.**(t1_power/10.)/1000. # W

t1_data = dnp.load(t1_file, dim = 'power', coord = t1_power_w) # load raw data and create 'power' dimension
t1_data = dnp.apodize(t1_data, dim = 't2', lw = 4.) # you can change apodize method and linewidth if necessary
t1_data = dnp.fourier_transform(t1_data, dim = 't2') # fourier_transform
t1_data = dnp.phase_cycle(t1_data, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
t1_data = dnp.ndalign(t1_data) # aligning spectra using correlation method, the reference is optional
# t1_data = dnp.autophase(t1_data) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
t1_phase = np.loadtxt('../../data/prospa/1mM_TEMPO_Water/t1_phase.csv', delimiter = ',', skiprows = 1)
for m, p in enumerate(t1_power_w):
    for n, theta in enumerate(t1_phase[:, m+1]):
        t1_data['power', m, 't1', n] *= np.exp(-1j * 180/ np.pi * theta)

t1_data = t1_data.sum('Average') # averaging the 'Average' dimension
t1_data.coords['f2'] -= t1_data.abs.argmax('f2').values[0][0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['odnp_ir_proc'] = t1_data # store processed dataset in dictionary

# %%
# Plot processed ODNP_enhanced T1 Inversion Recovery Spectra

for m, p in enumerate(t1_power_w):
    dnp.plt.figure('Processing Result, MW = %0.01f' %p)
    dnp.plot(t1_data['power', m])
    dnp.plt.xlabel('Chemical Shift (ppm)')
    dnp.plt.ylabel('Signal (a.u.)')
    dnp.plt.title('Results\nMW = %0.01f W' %p)
    dnp.plt.xlim(20.,-10.)
    dnp.plt.grid(':')
    dnp.plt.title(sampleTag)
dnp.plt.show()

# %% 
# Load And Process T1 Inversion Recovery Experiment Data (no radical), so as called T10
# --------------------------------------------------
# The processing steps are similar with that of ODNP-enhanced NMR Experiment Data

t10_file = '../../data/prospa/1mM_TEMPO_Water/T1-IR-FID_no_radical_20220429/1/data.3d' # path to file data.3d file directly
t10_data = dnp.load(t10_file, dim = 'power', coord = t1_power_w) # load raw data and create 'power' dimension
t10_data = dnp.apodize(t10_data, dim = 't2', lw = 4.) # you can change apodize method and linewidth if necessary
t10_data = dnp.fourier_transform(t10_data, dim = 't2') # fourier_transform
t10_data = dnp.phase_cycle(t10_data, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
t10_data = dnp.ndalign(t10_data)
# t10_data = dnp.autophase(t10_data) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
t10_phase = [0.035, 0.035, 0.040, 0.040, 0.035, 0.030, 0.030, 0.023, 0.026, 0.023, 0.023, 0.023, 0.020, 0.016, 0.014]
for n, theta in enumerate(t10_phase):
    t10_data['t1', n] *= np.exp(1j * 180/ np.pi * theta)

t10_data = t10_data.sum('Average') # averaging the 'Average' dimension
t10_data.coords['f2'] -= t10_data.abs.argmax('f2').values[0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['ir_proc'] = t10_data # store processed dataset in dictionary

# %%
# Plot processed T1 Inversion Recovery Spectrum
t10_data.attrs['experiment_type'] = 'nmr_spectrum'
dnp.plt.figure()
dnp.fancy_plot(t10_data, xlim = [20, -10], title = sampleTag)
dnp.plt.show()

# %% 
# To perform hydration analysis, radical concentration (M) is required in h5 file
# Other experiment details are recommended added to the dictionary

hydration_info['sample_information'] = {'radical_concentration': 1e-3, 'smax': 0.39, 'sample': '1 mM TEMPOL in MeCN/water, chi = 0.25'}

# %%
# Saving h5 file
dnp.save(hydration_info, '../../data/prospa/1mM_TEMPO_Water/h5/hydration_info.h5', overwrite = True)