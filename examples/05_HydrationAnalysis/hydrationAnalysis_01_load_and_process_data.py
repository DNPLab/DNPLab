# %%
"""
.. _hydrationAnalysis_01_load_and_process_data:

===============================================================
Load ODNP-enhanced NMR and T1 Inversion-Recovery Experiments
===============================================================

This example demonstrates how to load, process and save the dataset that is used to perform hydration analysis of Prospa data from ODNP-enhanced NMR and T1 Inversion-Recovery Experiments through a function

The sample is 1 mM TEMPO in MeCN and water mixture, chi = 0.25

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

# %%
# Load And Process ODNP_enhanced NMR Experiment Data
# --------------------------------------------------
# In this example we use 1D NMR spectra acquired using Prospa. Example data is located in the data folder. 
# The dataset is splitted into 2: average of 4 and average of 16 because the sample has poor signal-to-noise level at low microwave power.
# The Load function cannot process the dataset with different dimension, so, in this example, the sub-datasets are processed separately and merged after summing the 'Average' dimension.
# If the dataset has the unchanging dimension, no separate dataset processing is required.

# sub_dataset: Average of 16
enh_file_avg16 = [
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/1',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/2',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/3',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/4',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/5',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg16_20220412/6',
]

# Microwave Power for this sub_dataset
enh_power_avg16 = np.r_[-99., 0., 5., 10., 15., 20.,] # dBm
enh_power_w_avg16 = 10.**(enh_power_avg16/10.)/1000. # W

enh_data_avg16 = dnp.load(enh_file_avg16, dim = 'power', coord = enh_power_w_avg16) # load raw data and create 'power' dimension
hydration_info['odnp_raw_avg16'] = enh_data_avg16 # store raw data in dictionary
enh_data_avg16 = dnp.apodize(enh_data_avg16, dim = 't2', kind = 'lorentz_gauss', exp_lw = 4., gauss = 8.) # you can change apodize method and linewidth if necessary
enh_data_avg16 = dnp.fourier_transform(enh_data_avg16, dim = 't2') # fourier_transform
enh_data_avg16 = dnp.phase_cycle(enh_data_avg16, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
enh_data_avg16 = dnp.ndalign(enh_data_avg16) # aligning spectra using correlation method
# enh_data_avg16 = dnp.autophase(enh_data_avg16) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
enh_phase_avg16 = [0.980, 0.984, 0.980, 0.992, 0.985, 0.985]
for n, theta in enumerate(enh_phase_avg16):
    enh_data_avg16['power', n] *= np.exp(1j * 180/ np.pi * theta)

enh_data_avg16 = enh_data_avg16.sum('Average')/16. # averaging the 'Average' dimension
enh_data_avg16.coords['f2'] -= enh_data_avg16.abs.argmax('f2').values[0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['odnp_proc'] = enh_data_avg16 # store processed dataset in dictionary

# %%
# The following processing for the sub_dataset may be optional if the dimension is consistent.
# sub_dataset: Average of 4
enh_file_avg4 = [
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/1',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/2',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/3',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/4',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/5',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/6',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/7',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/8',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/9',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/1Pulse_avg4_20220412/10',
]

# Microwave Power for this sub_dataset
enh_power_avg4 = np.r_[25., 26., 27., 28., 29., 30., 31., 32., 33., 25.] # dBm
enh_power_avg4 = 10.**(enh_power_avg4/10.)/1000. # W

enh_data_avg4 = dnp.load(enh_file_avg4, dim = 'power', coord = enh_power_w_avg16) # load raw data and create 'power' dimension
hydration_info['odnp_raw_avg4'] = enh_data_avg4 # store raw data in dictionary
enh_data_avg4 = dnp.apodize(enh_data_avg4, dim = 't2', kind = 'lorentz_gauss', exp_lw = 4., gauss = 8.) # you can change apodize method and linewidth if necessary
enh_data_avg4 = dnp.fourier_transform(enh_data_avg4, dim = 't2') # fourier_transform
enh_data_avg4 = dnp.phase_cycle(enh_data_avg4, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
enh_data_avg4 = dnp.ndalign(enh_data_avg4) # aligning spectra using correlation method
# enh_power_avg4 = dnp.autophase(enh_power_avg4) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
enh_phase_avg4 = [0.986, 0.986, 0.985, 0.985, 0.987, 0.986, 0.985, 0.985, 0.983, 0.981]
for n, theta in enumerate(enh_phase_avg4):
    enh_data_avg4['power', n] *= np.exp(1j * 180/ np.pi * theta)

enh_data_avg4 = enh_data_avg4.sum('Average')/4. # averaging the 'Average' dimension
enh_data = enh_data_avg16.copy()
enh_data.concatenate(enh_data_avg4, dim = 'power') # combine two sub_datasets
enh_data = dnp.ndalign(enh_data) # aligning spectra again
enh_data.coords['f2'] -= enh_data.abs.argmax('f2').values[0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['odnp_proc'] = enh_data # store processed dataset in dictionary

# %%
# # Load And Process ODNP_enhanced T1 Inversion Recovery Experiment
# --------------------------------------------------
# The processing steps are similar with that of ODNP-enhanced NMR Experiment Data
t1_file = [
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/1',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/2',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/3',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/4',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/5',
    '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_20220412/6',
]

# Microwave Power for T1 Inversion Recovery Dataset
t1_power = np.r_[-99., 26., 29., 30.8, 32., 33.] # dBm
t1_power_w = 10.**(t1_power/10.)/1000. # W

t1_data = dnp.load(t1_file, dim = 'power', coord = t1_power_w) # load raw data and create 'power' dimension
t1_data = dnp.apodize(t1_data, dim = 't2', kind = 'lorentz_gauss', exp_lw = 4., gauss = 8.) # you can change apodize method and linewidth if necessary
t1_data = dnp.fourier_transform(t1_data, dim = 't2') # fourier_transform
t1_data = dnp.phase_cycle(t1_data, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
t1_data = dnp.ndalign(t1_data, reference = t1_data['Power', 1, 't1', 9, 'Average', 0].sum('Power').sum('t1').sum('Average')) # aligning spectra using correlation method, the reference is optional
# t1_data = dnp.autophase(t1_data) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
t1_phase = np.loadtxt('../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/t1_phase.csv', delimiter = ',', skiprows = 1)
for n, theta in enumerate(t1_phase):
    t1_data['power', n] *= np.exp(1j * 180/ np.pi * theta)

t1_data = t1_data.sum('Average') # averaging the 'Average' dimension
t1_data.coords['f2'] -= t1_data.abs.argmax('f2').values[0][0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['odnp_ir_proc'] = t1_data # store processed dataset in dictionary

# %%
# Load And Process T1 Inversion Recovery Experiment Data (no radical), so as called T10
# --------------------------------------------------
# The processing steps are similar with that of ODNP-enhanced NMR Experiment Data

t10_file = '../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/T1-IR-FID_no_radical_20220411/1/data.3d' # path to file data.3d file directly

t10_data = dnp.load(t10_file, dim = 'power', coord = t1_power_w) # load raw data and create 'power' dimension
t10_data = dnp.apodize(t10_data, dim = 't2', kind = 'lorentz_gauss', exp_lw = 4., gauss = 8.) # you can change apodize method and linewidth if necessary
t10_data = dnp.fourier_transform(t10_data, dim = 't2') # fourier_transform
t10_data = dnp.phase_cycle(t10_data, dim = 'Average', receiver_phase = [0,1,2,3]) # phase_cycle for Prospa dataset
t10_data = dnp.ndalign(t10_data)
# t10_data = dnp.autophase(t10_data) # auto phase the spectra so that the signal is negative maximum

# If autophase is not working, define phases and phase spectra manually
t10_phase = [1.086, 1.092, 1.085, 1.082, 1.082, 1.082, 1.080, 1.080, 1.080, 1.080]
for n, theta in enumerate(t10_phase):
    t10_data['power', n] *= np.exp(1j * 180/ np.pi * theta)

# If ndalign is not working properly because of poor signal-to-noise, align the spectra manually based on the following script
values = t10_data['t1', 1, 'Average', 0:2].values
values = np.roll(values, -23)
t10_data['t1', 1, 'Average', 0:2] = values

t10_data = t10_data.sum('Average') # averaging the 'Average' dimension
t10_data.coords['f2'] -= t10_data.abs.argmax('f2').values[0] - 4.6 # moving the the maximum spectra to water's chemical reference 4.6 ppm
hydration_info['ir_proc'] = t10_data # store processed dataset in dictionary

# %% 
# To perform hydration analysis, radical concentration (M) is required in h5 file
# Other experiment details are recommended added to the dictionary

hydration_info['experiment_details'] = {'radical_concentration': 1e-3, 'sample': '1 mM TEMPOL in MeCN/water, chi = 0.25'}

# %%
# Saving h5 file
dnp.save('../../data/1mM_TEMPOL_in_MeCN_Water_Mixture_Chi_0.25/h5/hydration_info.h5', overwrite = True)