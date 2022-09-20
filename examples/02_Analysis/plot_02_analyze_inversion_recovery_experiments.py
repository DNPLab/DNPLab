# %%
"""
.. _plot_02_analyse_inversion_recovery_experiments:

=========================================
Analyze T1 Inversion-Recovery Experiments
=========================================

This example demonstrates how to import TopSpin data from an inversion recovery NMR experiment and determine the T1 relaxation rate through a fit.
"""
# %%
# Load Inversion Recovery Spectra
# -------------------------------
# Start with importing data and creating the DNPLab workspace. Note that the data are added to the workspace as 'raw' then copied to the processing buffer 'proc' to preserve the raw data during processing.
# from numpy import size

import numpy as np
import dnplab as dnp

file_name_path = "../../data/topspin/304"
data = dnp.load(file_name_path, assign_vdlist='t1')

# data = dnp.load(file_name_path, verbose=True)
# values = dnp.io.topspin.topspin_vdlist(file_name_path)
# data.coords['t1'] = values


data.attrs["experiment_type"] = "nmr_spectrum"

# %%
# Next, we process the FIDs to obtain the frequency domain NMR spectrum

data = dnp.remove_background(data)

# dnp.remove_background(data)
# data = dnp.apodize(data, lw=100)
# data = dnp.fourier_transform(data)

# data = dnp.autophase(data, method="manual", phase=90)

# dnp.fancy_plot(data, xlim=[-50, 80], title="Inversion Recovery")
# dnp.plt.show()


# %%
# Align Inversion Recovery Spectra
# --------------------------------
# The spectra shown here are the individual NMR spectra for different evolution times. Note how the spectra are not perfectly aligned? This can be fixed using the aligning function. A more detailed description of the aligning routine of DNPLab is given in the tutorial :ref:`plot_01_align_nmr_spectra`.

# data = dnp.align(data)

# dnp.fancy_plot(data, xlim=[-50, 80], title="Inversion Recovery, aligned")
# dnp.plt.show()


# %%
# Integrate Spectra
# -----------------
# To determine the T1 relaxation time, we first have to integrate the peak intensity across the entire spectrum. After integration the workspace will have a new dnpdata object called "integrals" where the integral values and indirect axis are kept. To use ``fancy_plot`` the attribute experiment_type has to be changed to "inversion_recovery".


# dnp.plt.figure()
# dnp.waterfall(data, dx = 50, dy = 100000)
# dnp.plt.plot(data.values[:,7].real)
# dnp.plt.show()
# #


# dnp.plt.figure()
# dnp.plt.plot(data.values[:, 1].real)
# dnp.plt.show()

# integrals = dnp.integrate(data)

# print(integrals.values[:,0].real )
# print(integrals)


# integrals.attrs["experiment_type"] = "inversion_recovery"

# dnp.fancy_plot(integrals)
# dnp.plt.show()

# initial_guess = (2., -4000, 4000) # initial guess for: T1, M_0, M_inf

# out = dnp.fit(dnp.math.relaxation.t1, integrals, dim = 't1', p0 = initial_guess)

# fit = out['fit']
# popt = out['popt']
# err = out['err']

# dnp.fancy_plot(integrals)
# dnp.plot(fit, '-')
# print(popt)
# T1 = popt['popt',0]
# M_0 = popt['popt',1]
# M_inf = popt['popt',2]

# print(T1.values)

# dnp.plt.show()
