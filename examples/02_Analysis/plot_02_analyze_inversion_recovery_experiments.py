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

import numpy as np
import dnplab as dnp

file_name_path = "../../data/topspin/304"
data = dnp.load(file_name_path, assign_vdlist="t1")
data.attrs["experiment_type"] = "nmr_spectrum"

# %%
# Next, the FID has to be processed. Here, first a background correction is applied to remove any DC offsets, followed by windowing and FFT to get the frequency domain NMR spectrum.

data = dnp.remove_background(data)
data = dnp.apodize(data, lw=100)
data = dnp.fourier_transform(data)
data = dnp.autophase(data, method="manual", phase=90)

dnp.fancy_plot(data, xlim=[-50, 80], title="Inversion Recovery")
dnp.plt.show()


# %%
# Align Inversion Recovery Spectra
# --------------------------------
# The spectra shown here are the individual NMR spectra for different T1 recovery times. Note how these spectra are not perfectly aligned? This can be fixed using the aligning function. A more detailed description of the aligning routine of DNPLab is given in the tutorial :ref:`plot_01_align_nmr_spectra`.

data = dnp.ndalign(data)

dnp.fancy_plot(data, xlim=[-50, 80], title="Inversion Recovery, aligned")
dnp.plt.show()


# %%
# Integrate Spectra
# -----------------
# To determine the T1 relaxation time, we first have to integrate the peak intensity across the entire spectrum. After integration the workspace will have a new dnpdata object called "integrals" where the integral values and indirect axis are kept. To use ``fancy_plot`` the attribute experiment_type has to be changed to "inversion_recovery".

integrals = dnp.integrate(data)
integrals.attrs["experiment_type"] = "inversion_recovery"
dnp.fancy_plot(integrals)
dnp.plt.show()


# Depending on the quality of the data it is sometimes better to not integrate over the entire spectrum but to pick a peak region and only integrate over this region. To integrate over a region from 0 to 20 ppm us the following command:

integrals = dnp.integrate(data, regions=[(0, 20)])
integrals.attrs["experiment_type"] = "inversion_recovery"
dnp.fancy_plot(integrals)
dnp.plt.show()


# %%
# Fit Data
# --------
# To get the T1 value an inversion recovery function is fitted to the data sets. The fit requires an initial guess.

initial_guess = (2.0, -4000, 4000)  # initial guess for: T1, M_0, M_inf
out = dnp.fit(dnp.math.relaxation.t1, integrals, dim="t1", p0=initial_guess)

fit = out["fit"]
popt = out["popt"]
err = out["err"]

dnp.fancy_plot(integrals, title="Inversion Recovery")
dnp.plot(fit, "-")
dnp.plt.show()

T1 = popt["popt", 0]
M_0 = popt["popt", 1]
M_inf = popt["popt", 2]

print(T1.values)
print(M_0.values)
print(M_inf.values)
