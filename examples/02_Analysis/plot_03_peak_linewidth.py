# %%
"""
.. _plot_03_peak_linewidth:

==============
Peak Linewidth
==============

This example demonstrates how to determine the linewidth of peaks in an NMR spectrum.
"""
# %%
# Load and Process Data
# ---------------------
# Let's start with importing a NMR data set and apply some standard processing. Here we use an ODNP-enhanced low-frequency (14.5 MHz) NMR spectrum of 10 mM TEMPONE in toluene. The spectrum can be found in the example data folder of DNPLab.

import numpy as np
import dnplab as dnp

file_name_path = "../../data/prospa/toluene_10mM_Tempone/40"
data = dnp.load(file_name_path)
data.attrs["experiment_type"] = "nmr_spectrum"

data = dnp.apodize(data, lw=2)
data = dnp.fourier_transform(data)
data = dnp.phase(data, p0=175)

dnp.plt.figure()
dnp.fancy_plot(data, xlim=[-20, 20])
dnp.plt.title("ODNP-Enhanced NMR Spectrum of 10 mM TEMPONE in Toluene")
dnp.plt.show()

# %%
# Get Peak Linewidth
# ---------------------
# Next we apply the DNPLab function *find_peaks* to get a list of all the peaks the algorithm can find in the spectrum and print the list in the terminal.

peaks = dnp.find_peaks(data)
print(peaks.values)

# %%
# The function returns a dnpdata object with a 3xn arrays of values , with n being the number of peaks found in the spectrum, corresponding to the columns of the array. In this examples, the peak picking function finds two peaks, one for the methyl group and one for all aromatic protons. When printing the values of the dnpdata object the first row corresponds to the index of the peak in the spectrum, the linewidth in Hz is given in the second row, and the relative peak height is given in the third row. The values can be printed in a human readable form using the *peak_info* function

dnp.peak_info(peaks)

# %%
# Setting the Threshold
# ---------------------
# By default, *find_peaks* will identify every feature in the spectrum with a amplitude > 5 % of the maximum intensity as a peak. For this, the spectrum is first normalized to a maximum amplitude of 1. This default value can be change as shown in the next line

peaks = dnp.find_peaks(data, threshold=0.6)
dnp.peak_info(peaks)

# %%
# Now, only peaks with a relative intensity > 0.6 are included in the peak list.
