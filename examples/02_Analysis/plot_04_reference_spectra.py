# %%
"""
.. _plot_04_reference_spectra:

==================
Reference  Spectra
==================

This example demonstrates how to reference spectra using Toluene as an example.
"""
# %%
# Load Toluene Spectra
# -------------------------------
# Start with importing data and creating the DNPLab workspace.
# The data has been processed and saved in .h5 file.

import numpy as np
import dnplab as dnp

sampleTag = "10 mM TEMPO in Toluene"
# file_name_path = "../../data/h5/PowerBuildUpAligned.h5"
# data = dnp.load(file_name_path)

# %%
# Plot Toluene Spectra
# --------------------
# Once the data are imported and processed, it is time to plot the 1D spectra.

# dnp.fancy_plot(data, xlim=[-20, 20], title=sampleTag)
# dnp.plt.show()


# %%
# Find Toluene Peaks
# ------------------
# In the previous example, we have demonstrated how to find peaks. Now, we apply the DNPLab function *find_peaks* to get a DNPData object of all peaks.


# peaks = dnp.find_peaks(data)
# dnp.peak_info(peaks)

# %%
# Reference Proton Peak
# ---------------------
# We can see that chemical shift for proton peak is roughly at 2.79 ppm.
# Let's reference this peak to 7.70 ppm applying the DNPLab function *reference*.

# data = dnp.reference(data, old_ref = 2.79, new_ref = 7.70)

# We can plot the spectra after reference.
# dnp.fancy_plot(data, xlim=[-10, 30], title=sampleTag)
# dnp.plt.show()

# %%
# Check the Proton Peak
# ---------------------
# After we referenced the peaks, we can apply *find_peaks* functions again to get a DNPData object of all peaks.

# peaks = dnp.find_peaks(data)
# dnp.peak_info(peaks)
