# %%
"""
.. _08-calculate-dnp-enhancements-II:

====================================
08 - Calculate DNP Enhancements (II)
====================================

This is a slightly more advanaced example demonstrating how to ingegrate multiple regions of an NMR spectrum and calculate the DNP enhancement for two different regions.
"""
# %%
# Load NMR Spectra
# ----------------
# In this example, we will load the aligned NMR spectra saved in one of the previous examples (:ref:`07_align_nmr_spectra`) and save it in the workspace (ws). This now includes the objects "proc", "NotAligned", and "Aligned", as well as the "raw" data.

import dnplab as dnp

file_name_path = "../data/h5/PowerBuildUpAligned.h5"
ws = dnp.dnpImport.load(file_name_path)

# %%
# Integrate Calculate DNP Enhancement Factors
# -------------------------------------------
# Without the integral region the calculate_enhancement function will calculate the enhancement over the entire spectrum
# First lets take a look at the data that is in the workspace

print(ws.keys())

# %%
# We can see that there are four keys in the workspace: "Aligned", "NotAligned"," proc", and "raw". 
# %%
# Now integrate First define the integral region


ws.processing_buffer = "Aligned"

dnp.dnpTools.integrate(ws, integrate_center = [2, 7], integrate_width = [1, 1])

ws.processing_buffer = "proc"


print(ws.keys())



# dnp.dnpNMR.calculate_enhancement(ws, off_spectrum = 0)

# # %%
# # .. note::
# #     The default behavior of the ``calculate_enhancement`` function is to use the first spectrum as the Off signal. If this is the case, the argument ``off_spectrum`` is not necessary unless you want to specify the slice that contains the off spectrum.
# #     The ``calculate_enhancement``` function can also calculate the enhancement for specific regions of the spectrum. THis behavior will be discussed in the next example (:ref:`07_align_nmr_spectra`).

# # %%
# # If needed, access your array of enhancements as:
# enhancements = ws["enhancement"].values

# # %%
# # Plot Enhancement Data
# # ---------------------
# # Finally, we can plot the enhancement data versus the microwave power.

ws["integrals"].reorder(["power"])

dnp.dnpResults.figure()
# dnp.dnpResults.plot(ws["enhancement"],linestyle = '-', marker = 'o', fillstyle = 'none')
dnp.dnpResults.plot(ws["integrals"])
# dnp.dnpResults.plt.xlabel("Microwave Power (dBm)")
# dnp.dnpResults.plt.ylabel("ODNP Enhancement Factor")
dnp.dnpResults.plt.title("10 mM TEMPO in Toluene")
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()
