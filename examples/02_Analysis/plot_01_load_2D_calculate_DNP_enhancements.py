# %%
"""
.. _plot_01_load_2D_calculate_DNP_enhancements:

===================================================
Load a 2D dnpdata object and calculate enhancements
===================================================

This example demonstrates how to import DNP-NMR data in form of a 2D dnpdata object from an hdf5 file, calculate the DNP enhancement factors and plot the enhancment vs. the applied microwave power. This example uses the 2D data object that was created in a previous tutorial (:ref:`plot_04_create_dnpdata_object_from_individual_files`).

"""
# %%
# Load NMR Spectra
# ----------------
# In this example, we will calculate the DNP enhancement factor for each individual DNP spectrum and create a figure showing the DNP enhancement versus the applied microwave power. For this, we will import the 2D dnpdata object created in the previous sample. If you are not yet familiar with how to concatenate individual spectra into the 2D danpdata object, check out this tutorial: :ref:`plot_04_create_dnpdata_object_from_individual_files`.
#
# First, load the 2D dnplab data object:

import dnplab as dnp

# file_name_path = "../data/h5/PowerBuildUp.h5"
# ws = dnp.dnpImport.load(file_name_path)

# %%
# Calculate DNP Enhancement Factors
# ---------------------------------
# DNPLab provides a convenient way to calculate the DNP enhancement factors by using the ``calculate_enhancement`` function. Enhancement factors are calculated using integrals. Integrals can be calculated over the entire spectrum, multiple regions, or can be just a single point. However, without calculating integrals first, the ``calculate_enhancement`` function will return an error.

# dnp.dnpTools.integrate(ws)
# dnp.dnpNMR.calculate_enhancement(ws)

# %%
# In this case, the integral is calculated over the entire spectrum followed by calculating the enhancement factors.

# %%
# .. note::
#     The default behavior of the ``calculate_enhancement`` function is to use the first spectrum as the Off signal. If this is the case, the argument ``off_spectrum`` is not necessary. If you want to specify a particular spectrum that contains the off signal, use the ``off_spectrum`` argument.

# #     The ``calculate_enhancement``` function can also calculate the enhancement for specific regions of the spectrum. THis behavior will be discussed in the next example (:ref:`07_align_nmr_spectra`).


# # %%
# # If needed, access your array of enhancements as:
# enhancements = ws["enhancements"].values

# # %%
# # Plot Enhancement Data
# # ---------------------
# # Finally, we can plot the enhancement data versus the microwave power.

# dnp.dnpResults.figure()
# dnp.dnpResults.plot(ws["enhancements"],linestyle = '-', marker = 'o', fillstyle = 'none')
# dnp.dnpResults.plt.xlabel("Microwave Power (dBm)")
# dnp.dnpResults.plt.ylabel("ODNP Enhancement Factor")
# dnp.dnpResults.plt.title("10 mM TEMPO in Toluene")
# dnp.dnpResults.plt.grid(True)
# dnp.dnpResults.show()
