# %% [markdown]
"""
Calculate DNP enhancement
=========================

This example demonstrates how to import DNP-NMR data and calculate enhancement
"""
# %%
import numpy as np
import dnplab as dnp

# %%

# %% [markdown]
# The enhancement function can use 2D datasets where one spectrum is the "Off" signal and the
# remainder are the "On" signals. Since none of our example data are of this nature, below we
# will import individual 1D spectra sequentially and calculate the enhancements. But first,
# if you have 2D data use the code:
data = dnp.dnpImport.load("..", data_type="")
ws = dnp.create_workspace()
ws.add("raw", data)
ws.copy("raw", "proc")

dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)

# %% [markdown]
# With a 2D set ws["proc"].values, calculate the enhancements all at once with off_spectrum=1
# specifying that the first slice of your 2D array is the off spectrum. **Note** This is the
# default behavior and the argument off_spectrum is not necessary unless the off spectrum is
# not the first slice.
dnp.dnpNMR.calculate_enhancement(ws, off_spectrum=1)
# %%

# %% [markdown]
# If needed, access your array of enhancements as:
enhancements = ws["enhancement"].values
# %%

# %% [markdown]
# **Example of handling collection of separate 1D spectra**
# ----------------------------------------------------------
#
# If you have several 1D spectra to compare, you could do the following:
# Start by importing and processing the 1D "Off" spectrum:
data = dnp.dnpImport.load("../data/prospa/toluene_10mM_Tempone/1", data_type="prospa")
ws_off = dnp.create_workspace()
ws_off.add("rawOff", data)
ws_off.copy("rawOff", "proc")

dnp.dnpNMR.remove_offset(ws_off)
dnp.dnpNMR.window(ws_off, linewidth=10)
dnp.dnpNMR.fourier_transform(ws_off, zero_fill_factor=2)

# %% [markdown]
# Make a list of file paths for your 1D "On" spectra
filenames = [
    "../data/prospa/toluene_10mM_Tempone/29",
    "../data/prospa/toluene_10mM_Tempone/30",
    "../data/prospa/toluene_10mM_Tempone/31",
]
# %%

# %% [markdown]
# Create a workspace for the "On" spectra
ws_on = dnp.create_workspace()
# %%

# %% [markdown]
# This loop imports each spectrum from the above list and calculates enhancement for each,
# appending them to an array called 'enhancements'.
#
# Create the empty array to append to, then loop over the list of filenames processing and
# calculating the enhancement for each, then appending the enhancement value to the array.
enhancements = []
for index, path in enumerate(filenames):
    data = dnp.dnpImport.load(path, data_type="prospa")  # import
    ws_on.add("proc", data)  # add directly to the processing buffer "proc"

    dnp.dnpNMR.remove_offset(ws_on)
    dnp.dnpNMR.window(ws_on, linewidth=10)
    dnp.dnpNMR.fourier_transform(ws_on, zero_fill_factor=2)

    dnp.dnpNMR.calculate_enhancement(ws_on, off_spectrum=ws_off, on_spectra=ws_on)

    enhancements.append(ws_on["enhancement"].values)
# %%

# %% [markdown]
# This will add the object ["enhancement"] to the workspace that is the first argument.
# You may add to any workspace by putting it in the above code in place of ws_on.
# %%
