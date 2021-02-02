"""
Calculate DNP enhancement
============================

This example demonstrates how to import DNP-NMR data and calculate enhancement
"""

# Import any needed standard packages
import numpy as np

# Import DNPLab as dnp
import dnplab as dnp

# the enhancement function can use 2D datasets where one spectrum is the "Off" signal and the remainder are the "On" signal. Since none of our example data are of this nature, below we will import individual 1D spectra sequentially and calculate the enhancements. But first, if you have 2D data use the code:

# find and import some 2D data
data = dnp.dnpImport.load("..", data_type="")
# create a workspace and name it "ws":
ws = dnp.create_workspace()
ws.add("raw", data)
# add the raw data to the processing buffer "proc"
ws.copy("raw", "proc")

# do some basic processing
dnp.dnpNMR.remove_offset(ws)
dnp.dnpNMR.window(ws, linewidth=10)
dnp.dnpNMR.fourier_transform(ws, zero_fill_factor=2)

# With a 2D set ws["proc"].values, calculate the enhancements all at once with off_spectrum=1 specifying that the first slice of your 2D array is the off spectrum. **Note** This is the default behavior and the argument off_spectrum is not necessary unless the off spectrum is not the first slice.
dnp.dnpNMR.calculate_enhancement(ws, off_spectrum=1)

# access your array of enhancements as:
enhancements = ws["enhancements"].values

#################################################################################################
# If you have several 1D spectra to compare, you could do the following:

# Start by importing the 1D "Off" spectrum:
data = dnp.dnpImport.load("../data/prospa/toluene_10mM_Tempone/1", data_type="prospa")
ws_off = dnp.create_workspace()
ws_off.add("rawOff", data)
# copy the raw object to the processing buffer "proc"
ws_off.copy("rawOff", "proc")

# do some basic processing
dnp.dnpNMR.remove_offset(ws_off)
dnp.dnpNMR.window(ws_off, linewidth=10)
dnp.dnpNMR.fourier_transform(ws_off, zero_fill_factor=2)

# Make a list of file paths for your 1D "On" spectra
filenames = [
    "../data/prospa/toluene_10mM_Tempone/29",
    "../data/prospa/toluene_10mM_Tempone/30",
    "../data/prospa/toluene_10mM_Tempone/31",
]

# create a workspace for the "On" spectra
ws_on = dnp.create_workspace()

# this loop imports each spectrum from the above list and calculates enhancement for each, appending them to an array called 'enhancements'

# create empty array to append to
enhancements = []
for index, path in enumerate(filenames):
    # import
    data = dnp.dnpImport.load(path, data_type="prospa")
    # add directly to the processing buffer "proc"
    ws_on.add("proc", data)
    # process each "On" spectrum
    dnp.dnpNMR.remove_offset(ws_on)
    dnp.dnpNMR.window(ws_on, linewidth=10)
    dnp.dnpNMR.fourier_transform(ws_on, zero_fill_factor=2)
    # calculate enhancement and add to the ws_on["enhancement"] object. **Note** when entered this way you must choose the ws in which to add the object "enhancement". For example, to add the enhancements to the ws_off workspace, the first argument should be ws_off rather than ws_on.
    dnp.dnpNMR.calculate_enhancement(ws_on, off_spectrum=ws_off, on_spectra=ws_on)
    # append the enhancement to the array "enhancements"
    enhancements.append(ws_on["enhancement"].values)
