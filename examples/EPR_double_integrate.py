# %% [markup]
"""
DNPLab EPR data tutorial
========================

This example demonstrates how to use DNPLab to load and double integrate EPR data.

"""
# %%
import dnplab as dnp

# %% [markup]
# Use some 1D EPR data
data = dnp.dnpImport.load("../data/bes3t/1D_CW.DTA")
# %%

# %% [markup]
# Create a workspace and add the loaded data
ws = dnp.create_workspace()
ws.add("raw", data)
# %%

# %% [markup]
# Copy the data to the processing buffer 'proc'
ws.copy("raw", "proc")
# %%

# %% [markup]
# Raw data are now in both ws["proc"] and ws["raw"]. The following processing steps
# will modify ws["proc"] but leave ws["raw"] untouched so that you can always return
# to the original data.
#
# To see what attributes are loaded with the data use:
print("BEFORE PROCESSING: " + str(ws["raw"].attrs.keys()))
# %%

# %% [markup]
# Rename the x axis to NMR style notation to more seamlessly interact with some
# functions. "f2" is the default for argument dim.
ws["proc"].rename("G", "f2")
# %%

# %% [markup]
# Check the current axis names with
print("raw axis names: " + str(ws["raw"].dims))
print("proc axis names: " + str(ws["proc"].dims))
# %%

# %% [markup]
# Baseline correct with a zeroth order polynomial to remove the offset
dnp.dnpTools.baseline(ws, type="polynomial", order=0)
# %%

# %% [markup]
# Keep a copy of the processed spectrum before integrating
ws.copy("proc", "spec")
# %%

# %% [markup]
# Double integrate
dnp.dnpTools.integrate(ws, type="double")
# %%

# %% [markup]
# If needed, access your processed spectrum as follows:
x_axis = ws["raw"].coords["G"]  # imported field axis
spectrum = ws["spec"].values  # spectrum after any processing steps, before integration
first_integral = ws["proc"].attrs["first_integral"]  # first integral
double_integral = ws["proc"].values  # double integral value
# %%

# %% [markup]
# Your originally imported spectrum is ws["raw"].values.
#
# To see if any attributes have been added after processing
print("AFTER PROCESSING: " + str(ws["proc"].attrs.keys()))
# %%

# %% [markup]
# Plot the EPR spectrum using dnpResults
dnp.dnpResults.plot(ws["raw"], label="EPR Spectrum")
dnp.dnpResults.legend()
dnp.dnpResults.show()
# %%

# %% [markup]
# Hint: Use matplotlib to plot arrays such as the attributes "first_integral" or
# "baseline" over the spectrum, in this case access the spectrum directly as '.values'
import matplotlib.pyplot as plt

plt.figure
plt.plot(ws["raw"].coords["G"], ws["spec"].values, label="corrected EPR Spectrum")
plt.plot(ws["raw"].coords["G"], ws["raw"].values, label="original EPR Spectrum")
plt.plot(ws["raw"].coords["G"], ws["proc"].attrs["baseline"], label="Baseline")
plt.plot(
    ws["raw"].coords["G"], ws["proc"].attrs["first_integral"], label="First Integral"
)
plt.legend()
plt.show()
# %%
