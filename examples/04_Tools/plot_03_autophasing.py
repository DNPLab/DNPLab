# %%
"""
.. _plot_03_autophasing:

=========================================================
Load a 2D dnpdata object and show the autphasing function
=========================================================

This example demonstrates how to import DNP-NMR data in form of a 2D dnpdata object from ovnmrj folder and use the autophasing function to have a fixed phase relationship between the spectra.

"""
# %%
# Load NMR Spectra
# ----------------
# Some 2D ODNP data is loaded from the examples folder
#
# First, load the ovnmrj data into a DNPData object, by loading all folders. This is faciliated by using the standard library pathlib
# The data is located relative to the examples -> 04_Tools folder in the dnplab structure
import dnplab as dnp
import pathlib

basicFolder = pathlib.Path("../../data/vnmrj/autophase_example_data")

exp_number = lambda x: int(x.split("_")[-1].split(".")[0])
files = [str(k) for k in basicFolder.glob("Proton*") if k.is_dir()]
ind = [exp_number(k) for k in files]


files = sorted(files, key=exp_number)
experiment_number = [k for k in range(len(files))]
dnpData = dnp.load(files, dim="experiment_number", coord=experiment_number)
dnpData.attrs["experiment_type"] = "nmr_experiment"

# %%
# Apodize and fourier trasnform the data
dnpData = dnp.apodize(dnpData, lw=50)
dnpData = dnp.fourier_transform(dnpData)

# %%
# Looking at the data we see that the data is not correctly phased.

dnp.fancy_plot(dnpData, title="Raw example data")
dnp.plt.gca().set_xlabel("Offset (a.u.)", fontsize=22)
dnp.plt.gca().set_ylabel("Real part of spectrum", fontsize=22)
dnp.show()

# %%
# Let's autophase the data, each experiment indivdually

from dnplab.processing.phase import autophase

dnpAutophase = autophase(dnpData, dim="f2", deriv=1)

# %%
# then we plot each autophased spectra for the individually phased spectra

dnp.fancy_plot(dnpAutophase, title="Individually autophased data")
dnp.plt.gca().set_xlabel("Offset (a.u.)", fontsize=22)
dnp.plt.gca().set_ylabel("Real part of spectrum", fontsize=22)

dnp.show()

# %%
# now with the last experiment (3) as reference slice element.

dnpAutophase_ref = autophase(
    dnpData, dim="f2", reference_slice=("experiment_number", 2), deriv=1
)
dnp.fancy_plot(dnpAutophase_ref, title="Individually autophased data")
dnp.plt.gca().set_xlabel("Offset (a.u.)", fontsize=22)
dnp.plt.gca().set_ylabel("Real part of spectrum", fontsize=22)
dnp.show()

# %%
# From the processing attributes it is clear that the adjusted phase is different between the two processed data!
dnpAutophase.proc_info()
dnpAutophase_pivot.proc_info()
