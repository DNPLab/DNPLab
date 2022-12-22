"""Modules to calculate DNP enhancement profiles"""

import numpy as _np
from scipy import optimize
import dnplab as _dnp
import warnings


def sim_dnp_profile(
    data,
    B0,
    nucleus="1H",
    dnp_process="SE",
    add_details = False,
    remove_background=True,
    normalize=True,
    integrate=True,
):
    """Simulate DNP enhancment profile

    Simulate DNP enhancement profile based on the EPR spectrum. For more details:

    Banerjee, D., D. Shimon, A. Feintuch, S. Vega, and D. Goldfarb. “The Interplay between the Solid Effect and the Cross Effect Mechanisms in Solid State (1)(3)C DNP at 95 GHz Using Trityl Radicals.” Journal of Magnetic Resonance 230 (May 2013): 212–19.
    https://doi.org/10.1016/j.jmr.2013.02.010.

    Args:
        data (DNPdata): EPR spectrum
        B0 (float): Field position for the DNP experiment in (T)
        nucleus (int): Nucleus for DNP experiment
        dnp_process (int): Select DNP mechanism, SE - Solid Effect, CE/TM - Cross-Effect/Thermal Mixing
        add_details (boolean): Add individual spectra to proc_attrs. Default is False
        remove_background (boolean): Remove 0th order background from EPR spectrum. Default is True
        normalize (boolean): Normalize EPR spectrum to maximum amplitude of 1. Default is True
        integrate (boolean): Integrate EPR spectrum. Default is True

    Returns:
        data (DNPdata): Simulated DNP enhancement profile

    .. math::

    """
    out = data.copy()

    # Some error checks:
    if out.attrs.get("experiment_type") == None:
        print("Error: Key experiment_type not present")
        return

    if out.attrs["experiment_type"] != "epr_spectrum":
        print("Error: EPR spectrum required as input.")
        return

    if len(out.dims) > 1:
        print("Error: This function requires a 1D EPR spectrum as input.")
        return

    if remove_background == True:
        out = _dnp.remove_background(out, dim="B0", deg=0)  # Remove background

    if normalize == True:
        out = _dnp.normalize(out)

    if integrate == True:
        out = _dnp.cumulative_integrate(out, dim="B0")  # Calculate cumsum

    ## Calculate number of points to shift
    dnpLarmorFrequency = _dnp.mr_properties(nucleus, B0)  # Nuclear Larmor Frequency
    dnpLarmorFrequency_G = dnpLarmorFrequency / (
        1000 * _dnp.mr_properties("0e") / 2 / _np.pi
    )  # Nuclear Larmor Frequency in [G]
    deltaB0_G = (out.coords["B0"][1] - out.coords["B0"][0]) * 10
    points_to_shift = round(dnpLarmorFrequency_G / deltaB0_G)

    ## Shift EPR spectra using mumpy's roll function
    EPRdataPos = out.copy()
    EPRdataPos.values = _np.roll(EPRdataPos.values, points_to_shift)

    EPRdataNeg = out.copy()
    EPRdataNeg.values = _np.roll(out.values, (-1) * points_to_shift)
    EPRdataNeg.values = (-1) * EPRdataNeg.values

    match dnp_process:
        case "SE":
            out = EPRdataPos + EPRdataNeg

        case "CE/TM":
            out = out.values * (EPRdataPos + EPRdataNeg)

    proc_attr_name = "sim_dnp_profile"
    proc_parameters = {
        "nucleus": nucleus,
        "dnp_process": dnp_process,
        "remove_background": remove_background,
        "normalize": normalize,
        "integrate": integrate,
    }

    out.add_proc_attrs(proc_attr_name, proc_parameters)

    if (add_details == True):
        sim_data = _np.array([out.values, EPRdataNeg.values, EPRdataPos.values])

        proc_attr_name = "sim_dnp_profile"
        proc_parameters = {
            "sim_data": sim_data,
        }
    
        out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


