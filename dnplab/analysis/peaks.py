import numpy as _np
import scipy.signal as _spsig
import dnplab as _dnp
import warnings


def find_peaks(
    data,
    dims="f2",
    normalize=True,
    threshold=0.05,
    height=0.5,
    regions=None,
    peak_info=False
):
    """Find peaks in spectrum

    Find peaks in spectrum (dnpdata object) and returns index, width, and relative peak height. The function uses the SciPy functions "find_peaks" and "peak_widths". Currently this function has only been tested for 1D spectra.

    Args:
        data (DNPData):         Data object
        dims (str):             Dimension to find peaks
        normalize (boolean):    Normalize data to a maximum value of 1. Default is True
        threshold (float):      Threshold of minimum peak height to be counted. Default is 0.05. This value is impacted by the normalize argument
        height (float):         Relative height at which peak width is measured. Default is 0.5 for FWHH
        regions (None, list):   List of tuples defining the region to find peaks (not implemented yet)
        peak_info (boolean):    If True print output to terminal

    Returns:
        data (DNPData):         nd array of peak index, peak width and relative peak height. The linewidth is returned in (Hz), based on the spectrometer frequency

    Examples:
        Find peaks in entire data region:

            >>> peak_list = dnp.find_peaks(data)

        Find peaks with an amplitude > 0.01 (after normalization):

            >>> peak_list = dnp.find_peaks(data, peak_height = 0.05)

        Find peaks with an amplitude > 500 (data not normalized):

            >>> peak_list = dnp.find_peaks(data, peak_height = 500, normalize = False)

    """

    out = data.copy()
    out.attrs["experiment_type"] = "peak_list"
    out.attrs["data_type"] = "peak_list"

    coords = []

    resolution = _np.sum(_np.diff(out.coords)) / _np.size(out.coords)
    frequency = out.attrs["nmr_frequency"]

    if normalize == True:
        out = _dnp.normalize(out)

    if regions == None:
        peak_index, _ = _spsig.find_peaks(out.values.real, height=threshold)
        peak_width_height = _spsig.peak_widths(
            out.values.real, peaks=peak_index, rel_height=height
        )

        peak_width = peak_width_height[0] * resolution * 1e-6 * frequency
        peak_height = peak_width_height[1]

        out.values = _np.vstack((_np.vstack((peak_index, peak_width)), peak_height))

    out = _dnp.update_axis(out, new_dims="index", start_stop=(0, len(out.values)))

    # else:
    #     data_list = []
    #     Not yet implemented for dimensions > 1

    proc_attr_name = "peak_list"
    proc_parameters = {
        "dims": dims,
        "regions": regions,
        "normalize": normalize,
        "threshold": threshold,
        "height": height,
    }

    out.add_proc_attrs(proc_attr_name, proc_parameters)

    if peak_info == True:

        _dnp.peak_info(out)

    return out


def peak_info(data):
    """
    Print peak list in human readable form

    Function to print the peak list in a human readable form. You first have to run find_peaks to create a dnpdata object that includes a peak list.

    Args:
        data (DNPData):     DNPData object created by find_peaks

    Returns:
        Output (str):       Peak list table
    """

    if data.attrs["experiment_type"] != "peak_list":
        print("Peak list required as input")

        return

    array_size = _np.shape(data.values)

    print("----------")
    print("Peak Table")
    print("----------")

    k = 1

    while k < array_size[1] + 1:
        print(
            "Peak %3d: Index: %5d, Width (Hz): %4.2f, Height (rel.): %2.2f"
            % (k, data.values[0][k - 1], data.values[1][k - 1], data.values[2][k - 1])
        )
        k += 1
