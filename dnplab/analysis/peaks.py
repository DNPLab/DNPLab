import numpy as _np
import scipy.signal as _spsig
import dnplab as _dnp


def find_peaks(
    data,
    dims="f2",
    normalize=True,
    peak_height=0.05,
    regions=None,
):
    """Find peaks in spectrum

    Find peaks in spectrum (dnpdata object). This is a wrapper for the find_peaks function of SciPy.

    Args:
        data (DNPData): Data object
        dims (str): Dimension in which to find peaks
        normalize (boolean): Normalize data to a maximum value of 1. Default is True
        peak_height (float): Threshold of minimum peak height to be counted. Default is 0.05. This value is impacted by the normalize argument
        regions (None, list): List of tuples defining the region to find peaks

    Returns:
        data (DNPData): Peak list

    Examples:
        Find peaks in entire data region:

            >>> peak_list = dnp.find_peaks(data)

        Find peaks with an amplitude > 0.01 (after normalization):

            >>> peak_list = dnp.find_peaks(data, peak_height = 0.05)

        Find peaks with an amplitude > 500 (data not normalized):

            >>> peak_list = dnp.find_peaks(data, peak_height = 500, normalize = False)

    """
    data = data.copy()
    data.attrs["experiment_type"] = "peak_list"

    coords = []

    if normalize == True:
        data = _dnp.normalize(data)

    if regions == None:
        data.values, _ = _spsig.find_peaks(data.values.real, height=peak_height)
        data.coords.pop(dims)
        data.dims = ["x0"]

        new_coords = _np.arange(0, len(data.values))
        data.coords = new_coords

    # else:
    #     data_list = []
    #     Not yet implemented for dimensions > 1

    proc_attr_name = "peak_list"
    proc_parameters = {
        "dim": dims,
        "regions": regions,
    }

    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data
