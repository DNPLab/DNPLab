from warnings import warn

import re

import numpy as np

__all__ = ["fourier_transform", "inverse_fourier_transform"]


def convert_to_ppm(freq_coord):
    return NotImplemented


def rename_ft_dim(dim, old_string, new_string):
    if re.fullmatch("%s[0-9]*" % old_string, dim) is not None:
        dim = dim.replace(old_string, new_string)

    return dim


def fourier_transform(
    data,
    dim="t2",
    zero_fill_factor=1,
    shift=True,
    convert_to_ppm=True,
):
    """Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    .. math::


    Args:
        all_data (dnpdata, dict): Data container

    +------------------+------+-----------+--------------------------------------------------+
    | parameter        | type | default   | description                                      |
    +==================+======+===========+==================================================+
    | dim              | str  | 't2'      | dimension to Fourier transform                   |
    +------------------+------+-----------+--------------------------------------------------+
    | zero_fill_factor | int  | 2         | factor to increase dim with zeros                |
    +------------------+------+-----------+--------------------------------------------------+
    | shift            | bool | True      | Perform fftshift to set zero frequency to center |
    +------------------+------+-----------+--------------------------------------------------+
    | convert_to_ppm   | bool | True      | Convert dim from Hz to ppm                       |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after FT
    """

    # handle zero_fill_factor
    zero_fill_factor = int(zero_fill_factor)
    if zero_fill_factor <= 0:
        zero_fill_factor = 1

    proc_parameters = {
        "dim": dim,
        "zero_fill_factor": zero_fill_factor,
        "shift": shift,
        "convert_to_ppm": convert_to_ppm,
    }

    index = data.dims.index(dim)

    dt = data.coords[dim][1] - data.coords[dim][0]
    n_pts = zero_fill_factor * len(data.coords[dim])
    f = (1.0 / (n_pts * dt)) * np.r_[0:n_pts]
    if shift == True:
        f -= 1.0 / (2 * dt)

    if convert_to_ppm:
        if "nmr_frequency" not in data.attrs.keys():
            warn(
                "NMR frequency not found in the attrs dictionary, coversion to ppm requires the NMR frequency. See docs."
            )
        else:
            nmr_frequency = data.attrs["nmr_frequency"]
            f /= nmr_frequency / 1.0e6  # updated

    data.values = np.fft.fft(data.values, n=n_pts, axis=index)

    if shift:
        data.values = np.fft.fftshift(data.values, axes=index)

    data.coords[dim] = f

    new_dim = rename_ft_dim(dim, "t", "f")
    data.rename(dim, new_dim)

    proc_attr_name = "fourier_transform"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def inverse_fourier_transform(
    data,
    dim="f2",
    zero_fill_factor=1,
    shift=True,
    convert_from_ppm=True,
):
    """Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    Args:
        all_data (dnpdata, dict): Data container

    +------------------+------+-----------+--------------------------------------------------+
    | parameter        | type | default   | description                                      |
    +==================+======+===========+==================================================+
    | dim              | str  | 't2'      | dimension to Fourier transform                   |
    +------------------+------+-----------+--------------------------------------------------+
    | zero_fill_factor | int  | 2         | factor to increase dim with zeros                |
    +------------------+------+-----------+--------------------------------------------------+
    | shift            | bool | True      | Perform fftshift to set zero frequency to center |
    +------------------+------+-----------+--------------------------------------------------+
    | convert_from_ppm | bool | True      | Convert dim from Hz to ppm                       |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after inverse FT
    """

    # handle zero_fill_factor
    zero_fill_factor = int(zero_fill_factor)
    if zero_fill_factor <= 0:
        zero_fill_factor = 1

    proc_parameters = {
        "dim": dim,
        "zero_fill_factor": zero_fill_factor,
        "shift": shift,
        "convert_from_ppm": convert_from_ppm,
    }

    index = data.dims.index(dim)

    df = data.coords[dim][1] - data.coords[dim][0]
    if convert_from_ppm:
        if "nmr_frequency" not in data.attrs.keys():
            warn(
                "NMR frequency not found in the attrs dictionary, coversion from ppm requires the NMR frequency. See docs."
            )
        else:
            nmr_frequency = data.attrs["nmr_frequency"]
            df /= 1 / (nmr_frequency / 1.0e6)  # updated

    n_pts = zero_fill_factor * len(data.coords[dim])
    t = (1.0 / (n_pts * df)) * np.r_[0:n_pts]

    if shift:
        data.values = np.fft.fftshift(data.values, axes=index)

    data.values = np.fft.ifft(data.values, n=n_pts, axis=index)
    data.coords[dim] = t

    new_dim = rename_ft_dim(dim, "f", "t")
    data.rename(dim, new_dim)

    proc_attr_name = "inverse_fourier_transform"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def zero_fill():
    return NotImplemented
