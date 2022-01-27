import numpy as np

from ..core.data import DNPData

def calculate_enhancement(
    all_data,
    off_spectrum=1,
    on_spectra="all",
    integrate_center=0,
    integrate_width="full",
    method="integrate",
    dim="f2",
):
    """Calculate enhancement from DNP data

    Args:
        all_data (dnpdata, dict): data container

    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | parameter        | type                       | default     | description                                                          |
    +==================+============================+=============+======================================================================+
    | off_spectrum     | int or dnpdata             | 1           | slice of 2D data to be used as p = 0 spectrum, or dnpdata            |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | on_spectra       | str or dnpdata             | "all"       | "all"  unless dnpdata given                                          |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_center | str, int, list, or ndarray | 0           | "max", center of integration window, or index used to find amplitude |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_width  | str, int, list, or ndarray | "full"      | "full" or width of integration window                                |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | method           | str                        | "integrate" | either "integrate" or "ampltiude"                                    |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | dim              | str                        | "f2"        | dimension to integrate down or search down for max                   |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+

    Returns:
        dnpdata: data object with "enhancements" key added

    """
    if off_spectrum == 0:
        raise ValueError(
            "please use indices from 1 to number of slices, i.e. use 1 instead of 0"
        )

    _, isDict = return_data(all_data)

    if isDict and "integrals" in all_data.keys() and isinstance(off_spectrum, int):
        enh = np.array(
            all_data["integrals"].values.real
            / all_data["integrals"].values.real[off_spectrum - 1]
        )
        enhancement_data = dnpdata(
            enh,
            [all_data["integrals"].coords[x] for x in all_data["integrals"].dims],
            all_data["integrals"].dims,
        )

    elif isinstance(off_spectrum, dnpdata) and isinstance(on_spectra, dnpdata):

        data_off, _ = return_data(off_spectrum)
        if len(data_off.shape) != 1:
            raise TypeError("off_spectrum should be 1D")
        data_on, _ = return_data(on_spectra)
        index_on = data_on.dims.index(dim)

        if integrate_width == "full":
            int_width_off = max(data_off.coords[dim]) - min(data_off.coords[dim])
            int_width_on = max(data_on.coords[dim]) - min(data_on.coords[dim])
        elif (
            isinstance(integrate_width, (list, np.ndarray))
            and len(integrate_width) == 2
        ):
            int_width_off = integrate_width[0]
            int_width_on = integrate_width[1]
        elif isinstance(integrate_width, int):
            int_width_off = integrate_width
            int_width_on = integrate_width
        else:
            raise ValueError(
                "integrate_width must be an integer, a list/array with len = 2, or 'full'"
            )

        if integrate_center == "max":
            on_maxs = np.argmax(data_on.values.real, axis=index_on)
            int_center_off = data_off.coords[dim][np.argmax(data_off.values.real)]
            if on_maxs.size == 1:
                int_center_on = data_on.coords[dim][on_maxs]
            else:
                int_center_on = [data_on.coords[dim][x] for x in on_maxs]
        elif (
            isinstance(integrate_center, (list, np.ndarray))
            and len(integrate_center) == 2
        ):
            int_center_off = integrate_center[0]
            int_center_on = integrate_center[1]
        elif isinstance(integrate_center, int):
            int_center_off = integrate_center
            int_center_on = integrate_center
        else:
            raise ValueError(
                "integrate_center must be an integer, a list/array with len = 2, or 'max'"
            )

        if method == "integrate":
            off_data = dnpTools.integrate(
                data_off,
                dim=dim,
                integrate_center=int_center_off,
                integrate_width=int_width_off,
            )

            on_data = dnpTools.integrate(
                data_on,
                dim=dim,
                integrate_center=int_center_on,
                integrate_width=int_width_on,
            )

            enh = np.array(on_data.values.real / off_data.values.real)
            enhancement_data = DNPData(
                enh, [on_data.coords[x] for x in on_data.dims], on_data.dims
            )

        elif method == "amplitude":
            on_maxs = np.argmax(abs(data_on.values.real), axis=index_on)
            if integrate_center == "max":
                off_data = data_off.values.real[np.argmax(abs(data_off.values.real))]
                if on_maxs.size == 1:
                    on_data = data_on.values.real[on_maxs]
                else:
                    on_data = [
                        data_on.values.real[x, indx] for indx, x in enumerate(on_maxs)
                    ]
            else:
                off_data = data_off.values.real[int_center_off]
                if on_maxs.size == 1:
                    on_data = data_on.values.real[int_center_on]
                else:
                    on_data = [
                        data_on.values.real[int_center_on, indx]
                        for indx, _ in enumerate(on_maxs)
                    ]

            if (isinstance(on_data, list) and len(on_data) == 1) or (
                isinstance(on_data, float) and on_data.size == 1
            ):
                enh = np.array([on_data / off_data])
            else:
                enh = np.array(on_data / off_data)

            remaining_dims = [x for x in data_on.dims if x != dim]
            if len(remaining_dims) == 0:
                remaining_dims = ["index"]
                remaining_coords = [np.array([0])]
            else:
                remaining_coords = [data_on.coords[x] for x in data_on.dims if x != dim]

            enhancement_data = dnpdata(enh, remaining_coords, remaining_dims)

    else:
        raise TypeError(
            "Either use the integrate function first and define the index of the off spectrum, or pass dnpata objects for off_spectrum and on_spectra"
        )

    return enhancement_data

def signal_to_noise(
    data,
    dim="f2",
    signal_center=0,
    signal_width="full",
    noise_center="default",
    noise_width="default",
):
    """Find signal-to-noise ratio

    .. note::

        S/N = signal / stdd(noise)

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+-------+-----------+------------------------------+
    | parameter        | type  | default   | description                  |
    +==================+=======+===========+==============================+
    | dim              | str   | 'f2'      | dimension                    |
    +------------------+-------+-----------+------------------------------+
    | signal_center    | float | 0         | center of signal             |
    +------------------+-------+-----------+------------------------------+
    | signal_width     | float | "full"    | width of signal              |
    +------------------+-------+-----------+------------------------------+
    | noise_center     | float | "default" | center of noise region       |
    +------------------+-------+-----------+------------------------------+
    | noise_width      | float | "default" | width of noise region        |
    +------------------+-------+-----------+------------------------------+

    Returns:
        dnpdata: data object with attrs "s_n", "signal", and "noise" added
    """

    index = data.dims.index(dim)

    if signal_width == "full" and isinstance(signal_center, (int, float)):
        s_data = data[dim, :].real
    elif isinstance(signal_width, (int, float)) and isinstance(
        signal_center, (int, float)
    ):
        signalMin = signal_center - np.abs(signal_width) / 2.0
        signalMax = signal_center + np.abs(signal_width) / 2.0
        s_data = data[dim, (signalMin, signalMax)].real
    else:
        raise ValueError(
            "signal_center and signal_width must be int or float, signal_width may also be 'full'"
        )

    if noise_center == "default" and noise_width == "default":
        noise_width = 0.05 * (max(data.coords[dim]) - min(data.coords[dim]))
        noise_center = max(data.coords[dim]) - (np.abs(noise_width) / 2.0)
    elif isinstance(noise_center, (int, float)) and noise_width == "default":
        noise_width = 0.05 * (max(data.coords[dim]) - min(data.coords[dim]))
        if noise_center + (np.abs(noise_width) / 2.0) > max(data.coords[dim]):
            noise_width = 2 * (max(data.coords[dim]) - noise_center)
    elif isinstance(noise_width, (int, float)) and noise_center == "default":
        noise_center = max(data.coords[dim]) - (np.abs(noise_width) / 2.0)
    elif isinstance(noise_width, (int, float)) and isinstance(
        noise_center, (int, float)
    ):
        pass
    else:
        raise ValueError(
            "noise_center and noise_width must be int, float, or 'default'"
        )

    noiseMin = noise_center - np.abs(noise_width) / 2.0
    noiseMax = noise_center + np.abs(noise_width) / 2.0
    n_data = data[dim, (noiseMin, noiseMax)].real

    if len(data.dims) == 1:
        s_n = s_data.values[np.argmax(s_data.values)] / np.std(n_data.values)
    else:
        sn_maxs = np.argmax(s_data.values, axis=index)
        s_n = [
            s_data.values[x, ix] / np.std(n_data.values[:, ix], axis=index)
            for ix, x in enumerate(sn_maxs)
        ]

    data.attrs["s_n"] = s_n
    data.attrs["signal"] = s_data
    data.attrs["noise"] = n_data

    return data

def left_shift(data, dim="t2", shift_points=0):
    """Remove points from the left of data

    Args:
        all_data (dnpdata, dict): Data container for data

    +---------------+------+---------+--------------------------------------------------+
    | parameter     | type | default | description                                      |
    +===============+======+=========+==================================================+
    | dim           | str  | "t2"    | dimension to shift                               |
    +---------------+------+---------+--------------------------------------------------+
    | shift_points  | int  | 0       | Number of points to remove from left of data     |
    +---------------+------+---------+--------------------------------------------------+

    Returns:
        dnpdata: data object with left-shifted data
    """

    data = data[dim, shift_points:]

    proc_attr_name = "left_shift"
    proc_parameters = {
        "dim": dim,
        "points": shift_points,
    }
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data

def normalize():
    return NotImplemented

def reference():
    return NotImplemented