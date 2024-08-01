import numpy as _np
from scipy.signal import savgol_filter

from ..core.data import DNPData
from ..processing.integration import integrate
from ..processing.offset import remove_background as dnp_remove_background
from ..constants import constants as _const

# import dnplab as dnp

from scipy.special import jv as _jv
from scipy.fft import fft as _fft
from scipy.fft import ifft as _ifft
import copy as _copy


def calculate_enhancement(data, off_spectrum_index=0, return_complex_values=False):
    """Calculate enhancement of a power series. Needs integrals as input

    Args:
        integrals (DNPData):
        off_spectrum_index (int):
        return_complex_values (bool):

    Returns:
        enhancements (DNPData): Enhancement values
    """

    enhancements = data.copy()

    proc_parameters = {
        "off_spectrum_index": off_spectrum_index,
        "return_complex_values": return_complex_values,
    }

    if not "experiment_type" in data.attrs.keys():
        raise KeyError("Experiment type not defined")

    if data.attrs["experiment_type"] != "integrals":
        raise ValueError("dnpdata object does not contain integrals.")

    if data.dims[0] == "Power":
        enhancements.attrs["experiment_type"] = "enhancements_P"

        enhancements.values = (
            enhancements.values / enhancements.values[off_spectrum_index]
        )

    elif data.dims[0] == "B0":
        enhancements.attrs["experiment_type"] = "enhancements_B0"
        print("This is a DNP enhancement profile. Not implemented yet.")

    else:
        raise TypeError(
            "Integration axis not recognized. First dimension should be Power or B0."
        )

    proc_attr_name = "calculate_enhancement"
    enhancements.add_proc_attrs(proc_attr_name, proc_parameters)

    if return_complex_values == True:
        return enhancements

    elif return_complex_values == False:
        return enhancements.real


def _create_complexEXT(data, real, imag):
    # helper function for create_complex
    complexData = _np.vectorize(complex)(real, imag)

    dims = data.dims
    dims.pop(-1)

    coords = data.coords
    coords = list(coords)
    coords.pop(-1)

    attrs = data.attrs

    out = DNPData(complexData, dims, coords, attrs)

    return out


def _create_complexINT(dnpdata, dim, real=0, imag=1):
    try:
        if len(dnpdata.coords[dim]) != 2:
            raise ValueError(
                "create_complex: Dimension {0} has length {1} != 2".format(
                    dim, len(dnpdata.coords[dim])
                )
            )
    except KeyError:
        raise KeyError(
            "dim {} not in dims of dnpDataobject, available dims are: {1}".format(
                dim, dnpdata.dims
            )
        )
    out = dnpdata.copy()
    out._values = out._values.astype(complex)
    cut_position = 0
    for k in out.dims:
        if k == dim:
            break
        cut_position = cut_position + 1
    if real == 0 and imag == 1:
        out[dim, 0] = out[dim, 0]._values + 1j * out[dim, 1]._values
    elif real == 1 and imag == 0:
        out[dim, 0] = out[dim, 1]._values + 1j * out[dim, 0]._values
    else:
        raise ValueError(
            "create_complex: only real=0 and imag=1 or other way around allowed! you chose {}/{}".format(
                real, imag
            )
        )
    axis_int = 0
    for k in dnpdata.dims:
        if k == dim:
            break
        axis_int = axis_int + 1
    out._values = _np.delete(out._values, 1, axis=axis_int)
    out.coords[dim] = _np.array([0])

    shape = out.shape
    restore_dims = [
        out.dims[x] for x in range(len(shape)) if (shape[x] == 1) and out.dims[x] != dim
    ]
    restore_dims_index = [
        x for x in range(len(shape)) if (shape[x] == 1) and out.dims[x] != dim
    ]
    restore_coords = [out.coords[k][0] for k in restore_dims]  # have length 1
    dims_position = []
    restored_dim_order = []
    i = 0
    for ind, k in enumerate(out.dims):
        if k in restore_dims and ind < cut_position:
            dims_position.append(ind)
        if k in restore_dims and ind > cut_position:
            dims_position.append(ind - 1)

        # save indizes for old dimensions
        if ind in restore_dims_index:
            restored_dim_order.append(
                len(shape) - len(restore_dims) - 1 + i
            )  # -1 from the removed cmplx axis and -1 because counting starts from 0
            i = i + 1
        elif ind < cut_position:
            restored_dim_order.append(ind)
        elif ind > cut_position:
            restored_dim_order.append(ind - 1 - i)

    # remove the single dimensions, including the complex source
    out.squeeze()

    # move single dimensions back to position
    move_dims = 0
    for d, c in zip(restore_dims, restore_coords):
        out.new_dim(d, c)
        move_dims = move_dims + 1
    # move back to old position...
    out._values = _np.moveaxis(
        out._values, [k for k in range(-move_dims, 0, 1)], dims_position
    )
    out.coords.reorder_index(restored_dim_order)

    return out


def create_complex(data, real, imag=None, real_index=0, imag_index=1):
    """Create complex array from input

    This function can be used to concatenate a two dimensions of a DNPData object into a complex array. The unused dims and coords will be removed from the input DNPData object.
    When a String is provided as the second argument the index in that dimension given by real_index is assumed to be the real part of the dataset and the one by imag_index is the iamginary part.
    The dataset is then combined to form one complex dataset, imag is ignored. Note that dimension with size 1 are retained but will be placed at the end of the retuned DNPData object.

    Args:
        data (DNPData): DNPData input object
        real (array, String): Real data if array or when a String is provided the dimension that contains real and imaginary part (the dimension must have length 2)
        imag (array, None): Imaginary data or None, if None is provided a complex dataset is created with the imaginary part set to 0
        real_index (Integer): Index of real part in chosen dimension, default=0, must be 0 or 1 and be different from imag_index
        imag_index (Integer): Index of imaginary part in chosen dimension, default=1, must be 0 or 1 and be different from real_index

    Returns:
        data (DNPData): New DNPData object

    Examples:
    In this example, first a data set is loaded. The data set is of the size 4000 x 2 (ndarray, float32) and the dims are called 't2','x'

    With the first dimension ([...,0]) being the real data and the second ([...,1]) the imaginary data. Using the function create_complex the dnpdata object is converted into a complex data set.

        .. code-block:: python

            data = dnp.load("MyFile.exp")       # Load example data

            data_complex = dnp.create_complex(data, data.values[..., 0], data.values[..., 1])

    Or with the second variant;

        .. code-block:: python

            data = dnp.load("MyFile.exp")       # Load example data

            data_complex = dnp.create_complex(data,'x')


    """
    if (
        type(real) == __builtins__["str"]
    ):  # could be done with try .. except but for now force usage of strings
        dim = real
        return _create_complexINT(data, dim, real=real_index, imag=imag_index)
    if imag is None:
        imag = _np.zeros(real.shape)
    return _create_complexEXT(data, real, imag)


def signal_to_noise(
    data: DNPData,
    signal_region: list = slice(0, None),
    noise_region: list = (None, None),
    dim: str = "f2",
    remove_background: list = None,
    complex_noise=False,
    **kwargs
):
    """Find signal-to-noise ratio

    Simplest implementation: select largest value in a signal_region and divide this value by the estimated std. deviation of another noise_region. If the noise_region list contains (None,None) (the default) then all points except the points +10% and -10% around the maximum are used for the noise_region.

    Args:
        data: Spectrum data
        signal_region (list): list with a single tuple (start,stop) of a region where a signal should be searched, default is [slice(0,None)] which is the whole spectrum
        noise_region (list): list with tuples (start,stop) of regions that should be taken as noise, default is (None,None)
        dim (str): dimension of data that is used for snr calculation, default is 'f2'
        remove_background (list): if this is not None (a list of tuples, or a single tuple) this will be forwarded to dnp.remove_background, together with any kwargs
        complex_noise (bool): Flag that indicates whether the noise should be calculated on the real part of the noise or on the complex data (default = False)
        kwargs : parameters for dnp.remove_background

    Returns:
        SNR (DNPData): DNPData object that contains SNR values, the axis dim is replaced by an axis named "signal_region"

    Examples:

        A note for the usage: regions can be provided as (min,max), slices use indices.
        To use the standard values just use

            >>> snr = dnp.signal_to_noise(data)

        If you want to select a region for the noise and the signal:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)])

        With background subtracted:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)],remove_background=[(123.4,213.5)])

        This function allows to use a single tuple instead of a list with a single tuple for signal_region, noise_region and remove_background. This is for convenience, slices are currently only supoprted for signal_region and noise_region.

            >>> snr = dnp.signal_to_noise(data,(-1.23,300.4),noise_region=[(-400,-240.5),(123.4,213.5],remove_background=(123.4,213.5))

    """
    import warnings
    import scipy.optimize as _scipy_optimize

    # convenience for signal and noise region
    def _convenience_tuple_to_list(possible_region: list):
        if possible_region is None:
            return possible_region
        # we assume its iterable
        try:
            l = len(possible_region)
            if l != 2:
                return possible_region  # it seems to be iterable?
        except TypeError:
            return [possible_region]  # return as is in a list, might be slice
        try:
            # check whether we can interpret it as value
            a = int(possible_region[0])
            return [
                (possible_region[0], possible_region[1])
            ]  # make a list that contains a tuple
        except TypeError:
            if type(possible_region) == list:
                return possible_region
            return [possible_region]

    signal_region = _convenience_tuple_to_list(signal_region)
    if len(signal_region) > 1:
        snr = []
        for sr in signal_region:
            snr.append(
                _np.squeeze(
                    signal_to_noise(
                        data, sr, noise_region, dim, remove_background, **kwargs
                    )._values
                )
            )

        # return DNPData object with dims: signal_region and all other dimensions copied from data
        dims = ["signal_region" if x == dim else x for x in data.dims]
        coords_new = [
            _np.arange(len(signal_region)) if x == dim else data.coords[x]
            for x in data.dims
        ]
        data_new = _np.array(snr)
        snrData = DNPData(data_new, dims, coords_new)
        return snrData

    noise_region = _convenience_tuple_to_list(noise_region)
    remove_background = _convenience_tuple_to_list(remove_background)

    if not (dim in data.dims):
        raise ValueError("dim {0} not in data.dims ({1})".format(dim, data.dims))

    # remove background
    if remove_background is not None:
        deg = kwargs.pop("deg", 1)
        data = dnp_remove_background(data, dim, deg, remove_background)

    # unfold and calculate snr for each fold_index
    sdata = data
    sdata.unfold(dim)

    # currently only absolute value comparison
    signal = []
    for indx in range(sdata.shape[1]):
        signal.append(
            _np.max(_np.abs(sdata[dim, signal_region[0], "fold_index", indx]))
        )

    # now calculate noise
    noise = []
    for indx in range(sdata.shape[1]):
        idata = sdata[dim, :, "fold_index", indx]
        if (None, None) in noise_region:
            signal_arg = _np.argmax(idata[dim, signal_region[0], "fi", 0])
            datasize = idata[dim, :, "fi", 0].size
            noise_region = [
                slice(0, int(_np.maximum(2, int(signal_arg * 0.9)))),
                slice(int(_np.minimum(datasize - 2, int(signal_arg * 1.1))), None),
            ]

        # concatenate noise_regions
        noise_0 = idata[dim, noise_region[0], "fi", 0]
        for k in noise_region[1:]:
            noise_0.concatenate(idata[dim, k, "fi", 0], dim)
        if complex_noise:
            noise.append(_np.std(noise_0[dim, slice(0, None)]))
        else:
            noise.append(_np.std(_np.real(noise_0[dim, slice(0, None)])))

    sdata.fold()

    # return DNPData object
    dims = ["signal_region" if x == dim else x for x in sdata.dims]
    coords_new = [
        _np.arange(1) if x == dim else sdata.coords[x] for x in sdata.dims
    ]  # we know that only one sr is there
    data_new = (_np.array(signal) / _np.array(noise)).reshape(
        [x.size for x in coords_new]
    )
    snrData = DNPData(data_new, dims, coords_new)
    return snrData


def smooth(data, dim="t2", window_length=11, polyorder=3):
    """Apply Savitzky-Golay Smoothing

    This function is a wrapper function for the savgol_filter from the SciPy python package (https://scipy.org/). For a more detailed description see the SciPy help for this function.

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform smoothing
        window_length (int): Length of window (number of coefficients)
        polyorder (int): Polynomial order to fit samples

    Returns:
        data (DNPData): Data with Savitzky-Golay smoothing applied
    """
    out = data.copy()

    proc_parameters = {
        "dim": dim,
        "window_length": window_length,
        "polyorder": polyorder,
    }

    out.unfold(dim)

    out.values = savgol_filter(out.values, window_length, polyorder, axis=0)

    out.fold()

    proc_attr_name = "smooth"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def pseudo_modulation(data, modulation_amplitude, dim="B0", order=1, zero_padding=2):
    """Calculate the first derivative of an EPR spectrum due to field modulation

    Calculation is based on: Hyde et al., “Pseudo Field Modulation in EPR Spectroscopy.”,
    Applied Magnetic Resonance 1 (1990): 483–96.

    Args:
        data (DNPData): DNPData object (typically an absorption line EPR spectrum)
        modulation_amplitude: Peak to peak modulation amplitude. The unit is equal to the unit of the axis. E.g. if the spectrum axis is given in (T), the unit of the modulation amplitude is in (T) as well.
        dim: Dimension to pseudo modulate (default is B0)
        order: Harmonic of field modulation (default is 1, 1st derivative)
        zero_padding: Number of points for zero-padding (multiples of spectrum vector length). Default is 2. Increase this number for short signal vectors.

    Returns:
        data (DNPData): Pseudo modulated spectrum


    Examples:
        .. code-block:: python

            # Calculate pseudo_modulated spectrum (1st derivative). Field axis given in (T)
            spec_mod = dnp.pseudo_modulation(spec, modulation_amplitude = 0.001)

            # Calculate pseudo_modulated spectrum (2nd derivative). Field axis given in (T)
            spec_mod = dnp.pseudo_modulation(spec, modulation_amplitude = 0.001, order = 2)

    """

    out = data.copy()
    out.unfold(dim)

    proc_parameters = {
        "dim": dim,
        "modulation_amplitude": modulation_amplitude,
        "order": order,
        "zero_padding": zero_padding,
    }

    n = len(out.coords[dim])
    delta_B = out.coords[dim][2] - out.coords[dim][1]
    Zmin = 0
    Zmax = _const.pi * modulation_amplitude / delta_B
    Z = _np.linspace(Zmin, Zmax, zero_padding * n)

    spec = out.values
    spec = _np.squeeze(spec)

    fft_spec = _fft(spec, zero_padding * n)  # Zero pad data
    fft_spec[int(n) + 1 : zero_padding * n] = 0

    fft_spec_mod = fft_spec * _jv(order, Z)
    # Convolute fft spectrum with bessel function

    spec_mod = _ifft(fft_spec_mod)
    spec_mod = 1j**order * spec_mod[0:n]  # Pick the right dimension for higher orders
    spec_mod = _np.real(spec_mod)  # Only return real part

    out.values = spec_mod

    out.fold()

    proc_attr_name = "pseudo_modulation"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def left_shift(data, dim="t2", shift_points=0):
    """Remove points from the left

    Args:
        data (DNPData): Data object
        dim (str): Name of dimension to left shift, default is "t2"
        shift_points (int): Number of points to left shift, default is 0.

    Returns:
        data (DNPDdata): Shifted data object
    """

    out = data.copy()

    out = out[dim, shift_points:]

    proc_attr_name = "left_shift"
    proc_parameters = {
        "dim": dim,
        "points": shift_points,
    }
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def normalize(data, amplitude=True, dim="f2", regions=None):
    """Normalize spectrum

    The function is used to normalize the amplitude (or area) of a spectrum to a value of 1. The sign of the original data will be conserved.

    Args:
        data (DNPData):         Data object
        amplitude (boolean):    True: normalize amplitude, false: normalize area. The default is True
        dim (str):              The dimension to normalize
        regions (None, list):   List of tuples to specify range of normalize [(-99., 99.)]

    Returns:
        data (DNPDdata):        Normalized data object
    """

    out = data.copy()

    if amplitude == True:
        if regions:
            factor = _np.max(abs(out["f2", regions].values))
        else:
            factor = _np.max(abs(out.values))

        out.values = out.values / factor
    elif amplitude == False:
        out.values = out.values  # Normalize to area = 1, not implemented yet

    proc_attr_name = "normalized"
    proc_parameters = {
        "amplitude": amplitude,
    }

    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def reference(data, dim="f2", old_ref=0, new_ref=0):
    """Function for referencing NMR spectra

    Args:
        data (DNPData): Data for referencing
        dim (str): dimension to perform referencing down. By default this dimension is "f2".
        old_ref (float): Value of old reference
        new_ref (float): New reference value

    Returns:
        DNPData: referenced data
    """

    out = data.copy()

    out.coords[dim] -= old_ref - new_ref

    proc_attr_name = "reference"
    proc_parameters = {
        "dim": dim,
        "old_ref": old_ref,
        "new_ref": new_ref,
    }

    return out
