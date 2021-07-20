from . import dnpMath, dnpNMR, dnpdata, dnpdata_collection
import numpy as np
import scipy.integrate

from .mrProperties import gmrProperties, radicalProperties


def concat(data_list, dim, coord=None):
    """Concatenates list of data objects down another dimension

    args:
        data_list (list): List of dnpdata objects to concatentate
        dim (str): new dimension name
        coord: coords for new dimension

    Returns:
        data (dnpdata): concatenated data object

    """

    shape = data_list[0].shape
    values_list = [data.values for data in data_list]

    for values in values_list:
        this_shape = values.shape
        if this_shape != shape:
            raise IndexError(
                "Cannot concatenate data objects. Array shapes do not match.",
                this_shape,
                shape,
            )

    dims = data_list[0].dims
    coords = data_list[0].coords.coords
    attrs = data_list[0].attrs

    values = np.stack(values_list, axis=-1)

    dims.append(dim)

    if coord is None:
        coords.append(values_list)
    else:
        coords.append(coord)

    data = dnpdata(values, coords, dims, attrs)

    return data


def return_data(all_data):

    is_workspace = False
    if isinstance(all_data, dnpdata):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError("Type dict is not supported")
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError("No data in processing buffer")
    else:
        raise ValueError("Data type not supported")

    return data, is_workspace


def baseline(
    all_data,
    dim="f2",
    indirect_dim=None,
    type="polynomial",
    order=1,
    reference_slice=None,
):
    """Baseline correction of NMR spectra down given dimension

    Args:
        all_data (object) : dnpdata object

    +-----------------+------+---------------+---------------------------------------------------+
    | parameter       | type | default       | description                                       |
    +=================+======+===============+===================================================+
    | dim             | str  | 'f2'          | Dimension to apply baseline correction            |
    +-----------------+------+---------------+---------------------------------------------------+
    | indirect_dim    | str  | None          | indirect dimension                                |
    +-----------------+------+---------------+---------------------------------------------------+
    | type            | str  | 'polynomial'  | type of baseline fit                              |
    +-----------------+------+---------------+---------------------------------------------------+
    | order           | int  | 1             | polynomial order, or 1=mono 2=bi exponential      |
    +-----------------+------+---------------+---------------------------------------------------+
    | reference_slice | int  | None          | slice of 2D data used to define the baseline      |
    +-----------------+------+---------------+---------------------------------------------------+

    Returns:
        dnpdata: Baseline corrected data, with attr "baseline" added
    """

    data, isDict = return_data(all_data)

    if not indirect_dim:
        if len(data.dims) == 2:
            ind_dim = list(set(data.dims) - set([dim]))[0]
        elif len(data.dims) == 1:
            ind_dim = data.dims[0]
        else:
            raise ValueError(
                "you must specify the indirect dimension, use argument indirect_dim= "
            )
    else:
        ind_dim = indirect_dim

    if reference_slice is not None:
        if len(np.shape(data.values)) == 1:
            reference_slice = None
            warnings.warn("ignoring reference_slice, this is 1D data")
        else:
            reference_slice -= 1

    if len(np.shape(data.values)) == 2:
        if reference_slice is not None:
            bline = dnpMath.baseline_fit(
                data.coords[dim], data.values[:, reference_slice], type, order
            )
            for ix in range(len(data.coords[ind_dim])):
                data.values[:, ix] -= bline
        elif reference_slice is None:
            for ix in range(len(data.coords[ind_dim])):
                bline = dnpMath.baseline_fit(
                    data.coords[dim], data.values[:, ix], type, order
                )
                data.values[:, ix] -= bline
        else:
            raise TypeError("invalid reference_slice")

    elif len(np.shape(data.values)) == 1:
        bline = dnpMath.baseline_fit(data.coords[dim], data.values, type, order)
        data.values -= bline

    else:
        raise ValueError("1D or 2D only")

    proc_parameters = {
        "dim": dim,
        "type": type,
        "order": order,
        "reference_slice": reference_slice,
    }
    proc_attr_name = "baseline"
    data.add_proc_attrs(proc_attr_name, proc_parameters)
    data.attrs["baseline"] = bline

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def integrate(
    all_data, dim="f2", type="single", integrate_center=0, integrate_width="full"
):
    """Integrate data down given dimension

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+---------------+----------+-------------------------------+
    | parameter        | type          | default  | description                   |
    +==================+===============+==========+===============================+
    | dim              | str           | 'f2'     | dimension to integrate        |
    +------------------+---------------+----------+-------------------------------+
    | type             | str           | 'single' | 'single' or 'double' integral |
    +------------------+---------------+----------+-------------------------------+
    | integrate_center | float or list | 0        | center of integration window  |
    +------------------+---------------+----------+-------------------------------+
    | integrate_width  | float or list | "full"   | width of integration window   |
    +------------------+---------------+----------+-------------------------------+

    Returns:
        dnpdata: integrals of data
    """

    data, isDict = return_data(all_data)
    index = data.index(dim)

    if len(data.dims) == 2:
        ind_dim = list(set(data.dims) - set([dim]))[0]
        indirect_coords = data.coords[ind_dim]
    elif len(data.dims) == 1:
        ind_dim = "index"
        indirect_coords = [0]

    data_new = None
    if type == "double":
        first_int = scipy.integrate.cumtrapz(
            data.values, x=data.coords[dim], axis=index, initial=0
        )
        data.values = first_int

    if integrate_width == "full":
        pass
    elif (isinstance(integrate_width, int) or isinstance(integrate_width, float)) and (
        isinstance(integrate_center, int) or isinstance(integrate_center, float)
    ):
        integrateMin = integrate_center - np.abs(integrate_width) / 2.0
        integrateMax = integrate_center + np.abs(integrate_width) / 2.0
        data = data[dim, (integrateMin, integrateMax)]

    elif (
        (isinstance(integrate_width, list) and isinstance(integrate_center, list))
        and (all((isinstance(x, int) or isinstance(x, float)) for x in integrate_width))
        and (
            all((isinstance(x, int) or isinstance(x, float)) for x in integrate_center)
        )
    ):
        if len(integrate_width) != len(integrate_center):
            raise TypeError(
                "If integrate_center and integrate_width are both lists, they must be the same length"
            )

        integrateMin = []
        integrateMax = []
        for x, cent in enumerate(integrate_center):
            integrateMin.append((cent - np.abs(integrate_width[x]) / 2.0))
            integrateMax.append((cent + np.abs(integrate_width[x]) / 2.0))
        data_new = []
        for mx, mn in enumerate(integrateMin):
            data_new.append(data[dim, (mn, integrateMax[mx])])
    elif (
        (isinstance(integrate_width, int) or isinstance(integrate_width, float))
        and isinstance(integrate_center, list)
    ) and (all((isinstance(x, int) or isinstance(x, float)) for x in integrate_center)):
        integrateMin = []
        integrateMax = []
        for cent in integrate_center:
            integrateMin.append((cent - np.abs(integrate_width) / 2.0))
            integrateMax.append((cent + np.abs(integrate_width) / 2.0))
        data_new = []
        for mx, mn in enumerate(integrateMin):
            data_new.append(data[dim, (mn, integrateMax[mx])])

    else:
        raise ValueError(
            "integrate_width must be 'full', int, float, or list of int or float; integrate_center must be int, float, or list of ints or floats"
        )

    if data_new and isinstance(data_new, list):
        data_integrals = []
        for x in data_new:
            data_integrals.append(np.trapz(x.values, x=x.coords[dim], axis=index))

        data.values = np.array(data_integrals)
        int_coords = [integrate_center, indirect_coords]
        indirect_dim = ["center", ind_dim]

    else:
        data.values = np.trapz(data.values, x=data.coords[dim], axis=index)
        int_coords = [indirect_coords]
        indirect_dim = [ind_dim]

    integrate_data = dnpdata(data.values, int_coords, indirect_dim)
    if type == "double":
        integrate_data.attrs["first_integral"] = first_int

    if isDict:
        all_data["integrals"] = integrate_data
    else:
        return integrate_data


def mr_properties(nucleus, *args):
    """Return magnetic resonance property of specified isotope.

    This function is modeled after the Matlab function gmr written by Mirko Hrovat: https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties

    Also see: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818. Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants, http://physics.nist.gov/PhysRefData/codata86/codata86.html, or http://www.isis.rl.ac.uk/neutronSites/constants.htm. Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

    Args:

        nucleus:          '1H', '2H', '6Li', '13C', 14N', etc.
        numerical:        If only a numerical is given in addition to the nucleus it must be a B0 value in Tesla and the Larmor frequency will be returned

    +------------------+----------------------------------------------------------------------------+
    | args             |  returns                                                                   |
    +==================+============================================================================+
    | "gamma"          | Gyromagnetic Ration [radians/T/s]                                          |
    +------------------+----------------------------------------------------------------------------+
    | "spin"           | Spin number of selected nucleus [1]                                        |
    +------------------+----------------------------------------------------------------------------+
    | "qmom"           | Quadrupole moment [fm^2] (100 barns)                                       |
    +------------------+----------------------------------------------------------------------------+
    | "natAbundance"   | Natural abundance [%]                                                      |
    +------------------+----------------------------------------------------------------------------+
    | "relSensitivity" | Relative sensitiviy with respect to 1H at constant B0                      |
    +------------------+----------------------------------------------------------------------------+
    | "moment"         | Magnetic dipole moment, abs(u)/uN = abs(gamma)*hbar[I(I + 1)]^1/2/uN,      |
    +------------------+----------------------------------------------------------------------------+
    | "qlw"            | quadrupolar line-width factor, Qlw = Q^2(2I + 3)/[I^2(2I + 1)]             |
    +------------------+----------------------------------------------------------------------------+


    Examples:
        .. code-block:: python

            dnp.dnpTools.mrProperties('1H') = 26.7522128 # 1H Gyromagnetic Ratio (10^7r/Ts)

            dnp.dnpTools.mrProperties('1H', 0.35) = 14902114.17018196 # 1H Larmor Freq at .35 T (Hz)

            dnp.dnpTools.mrProperties('2H', 'qmom') = 0.286 # Nuclear Quadrupole Moment (fm^2)

            dnp.dnpTools.mrProperties('6Li', 'natAbundance') = 7.59 # % Natural Abundance

            dnp.dnpTools.mrProperties('6Li', 'relSensitivity') = 0.000645 # Relative sensitivity
    """

    if isinstance(nucleus, str):
        if nucleus in gmrProperties:
            gmr = gmrProperties.get(nucleus)[1]
        else:
            print("Isotope doesn't exist in list")
            return
    else:
        print("ERROR: String expected")

    if len(args) == 0:
        return gmr

    elif len(args) == 1:

        if isinstance(args[0], str):
            if args[0] == "gamma":
                return gmrProperties.get(nucleus)[1]

            if args[0] == "spin":
                return gmrProperties.get(nucleus)[0]

            elif args[0] == "qmom":
                return gmrProperties.get(nucleus)[2]

            elif args[0] == "natAbundance":
                return gmrProperties.get(nucleus)[3]

            elif args[0] == "relSensitivity":
                return gmrProperties.get(nucleus)[4]

            elif args[0] == "moment":
                return gmrProperties.get(nucleus)[5]

            elif args[0] == "qlw":
                return gmrProperties.get(nucleus)[6]

            else:
                print("Keyword not recognize")

        else:
            vLarmor = args[0] * gmr * 1e7 / 2 / np.pi
            return vLarmor

    elif len(args) == 2:

        if args[1] == True:
            print(" ")
            print("Nucleus                    : ", nucleus)
            print("Spin                       : ", gmrProperties.get(nucleus)[0])
            print(
                "Gyromagnetic Ratio [kHz/T] : %5.2f"
                % (gmrProperties.get(nucleus)[1] * 10 / 2 / np.pi)
            )
            print(
                "Natural Abundance      [%%] : %5.2f" % (gmrProperties.get(nucleus)[3])
            )
            print("")

    elif len(args) > 2:
        print("Too many input arguments")


def radical_properties(name):
    """Return properties of different radicals. At the minimum the g value is returned. If available, large hyperfine couplings to a nucleus are returned. Add new properties or new radicals to mrProperties.py

    Args:

    +-----------+---------------------------------------------------------------+
    | arg       |  returns                                                      |
    +===========+===============================================================+
    | "gfree"   | 2.00231930436153                                              |
    +-----------+---------------------------------------------------------------+
    | "tempo1"  | [[2.00980, 2.00622, 2.00220], "14N", [16.8, 20.5, 95.9]]      |
    +-----------+---------------------------------------------------------------+
    | "tempo2"  | [[2.00909, 2.00621, 2.00222], "14N", [20.2, 20.2, 102.1]]     |
    +-----------+---------------------------------------------------------------+
    | "bdpa"    | [[2.00263, 2.00260, 2.00257], "1H", [50.2, 34.5, 13.0]]       |
    +-----------+---------------------------------------------------------------+

    Returns:
        principle g values and hyperfine coupling tensor
    """

    name = name.lower()
    if isinstance(name, str):
        if name in radicalProperties:
            giso = radicalProperties.get(name)[0]
        else:
            print("Radical doesn't exist in dictonary")
            return
    else:
        print("ERROR: String expected")

    return giso


def show_dnp_properties(radical, mwFrequency, dnpNucleus):
    """Calculate DNP Properties

    Currently only implemented for liquid state experiments

    Args:
        radical:        Radical name, see mrProperties.py
        mwFrequency:    Microwave frequency in (Hz)
        dnpNucleus:     Nucleus for DNP-NMR experiments

    Example:
        .. code-block:: python

            dnp.dnpTools.show_dnp_poperties('gfree', 9.45e9, '1H')
    """

    # http://physics.nist.gov/constants
    mub = 9.27400968e-24
    planck = 6.62606957e-34

    # Get radical properties
    glist = radicalProperties.get(radical)[0]
    nucleus = radicalProperties.get(radical)[1]
    Alist = radicalProperties.get(radical)[2]

    # Get g-value
    g = np.array(glist)
    giso = np.sum(g) / g.size

    B0 = mwFrequency * planck / giso / mub

    # Get hyperfine coupling and calculate isotropic value
    A = np.array(Alist)
    AisoMHz = np.sum(A) / A.size

    gmr_e = mr_properties("0e")
    AisoT = AisoMHz / gmr_e / 2 / np.pi

    if nucleus != None:
        nucSpin = mr_properties(nucleus, "spin")
        n = 2 * nucSpin + 1
        ms = np.linspace(-1.0 * nucSpin, nucSpin, int(n))
        B = B0 + ms * AisoT

    else:
        nucSpin = 0
        B = B0

    print("")
    print("Input Parameters: ")
    print("Radical                  : ", radical)
    print("giso                     :  %8.6f" % giso)
    print("Nucleus                  : ", nucleus)
    print("Nuc Spin                 : ", nucSpin)
    print("Aiso               (MHz) :  %4.2f" % AisoMHz)
    print("")
    print("Predicted Field Values for DNP: ")
    m = 1
    for b in B:
        print("Transition: ", m)
        print("B                    (T) :  %6.4f" % b)
        nmr = mr_properties("1H") * b * 10 / 2 / np.pi
        print("NMR Frequency      (MHz) :  %6.3f" % nmr)
        print("")
        m += 1


def signal_to_noise(
    all_data,
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

    data, isDict = return_data(all_data)

    if signal_width == "full" and (
        isinstance(signal_center, int) or isinstance(signal_center, float)
    ):
        s_data = data.real
    elif (isinstance(signal_width, int) or isinstance(signal_width, float)) and (
        isinstance(signal_center, int) or isinstance(signal_center, float)
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
    elif (
        isinstance(noise_center, int) or isinstance(noise_center, float)
    ) and noise_width == "default":
        noise_width = 0.05 * (max(data.coords[dim]) - min(data.coords[dim]))
        if noise_center + (np.abs(noise_width) / 2.0) > max(data.coords[dim]):
            noise_width = 2 * (max(data.coords[dim]) - noise_center)
    elif (
        isinstance(noise_width, int) or isinstance(noise_width, float)
    ) and noise_center == "default":
        noise_center = max(data.coords[dim]) - (np.abs(noise_width) / 2.0)
    elif (isinstance(noise_width, int) or isinstance(noise_width, float)) and (
        isinstance(noise_center, int) or isinstance(noise_center, float)
    ):
        pass
    else:
        raise ValueError(
            "noise_center and noise_width must be int, float, or 'default'"
        )

    noiseMin = noise_center - np.abs(noise_width) / 2.0
    noiseMax = noise_center + np.abs(noise_width) / 2.0
    n_data = data[dim, (noiseMin, noiseMax)].real

    if data.ndim == 2:
        sig = []
        noi = []
        for ix in range(data.shape[1]):
            sig.append(s_data.values[np.argmax(s_data.values[:, ix], axis=0), ix])
            noi.append(np.std(n_data.values[:, ix], axis=0))

        s_n = np.array(sig) / np.array(noi)
    elif data.ndim == 1:
        s_n = s_data.values[np.argmax(s_data.values, axis=0)] / np.std(
            n_data.values, axis=0
        )
    else:
        raise TypeError("only 1D or 2D data currently supported")

    data.attrs["s_n"] = s_n
    data.attrs["signal"] = s_data
    data.attrs["noise"] = n_data

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def zero_fill(
    all_data,
    dim="t2",
    zero_fill_factor=2,
    shift=True,
    inverse=False,
    convert_from_ppm=True,
):

    """Perform zero filling down dim dimension

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
    | inverse          | bool | False     | True means zero-fill in frequency domain         |
    +------------------+------+-----------+--------------------------------------------------+
    | convert_from_ppm | bool | True      | True if frequency axis is in ppm                 |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after zero fill
    """

    data, isDict = return_data(all_data)

    # handle zero_fill_factor
    if not isinstance(zero_fill_factor, int) or zero_fill_factor <= 0:
        raise ValueError("zero_fill_factor must be type int greater than 0")

    proc_parameters = {
        "dim": dim,
        "zero_fill_factor": zero_fill_factor,
        "shift": shift,
        "inverse": inverse,
    }

    if inverse:
        df = data.coords[dim][1] - data.coords[dim][0]
        if convert_from_ppm:
            df /= -1 / (data.attrs["nmr_frequency"] / 1.0e6)

        n_pts = zero_fill_factor * len(data.coords[dim])
        data.coords[dim] = (1.0 / (n_pts * df)) * _np.r_[0:n_pts]

        proc_parameters["convert_from_ppm"] = convert_from_ppm

    else:
        dt = data.coords[dim][1] - data.coords[dim][0]
        n_pts = zero_fill_factor * len(data.coords[dim])
        data.coords[dim] = (1.0 / (n_pts * dt)) * np.r_[0:n_pts]
        if shift == True:
            data.coords[dim] -= 1.0 / (2 * dt)

    if len(data.shape) == 2:
        temp = np.zeros((n_pts, data.shape[1]), dtype=np.complex)
        temp[: data.shape[0], : data.shape[1]] = data.values
    elif len(data.shape) == 1:
        temp = np.zeros(n_pts, dtype=np.complex)
        temp[: data.shape[0]] = data.values

    data.values = temp

    proc_attr_name = "zero_fill"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data
