from warnings import warn

from . import return_data, dnpdata, dnpdata_collection
from . import dnpMath, dnpNMR
import numpy as np
import scipy.integrate
from scipy.special import wofz

from .mrProperties import gmrProperties, radicalProperties


def baseline(
    all_data,
    dim="f2",
    type="polynomial",
    order=1,
    p0=None,
    mode="subtract",
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
    | type            | str  | 'polynomial'  | type of baseline fit                              |
    +-----------------+------+---------------+---------------------------------------------------+
    | order           | int  | 1             | polynomial order, or 1=mono 2=bi exponential      |
    +-----------------+------+---------------+---------------------------------------------------+
    | p0              | list | None          | initial guess for exponential baseline fit        |
    +-----------------+------+---------------+---------------------------------------------------+
    | mode            | str  | 'subtract'    | either 'subtract' or 'divide' by baseline         |
    +-----------------+------+---------------+---------------------------------------------------+
    | reference_slice | int  | None          | slice of 2D data used to define the baseline      |
    +-----------------+------+---------------+---------------------------------------------------+

    Returns:
        dnpdata: Baseline corrected data, with attr "baseline" added
    """

    data, isDict = return_data(all_data)
    index = data.dims.index(dim)

    if reference_slice is not None:
        if len(np.shape(data.values)) == 1:
            reference_slice = None
            warn("ignoring reference_slice, this is 1D data")
        else:
            reference_slice -= 1

    if len(data.dims) == 1:
        bline = dnpMath.baseline_fit(data.coords[dim], data.values, type, order, p0=p0)
        if mode == "subtract":
            data.values -= bline
        elif mode == "divide":
            data.values /= bline
    else:
        ind_dim = list(set(data.dims) - set([dim]))[0]
        ind_shape = data.shape[data.index(ind_dim)]
        bline_array = np.zeros(shape=(data.shape[index], ind_shape))
        if reference_slice is not None:
            bline = dnpMath.baseline_fit(
                data.coords[dim],
                data[dim, :].values[:, reference_slice],
                type,
                order,
                p0=p0,
            )
            for ix in range(ind_shape):
                bline_array[:, ix] = bline.real
        elif reference_slice is None:
            for ix in range(ind_shape):
                bline = dnpMath.baseline_fit(
                    data.coords[dim],
                    data[dim, :].values[:, ix],
                    type,
                    order,
                    p0=p0,
                )
                bline_array[:, ix] = bline.real
        else:
            raise TypeError("invalid reference_slice")
        if mode == "subtract":
            data.values -= bline_array
        elif mode == "divide":
            data.values /= bline_array

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
    all_data,
    dim="f2",
    type="trapz",
    integrate_center=0,
    integrate_width="full",
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
    index = data.dims.index(dim)

    data_new = data.copy()
    if type == "double":
        first_int = scipy.integrate.cumtrapz(
            data.values, x=data.coords[dim], axis=index, initial=0
        )
        data_new.values = first_int

    if integrate_width == "full":
        pass
    elif isinstance(integrate_width, (int, float)) and isinstance(
        integrate_center, (int, float)
    ):
        integrateMin = integrate_center - np.abs(integrate_width) / 2.0
        integrateMax = integrate_center + np.abs(integrate_width) / 2.0
        data_new = data_new[dim, (integrateMin, integrateMax)]

    elif (
        isinstance(integrate_width, list)
        and isinstance(integrate_center, list)
        and all((isinstance(x, (int, float)) for x in integrate_width))
        and all((isinstance(x, (int, float)) for x in integrate_center))
    ):
        if len(integrate_width) != len(integrate_center):
            raise TypeError(
                "If integrate_center and integrate_width are both lists, they must be the same length"
            )

        integrateMinMax = [
            [
                (cent - np.abs(integrate_width[x]) / 2.0),
                (cent + np.abs(integrate_width[x]) / 2.0),
            ]
            for x, cent in enumerate(integrate_center)
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]
    elif (
        isinstance(integrate_width, (int, float))
        and isinstance(integrate_center, list)
        and all((isinstance(x, (int, float)) for x in integrate_center))
    ):
        integrateMinMax = [
            [
                (cent - np.abs(integrate_width) / 2.0),
                (cent + np.abs(integrate_width) / 2.0),
            ]
            for cent in integrate_center
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]

    elif (
        isinstance(integrate_center, (int, float))
        and isinstance(integrate_width, list)
        and all((isinstance(x, (int, float)) for x in integrate_width))
    ):
        integrateMinMax = [
            [
                (integrate_center - np.abs(wid) / 2.0),
                (integrate_center + np.abs(wid) / 2.0),
            ]
            for wid in integrate_width
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]

    else:
        raise ValueError(
            "integrate_width must be 'full', int, float, or list of int or float; integrate_center must be int, float, or list of ints or floats"
        )

    remaining_dims = [x for x in data.dims if x != dim]
    if (
        len(remaining_dims) == 0
        and isinstance(integrate_center, (int, float))
        and isinstance(integrate_width, (int, float))
    ):
        remaining_dims = ["index"]
        remaining_coords = [np.array([0])]
    else:
        remaining_coords = [data.coords[x] for x in remaining_dims]

    if isinstance(data_new, list):
        data_integrals = [
            np.trapz(x.values, x=x.coords[dim], axis=index) for x in data_new
        ]
        if not all([isinstance(x, list) for x in [integrate_center, integrate_width]]):
            if isinstance(integrate_center, list) and not isinstance(
                integrate_width, list
            ):
                remaining_coords = [np.array(integrate_center)] + remaining_coords
                remaining_dims = ["center"] + remaining_dims
            elif isinstance(integrate_width, list) and not isinstance(
                integrate_center, list
            ):
                remaining_coords = [np.array(integrate_width)] + remaining_coords
                remaining_dims = ["width"] + remaining_dims
            data_values = np.array(data_integrals)
        elif isinstance(integrate_center, list) and isinstance(integrate_width, list):
            ind_dim = list(set(data.dims) - set([dim]))[0]
            ind_shape = data.shape[data.index(ind_dim)]
            remaining_coords = [
                np.array(integrate_center),
                np.array(integrate_width),
            ] + remaining_coords
            remaining_dims = ["center", "width"] + remaining_dims
            data_values = np.array(
                tuple([data_integrals for _ in range(len(data_integrals))])
            ).reshape(len(data_integrals), len(data_integrals), ind_shape)
    else:
        data_values = np.trapz(data_new.values, x=data_new.coords[dim], axis=index)

    if not isinstance(data_values, (list, np.ndarray)):
        data_values = [data_values]

    integrate_data = dnpdata(np.array(data_values), remaining_coords, remaining_dims)

    integrate_data.attrs["integrate_center"] = integrate_center
    integrate_data.attrs["integrate_width"] = integrate_width

    if type == "double":
        integrate_data.attrs["first_integral"] = first_int
        integrate_data.attrs["dim_coords"] = data.coords[dim]

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

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def zero_fill(
    all_data,
    dim="t2",
    factor=2,
):
    """
    Perform zero-filling (append) down given dimension

    Args:
        all_data (dnpdata): data object

    +------------------+------+-----------+--------------------------------------------------+
    | parameter        | type | default   | description                                      |
    +==================+======+===========+==================================================+
    | dim              | str  | 't2'      | dimension to zero-fill                           |
    +------------------+------+-----------+--------------------------------------------------+
    | factor           | int  | 2         | factor to increase dim with zeros                |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after zero-fill
    """

    data, isDict = return_data(all_data)
    index = data.dims.index(dim)

    factor = int(factor)
    if factor <= 0:
        factor = 1

    n_pts = data.shape[index] * factor
    data.coords[dim] = np.linspace(data.coords[dim][0], data.coords[dim][-1], num=n_pts)

    shape = list(data.shape)
    shape[index] = n_pts - data.shape[index]
    data.values = np.concatenate(
        (data.values, np.zeros(shape=tuple(shape))), axis=index
    )

    proc_parameters = {
        "dim": dim,
        "factor": factor,
    }

    proc_attr_name = "zero_fill"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def voigtian(x, amp, x0, sigma, gamma):
    """
    Voigt is a combintaion of Gaussian and Lorentzian lineshapes
    """
    z = ((x0 - x) + 1j * gamma) / (sigma * np.sqrt(2.0))
    fit = amp * np.real(wofz(z))
    return fit


def gaussian(x, amp, x0, sigma):
    return amp * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))


def lorentzian(x, amp, x0, gamma):
    return amp * gamma ** 2 / ((x - x0) ** 2 + gamma ** 2)
