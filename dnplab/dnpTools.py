"""dnpTools
Collection of tools and functions useful to process DNP-NMR data
"""

from . import dnpNMR, dnpdata, dnpdata_collection
import numpy as np
import scipy.integrate
import copy

from .mrProperties import gmrProperties, radicalProperties


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


def exp_fit_func_1(x_axis, C1, C2, tau):
    return C1 + C2 * np.exp(-1.0 * x_axis / tau)


def exp_fit_func_2(x_axis, C1, C2, tau1, C3, tau2):
    return C1 + C2 * np.exp(-1.0 * x_axis / tau1) + C3 * np.exp(-1.0 * x_axis / tau2)


def baseline_fit(temp_coords, temp_data, type, order):

    if type == "polynomial":
        base_line = np.polyval(np.polyfit(temp_coords, temp_data, order), temp_coords)
    elif type == "exponential":
        temp_data = temp_data.real
        if order == 1:
            x0 = [temp_data[-1], temp_data[0], 1]
            out, cov = curve_fit(
                exp_fit_func_1, temp_coords, temp_data, x0, method="lm"
            )
            base_line = exp_fit_func_1(temp_coords, out[0], out[1], out[2])
        elif order == 2:
            x0 = [temp_data[-1], temp_data[0], 1, temp_data[0], 1]
            out, cov = curve_fit(
                exp_fit_func_2, temp_coords, temp_data, x0, method="lm"
            )
            base_line = exp_fit_func_2(
                temp_coords, out[0], out[1], out[2], out[3], out[4]
            )
        else:
            raise ValueError(
                "Use order=1 for mono-exponential, order=2 for bi-exponential"
            )

    else:
        raise TypeError("type must be either 'polynomial' or 'exponential'")

    return base_line


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
    +-----------------+------+---------------+---------------------------------------------------+
    | dim             | str  | 'f2'          | Dimension to apply baseline correction            |
    +-----------------+------+---------------+---------------------------------------------------+
    | type            | str  | 'polynomial'  | type of baseline fit                              |
    +-----------------+------+---------------+---------------------------------------------------+
    | order           | int  | 1             | polynomial order, or 1=mono 2=bi exponential      |
    +-----------------+------+---------------+---------------------------------------------------+
    | reference_slice | int  | None          | slice of 2D data used to define the baseline      |
    +-----------------+------+---------------+---------------------------------------------------+

    returns:
        all_data (dnpdata, dict): Baseline corrected data in container
        attributes: "baseline", baseline function
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
            bline = baseline_fit(
                data.coords[dim], data.values[:, reference_slice], type, order
            )
            for ix in range(len(data.coords[ind_dim])):
                data.values[:, ix] -= bline
        elif reference_slice is None:
            for ix in range(len(data.coords[ind_dim])):
                bline = baseline_fit(data.coords[dim], data.values[:, ix], type, order)
                data.values[:, ix] -= bline
        else:
            raise TypeError("invalid reference_slice")

    elif len(np.shape(data.values)) == 1:
        bline = baseline_fit(data.coords[dim], data.values, type, order)
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
        return all_data
    else:
        return data


def integrate(
    all_data, dim="f2", type="single", integrate_center=0, integrate_width="full"
):
    """Integrate data down given dimension

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+-------+---------+------------------------------+
    | parameter        | type  | default | description                  |
    +------------------+-------+---------+------------------------------+
    | dim              | str   | 't2'    | dimension to integrate       |
    +------------------+-------+---------+------------------------------+
    | integrate_center | float | 0       | center of integration window |
    +------------------+-------+---------+------------------------------+
    | integrate_width  | float | 100     | width of integration window  |
    +------------------+-------+---------+------------------------------+

    Returns:
        all_data (dnpdata,dict): Processed data

    Example::

        dnplab.dnpNMR.integrate(all_data)

    """

    data, isDict = return_data(all_data)
    index = data.index(dim)

    if type == "double":
        data.attrs["first_integral"] = scipy.integrate.cumtrapz(
            data.values, x=data.coords[dim], axis=index, initial=0
        )
        data.values = data.attrs["first_integral"]

    if integrate_width == "full":
        pass
    elif isinstance(integrate_width, int) or isinstance(integrate_width, float):
        integrateMin = integrate_center - np.abs(integrate_width) / 2.0
        integrateMax = integrate_center + np.abs(integrate_width) / 2.0
        data = data[dim, (integrateMin, integrateMax)]
    else:
        raise ValueError("integrate_width must be 'full', int, or float")

    data.values = np.trapz(data.values, x=data.coords[dim], axis=index)

    data.coords.pop(dim)

    proc_parameters = {
        "dim": dim,
        "integrate_center": integrate_center,
        "integrate_width": integrate_width,
    }

    proc_attr_name = "integrate"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def mr_properties(nucleus, *args):
    """Return magnetic resonance property of specified isotope.
    This function is model after the Matlab function gmr written by Mirko Hrovat
    https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties

    Reference: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818.
    Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants .
       or  http://physics.nist.gov/PhysRefData/codata86/codata86.html
       or  http://www.isis.rl.ac.uk/neutronSites/constants.htm
    Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

    Args:

        nucleus:        String defining the nucleus e.g. 1H, 13C, etc.

        args:           If numerical value is given, it is interpreted as the B0 value in Tesla and Larmor frequency is returned. As string the following values are valid:

        gamma:          Return Gyromagnetic Ration [radians/T/s]
        spin:           Spin number of selected nucleus [1]
        qmom:           Quadrupole moment [fm^2} (100 barns)
        natAbundance:   Natural abundance [%]
        relSensitivity: Relative sensitiviy with respect to 1H at constant B0
        moment:         Magnetic dipole moment in terms of the nuclear magneton, uN, |u|/uN = |gamma|*hbar[I(I + 1)]^1/2/uN , hbar=h/2pi.
        qlw:            quadrupolar line-width factor as defined by: Qlw = Q^2(2I + 3)/[I^2(2I � 1)]

    Returns:

        .. code-block:: python

            dnp.dnpTools.mrProperties('1H')
            26.7522128                          # 1H Gyromagnetic Ratio (10^7r/Ts)

            dnp.dnpTools.mrProperties('1H', 0.35)
            14902114.17018196                   # 1H Larmor Frequency at 0.35 T (Hz)

            dnp.dnpTools.mrProperties('2H', 'qmom')
            0.286                               # Nuclear Quadrupole Moment (fm^2)

            value = dnp.dnpTools.mrProperties('6Li', 'natAbundance')
            7.59                                # Natural Abundance (%)

            value = dnp.dnpTools.mrProperties('6Li', 'relSensitivity')
            0.000645                            # Relative sensitivity
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

    Returns:

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
        mwFreguency:    Microwave frequency in (Hz)
        dnpNuclues:     Nucleus for DNP-NMR experiments

        Returns:

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

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+-------+-----------+------------------------------+
    | parameter        | type  | default   | description                  |
    +------------------+-------+-----------+------------------------------+
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
        all_data (dnpdata,dict): data with signal_to_noise attribute added

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

    data.attrs["signal_to_noise"] = s_n

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data
