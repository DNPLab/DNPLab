"""dnpTools
Collection of tools and functions useful to process DNP-NMR data
"""

from . import dnpNMR, dnpdata, dnpdata_collection
import numpy as np
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


def calculate_enhancement(
    all_data,
    off_spectrum=1,
    on_spectra="all",
    integrate_center=0,
    integrate_width="full",
    method="integrate",
    dim="t2",
):
    """Calculate enhancement from DNP data

    Args:
        all_data (dnpdata, dict): data container

    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | parameter        | type                       | default     | description                                                          |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
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
    | dim              | str                        | "t2"        | dimension to integrate down or search down for max                   |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+

    Returns:
        all_data (dnpdata, dict): data object with "enhancement" added

    """

    orig_data, isDict = return_data(all_data)
    if dim == "t2":
        ind_dim = "t1"
    elif dim == "t1":
        ind_dim = "t2"
    elif dim == "f2":
        ind_dim = "f1"
    elif dim == "f1":
        ind_dim = "f2"

    if (
        isinstance(off_spectrum, dnpdata)
        or isinstance(off_spectrum, dnpdata_collection)
    ) and (
        isinstance(on_spectra, dnpdata) or isinstance(on_spectra, dnpdata_collection)
    ):

        data_off, is_ws_off = return_data(off_spectrum)
        data_on, is_ws_on = return_data(on_spectra)

        if integrate_width == "full":
            int_width_off = "full"
            int_width_on = "full"
        elif (
            isinstance(integrate_width, list) or isinstance(integrate_width, np.ndarray)
        ) and len(integrate_width) == 2:
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
            pass
        elif (
            isinstance(integrate_center, list)
            or isinstance(integrate_center, np.ndarray)
        ) and len(integrate_center) == 2:
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
            off_data1 = dnpNMR.integrate(
                data_off,
                dim=dim,
                integrate_center=int_center_off,
                integrate_width=int_width_off,
            )
            off_data = off_data1.values

            on_data1 = dnpNMR.integrate(
                data_on,
                dim=dim,
                integrate_center=int_center_on,
                integrate_width=int_width_on,
            )
            on_data = on_data1.values

        elif method == "amplitude":
            on_data = []
            if integrate_center == "max":
                off_data = data_off.values[np.argmax(abs(data_off.values))]
                if len(data_on.shape) == 1:
                    on_data.append(data_on.values[np.argmax(abs(data_on.values))])
                else:
                    for indx in range(data_on.shape[-1]):
                        on_data.append(
                            data_on.values[np.argmax(abs(data_on.values[indx])), indx]
                        )
            else:
                off_data = data_off.values[int_center_off]
                if len(data_on.shape) == 1:
                    on_data.append(data_on.values[int_center_on])
                else:
                    for indx in range(data_on.shape[-1]):
                        on_data.append(data_on.values[int_center_on, indx])

        if data_on.ndim == 2:
            enh_coords_on = data_on.coords[ind_dim]
        else:
            enh_coords_on = np.array(range(data_on.shape[-1]))

    elif isinstance(off_spectrum, int) and on_spectra == "all":

        enh_data = copy.deepcopy(all_data)
        if orig_data.ndim == 1:
            raise ValueError("data is 1D, enhancement will be equal to 1 !!")

        if (
            isinstance(integrate_width, list) or isinstance(integrate_width, np.ndarray)
        ) and len(integrate_width) > 1:
            raise ValueError("supply a single value for integrate_width, or use 'full'")
        elif isinstance(integrate_width, str) and integrate_width != "full":
            raise ValueError("the only allowed integrate_width string is 'full'")

        if (
            isinstance(integrate_center, list)
            or isinstance(integrate_center, np.ndarray)
        ) and len(integrate_center) > 1:
            raise ValueError("supply a single value for integrate_center, or use 'max'")
        elif isinstance(integrate_center, str) and integrate_center != "max":
            raise ValueError("the only allowed integrate_center string is 'max'")

        if off_spectrum == 0:
            off_spectrum = 1

        if method == "integrate":

            dnpNMR.integrate(
                enh_data,
                dim=dim,
                integrate_center=integrate_center,
                integrate_width=integrate_width,
            )
            data, _isDict = return_data(enh_data)
            data_1 = data.values

        elif method == "amplitude":
            data_1 = []
            if integrate_center == "max":
                for indx in range(orig_data.shape[-1]):
                    data_1.append(
                        orig_data.values[np.argmax(abs(orig_data.values[indx])), indx]
                    )
            else:
                for indx in range(orig_data.shape[-1]):
                    data_1.append(orig_data.values[integrate_center, indx])

        off_data = data_1[off_spectrum - 1]
        if off_spectrum == 1:
            on_data = data_1[1:]
            enh_coords_on = orig_data.coords[ind_dim][1:]
        elif off_spectrum > 1:
            on_data_1 = data_1[: off_spectrum - 1]
            on_coords_1 = orig_data.coords[ind_dim][: off_spectrum - 1]
            on_data_2 = data_1[off_spectrum:]
            on_coords_2 = orig_data.coords[ind_dim][off_spectrum:]
            on_data = np.concatenate((on_data_1, on_data_2))
            enh_coords_on = np.concatenate((on_coords_1, on_coords_2))
    else:
        raise TypeError(
            "the given combination of data, off_spectrum, and on_spectra is not valid"
        )

    enh = np.real(np.array(on_data) / np.array(off_data))
    enhancementData = dnpdata(enh, [enh_coords_on], [ind_dim])

    if isDict:
        all_data["enhancement"] = enhancementData
        return all_data
    else:
        return enhancementData


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
