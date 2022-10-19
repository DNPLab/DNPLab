import numpy as np
from . import mr_properties

#######################################
# EPR Properties of Selected Radicals #
#######################################
#
radicalProperties = {}

# When adding new radicals, make sure to add a literature reference for each entry. Keys in dictonary should all be lower case.

# Free electron
# http://physics.nist.gov/constants
radicalProperties["gfree"] = [2.00231930436153, None, [0]]


# TEMPO in Toluene
# Table 1, page 243 in Lebedev et al., Bioactive Spinlabels, Springer Verlag, 1992
radicalProperties["tempo1"] = [[2.00980, 2.00622, 2.00220], "14N", [16.8, 20.5, 95.9]]

# TEMPO in Methanol
# Table 1, page 243 in Lebedev et al., Bioactive Spinlabels, Springer Verlag, 1992
radicalProperties["tempo2"] = [[2.00909, 2.00621, 2.00222], "14N", [20.2, 20.2, 102.1]]

# BDPA in Polystyrene
# Bennati et al., JMR, 1999. (2H couplings were scaled to 1H)
radicalProperties["bdpa"] = [[2.00263, 2.00260, 2.00257], "1H", [50.2, 34.5, 13.0]]


def radical_properties(name):
    """Return properties of different radicals. At the minimum the g value is returned. If available, large hyperfine couplings to a nucleus are returned. Add new properties or new radicals to radicalProperties.py

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

    Args:
        name (str): Name of the radical

    Returns:
        radicalProperties (dict): Principle g values and hyperfine coupling tensor

    Examples:

        Return g value of a free electron

        >>> radical_properties("gfree")

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
