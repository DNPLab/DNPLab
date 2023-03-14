import numpy as _np

#####################################
# Gyromagnetic Properties of nuclei #
#####################################
#
#     Reference: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818.
#     Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants .
#        or  http://physics.nist.gov/PhysRefData/codata86/codata86.html
#        or  http://www.isis.rl.ac.uk/neutronSites/constants.htm
#     Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

gmrProperties = {}
# gmrProperties['label'] = ['spin', 'GMR 10^7r/Ts', 'Qmom fm^2', '%NatAbun', 'relH', 'moment', 'Qlw']
gmrProperties["0e"] = [0.5, 17608.59794, 0, 0, 0, 1838.28200037, 0]
gmrProperties["1H"] = [0.5, 26.7522128, 0, 99.9885, 1, 4.83735357, 0]
gmrProperties["2H"] = [1, 4.10662791, 0.286, 0.0115, 1.11e-06, 1.21260077, 0.41]
gmrProperties["3H"] = [0.5, 28.5349779, 0, 0, 0, 5.159714367, 0]
gmrProperties["3He"] = [0.5, -20.3801587, 0, 0.000137, 6.06e-07, -3.685154336, 0]
gmrProperties["6Li"] = [1, 3.9371709, -0.0808, 7.59, 0.000645, 1.1625637, 0.033]
gmrProperties["7Li"] = [1.5, 10.3977013, -4.01, 92.41, 0.271, 4.20407505, 21]
gmrProperties["9Be"] = [1.5, -3.759666, 5.288, 100, 0.0139, -1.520136, 37]
gmrProperties["10B"] = [3, 2.8746786, 8.459, 19.9, 0.00395, 2.0792055, 14]
gmrProperties["11B"] = [1.5, 8.5847044, 4.059, 80.1, 0.132, 3.4710308, 22]
gmrProperties["13C"] = [0.5, 6.728284, 0, 1.07, 0.00017, 1.216613, 0]
gmrProperties["14N"] = [1, 1.9337792, 2.044, 99.632, 1000, 0.57100428, 21]
gmrProperties["15N"] = [0.5, -2.71261804, 0, 0.368, 3.84e-06, -0.49049746, 0]
gmrProperties["17O"] = [2.5, -3.62808, -2.558, 0.038, 1.11e-05, -2.24077, 2.1]
gmrProperties["19F"] = [0.5, 25.18148, 0, 100, 0.834, 4.553333, 0]
gmrProperties["21Ne"] = [1.5, -2.11308, 10.155, 0.27, 6.65e-06, -0.854376, 140]
gmrProperties["23Na"] = [1.5, 7.0808493, 10.4, 100, 0.0927, 2.8629811, 140]
gmrProperties["25Mg"] = [2.5, -1.63887, 19.94, 10, 0.000268, -1.0122, 130]
gmrProperties["27Al"] = [2.5, 6.9762715, 14.66, 100, 0.207, 4.3086865, 69]
gmrProperties["29Si"] = [0.5, -5.319, 0, 4.6832, 0.000368, -0.96179, 0]
gmrProperties["31P"] = [0.5, 10.8394, 0, 100, 0.0665, 1.95999, 0]
gmrProperties["33S"] = [1.5, 2.055685, -6.78, 0.76, 1.72e-05, 0.8311696, 61]
gmrProperties["35Cl"] = [1.5, 2.624198, -8.165, 75.78, 3580, 1.061035, 89]
gmrProperties["37C1"] = [1.5, 2.184368, -6.435, 24.22, 0.000659, 0.8831998, 55]
gmrProperties["39K"] = [1.5, 1.2500608, 5.85, 93.2581, 0.000476, 0.50543376, 46]
gmrProperties["40K"] = [4, -1.5542854, -7.3, 0.0117, 6.12e-07, -1.4513203, 5.2]
gmrProperties["41K"] = [1.5, 0.68606808, 7.11, 6.7302, 5.68e-06, 0.27739609, 67]
gmrProperties["43Ca"] = [3.5, -1.803069, -4.08, 0.135, 8.68e-06, -1.494067, 2.3]
gmrProperties["45Sc"] = [3.5, 6.5087973, -22, 100, 0.302, 5.3933489, 66]
gmrProperties["47Ti"] = [2.5, -1.5105, 30.2, 7.44, 0.000156, -0.93294, 290]
gmrProperties["49Ti"] = [3.5, -1.51095, 24.7, 5.41, 0.000205, -1.25201, 83]
gmrProperties["50V"] = [6, 2.670649, 21, 0.25, 0.000139, 3.613757, 17]
gmrProperties["51V"] = [3.5, 7.0455117, -5.2, 99.75, 0.383, 5.8380835, 3.7]
gmrProperties["53Cr"] = [1.5, -1.5152, -15, 9.501, 8.63e-05, -0.61263, 300]
gmrProperties["55Mn"] = [2.5, 6.6452546, 33, 100, 0.179, 4.1042437, 350]
gmrProperties["57Fe"] = [0.5, 0.8680624, 0, 2.119, 7.24e-07, 0.1569636, 0]
gmrProperties["59Co"] = [3.5, 6.332, 42, 100, 0.278, 5.247, 240]
gmrProperties["61Ni"] = [1.5, -2.3948, 16.2, 1.1399, 4.09e-05, -0.96827, 350]
gmrProperties["63Cu"] = [1.5, 7.111789, -22, 69.17, 0.065, 2.8754908, 650]
gmrProperties["65Cu"] = [1.5, 7.60435, -20.4, 30.83, 0.0354, 3.07465, 550]
gmrProperties["67Zn"] = [2.5, 1.676688, 15, 4.1, 0.000118, 1.035556, 72]
gmrProperties["69Ga"] = [1.5, 6.438855, 17.1, 60.108, 0.0419, 2.603405, 390]
gmrProperties["71Ga"] = [1.5, 8.181171, 10.7, 39.892, 0.0571, 3.307871, 150]
gmrProperties["73Ge"] = [4.5, -0.9360303, -19.6, 7.73, 0.000109, -0.9722881, 28]
gmrProperties["75As"] = [1.5, 4.596163, 31.4, 100, 0.0254, 1.858354, 1300]
gmrProperties["77Se"] = [0.5, 5.1253857, 0, 7.63, 0.000537, 0.92677577, 0]
gmrProperties["79Br"] = [1.5, 6.725616, 31.3, 50.69, 0.0403, 2.719351, 1300]
gmrProperties["81Br"] = [1.5, 7.249776, 26.2, 49.31, 0.0491, 2.931283, 920]
gmrProperties["83Kr"] = [4.5, -1.0331, 5.9, 11.49, 0.000218, -1.07311, 50]
gmrProperties["85Rb"] = [2.5, 2.592705, 27.6, 72.17, 0.00767, 1.6013071, 240]
gmrProperties["87Rb"] = [1.5, 8.7864, 13.35, 27.83, 0.0493, 3.552582, 240]
gmrProperties["87Sr"] = [4.5, -1.1639376, 33.5, 7, 0.00019, -1.2090236, 83]
gmrProperties["89Y"] = [0.5, -1.3162791, 0, 100, 0.000119, -0.23801049, 0]
gmrProperties["91Zr"] = [2.5, -2.49743, -17.6, 11.22, 1070, -1.54246, 99]
gmrProperties["93Nb"] = [4.5, 6.5674, -32, 100, 0.488, 6.8217, 76]
gmrProperties["95Mo"] = [2.5, -1.751, -2.2, 15.92, 0.000521, -1.082, 1.5]
gmrProperties["97Mo"] = [2.5, -1.788, 5.5, 9.55, 0.000333, -1.105, 210]
gmrProperties["99Ru"] = [2.5, -1.229, 7.9, 12.76, 0.000144, -0.7588, 20]
gmrProperties["99Tc"] = [4.5, 6.046, -12.9, None, None, 6.281, 12]
gmrProperties["101Ru"] = [2.5, -1.377, 45.7, 17.06, 0.000271, -0.8505, 670]
gmrProperties["103Rh"] = [0.5, -0.8468, 0, 100, 3.17e-05, -0.1531, 0]
gmrProperties["105Pd"] = [2.5, -1.23, 66, 22.33, 0.000253, -0.76, 1400]
gmrProperties["107Ag"] = [0.5, -1.0889181, 0, 51.839, 3.5e-05, -0.19689893, 0]
gmrProperties["109Ag"] = [0.5, -1.2518634, 0, 48.161, 4.94e-05, -0.22636279, 0]
gmrProperties["111Cd"] = [0.5, -5.6983131, 0, 12.8, 0.00124, -1.0303729, 0]
gmrProperties["113Cd"] = [0.5, -5.9609155, 0, 12.22, 0.00135, -1.0778568, 0]
gmrProperties["113In"] = [4.5, 5.8845, 79.9, 4.29, 0.0151, 6.1124, 470]
gmrProperties["115Sn"] = [0.5, -8.8013, 0, 0.34, 0.000121, -1.5915, 0]
gmrProperties["115In"] = [4.5, 5.8972, 81, 95.71, 0.338, 6.156, 490]
gmrProperties["117Sn"] = [0.5, -9.58879, 0, 7.68, 0.00354, -1.73385, 0]
gmrProperties["119Sn"] = [0.5, -10.0317, 0, 8.59, 0.00453, -1.81394, 0]
gmrProperties["121Sb"] = [2.5, 6.4435, -36, 57.21, 0.0933, 3.9796, 410]
gmrProperties["123Te"] = [0.5, -7.059098, 0, 0.89, 0.000164, -1.276431, 0]
gmrProperties["123Sb"] = [3.5, 3.4892, -49, 42.79, 0.0199, 2.8912, 330]
gmrProperties["125Te"] = [0.5, -8.5108404, 0, 7.07, 0.00228, -1.538936, 0]
gmrProperties["127I"] = [2.5, 5.389573, -71, 100, 0.0954, 3.32871, 1600]
gmrProperties["129Xe"] = [0.5, -7.3999, 0, 26.44, 0.00572, -1.347494, 0]
gmrProperties["131Xe"] = [1.5, 2.209076, -11.4, 21.18, 0.000596, 0.8931899, 170]
gmrProperties["133Cs"] = [3.5, 3.5332539, -0.343, 100, 0.0484, 2.9277407, 0.016]
gmrProperties["135Ba"] = [1.5, 2.6755, 16, 6.592, 0.00033, 1.08178, 340]
gmrProperties["137Ba"] = [1.5, 2.99295, 24.5, 11.232, 0.000787, 1.21013, 800]
gmrProperties["138La"] = [5, 3.557239, 45, 0.09, 8.46e-05, 4.068095, 120]
gmrProperties["139La"] = [3.5, 3.8083318, 20, 99.91, 0.0605, 3.155677, 54]
gmrProperties["141Pr"] = [2.5, 8.1907, -5.89, 100, None, 5.0587, None]
gmrProperties["143Nd"] = [3.5, -1.457, -63, 12.2, None, -1.208, None]
gmrProperties["145Nd"] = [3.5, -0.898, -33, 8.3, None, -0.744, None]
gmrProperties["147Sm"] = [3.5, -1.115, -25.9, 14.99, None, -0.9239, None]
gmrProperties["149Sm"] = [3.5, -0.9192, 7.4, 13.82, None, -0.7616, None]
gmrProperties["151Eu"] = [2.5, 6.651, 90.3, 47.81, None, 4.1078, None]
gmrProperties["153Eu"] = [2.5, 2.9369, 241.2, 52.19, None, 1.8139, None]
gmrProperties["155Gd"] = [1.5, -0.82132, 127, 14.8, None, -0.33208, None]
gmrProperties["157Gd"] = [1.5, -1.0769, 135, 15.65, None, -0.4354, None]
gmrProperties["159Tb"] = [1.5, 6.431, 143.2, 100, None, 2.6, None]
gmrProperties["161Dy"] = [2.5, -0.9201, 250.7, 18.91, None, -0.5683, None]
gmrProperties["163Dy"] = [2.5, 1.289, 264.8, 24.9, None, 0.7958, None]
gmrProperties["165Ho"] = [3.5, 5.71, 358, 100, None, 4.732, None]
gmrProperties["167Er"] = [3.5, -0.77157, 356.5, 22.93, None, -0.63935, None]
gmrProperties["169Tm"] = [0.5, -2.218, 0, 100, 0.00057, -0.4011, 0]
gmrProperties["171Yb"] = [0.5, 4.7288, 0, 14.28, 0.000789, 0.85506, 0]
gmrProperties["173Yb"] = [2.5, -1.3025, 280, 16.13, None, -0.80446, None]
gmrProperties["175Lu"] = [3.5, 3.0552, 349, 97.41, None, 2.5316, None]
gmrProperties["176Lu"] = [7, 2.1684, 497, 2.59, None, 3.388, None]
gmrProperties["177Hf"] = [3.5, 1.086, 336.5, 18.6, 0.000261, 0.8997, 15000]
gmrProperties["179Hf"] = [4.5, -0.6821, 379.3, 13.62, 7.45e-05, -0.7085, 11000]
gmrProperties["181Ta"] = [3.5, 3.2438, 317, 99.988, 0.0374, 2.6879, 14000]
gmrProperties["183W"] = [0.5, 1.1282403, 0, 14.31, 1.07e-05, 0.20400919, 0]
gmrProperties["185Re"] = [2.5, 6.1057, 218, 37.4, 0.0519, 3.771, 15000]
gmrProperties["1870s"] = [0.5, 0.6192895, 0, 1.96, 2.43e-07, 0.1119804, 0]
gmrProperties["187Re"] = [2.5, 6.1682, 207, 62.6, 0.0895, 3.8096, 14000]
gmrProperties["1890s"] = [1.5, 2.10713, 85.6, 16.15, 0.000395, 0.85197, 9800]
gmrProperties["191Ir"] = [1.5, 0.4812, 81.6, 37.3, 1.09e-05, 0.1946, 8900]
gmrProperties["193Ir"] = [1.5, 0.5227, 75.1, 62.7, 2.34e-05, 0.2113, 7500]
gmrProperties["195Pt"] = [0.5, 5.8385, 0, 33.832, 0.00351, 1.0557, 0]
gmrProperties["197Au"] = [1.5, 0.47306, 54.7, 100, 2.77e-05, 0.191271, 4000]
gmrProperties["199Hg"] = [0.5, 4.8457916, 0, 16.87, 0.001, 0.87621937, 0]
gmrProperties["201Hg"] = [1.5, -1.788769, 38.6, 13.18, 0.000197, -0.7232483, 2000]
gmrProperties["203Tl"] = [0.5, 15.5393338, 0, 29.524, 0.0579, 2.80983305, 0]
gmrProperties["205Tl"] = [0.5, 15.6921808, 0, 70.476, 0.142, 2.83747094, 0]
gmrProperties["207Pb"] = [0.5, 5.58046, 0, 22.1, 0.00201, 1.00906, 0]
gmrProperties["209Bi"] = [4.5, 4.375, -51.6, 100, 0.144, 4.5444, 200]
gmrProperties["235U"] = [3.5, -0.52, 493.6, 0.72, None, -0.43, None]


def mr_properties(nucleus, *args):
    """Return magnetic resonance property of specified isotope.

    This function is modeled after the Matlab function gmr written by Mirko Hrovat: https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties

    Also see: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818. Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants, http://physics.nist.gov/PhysRefData/codata86/codata86.html, or http://www.isis.rl.ac.uk/neutronSites/constants.htm. Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

    Args:

        nucleus (str): '1H', '2H', '6Li', '13C', 14N', etc.
        B0 (float): (optional) B0 field in (mT)

        Additional flags (see examples below)

            gamma: Return Gyromagnetic Ration (radians/T/s)

            spin: Return spin number of selected nucleus

            qmom: Resturn quadrupole moment [fm^2] (100 barns)

            natAbundance: Return natural abundance (%)

            relSensitivity: Return relative sensitiviy with respect to 1H at constant B0

            moment: Return magnetic dipole moment, abs(u)/uN = abs(gamma)*hbar[I(I + 1)]^1/2/uN

            qlw: Return quadrupolar line-width factor, Qlw = Q^2(2I + 3)/[I^2(2I + 1)]

    Examples:
        .. code-block:: python

            dnp.dnpTools.mr_Properties('1H') = 26.7522128 # 1H Gyromagnetic Ratio (10^7r/Ts)

            dnp.dnpTools.mr_Properties('1H', 0.35) = 14902114.17018196 # 1H Larmor Freq at .35 T (Hz)

            dnp.dnpTools.mr_Properties('2H', 'qmom') = 0.286 # Nuclear Quadrupole Moment (fm^2)

            dnp.dnpTools.mr_Properties('6Li', 'natAbundance') = 7.59 # % Natural Abundance

            dnp.dnpTools.mr_Properties('6Li', 'relSensitivity') = 0.000645 # Relative sensitivity

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
            vLarmor = args[0] * gmr * 1e7 / 2 / _np.pi
            return vLarmor

    elif len(args) == 2:
        if args[1] == True:
            print(" ")
            print("Nucleus                    : ", nucleus)
            print("Spin                       : ", gmrProperties.get(nucleus)[0])
            print(
                "Gyromagnetic Ratio [kHz/T] : %5.2f"
                % (gmrProperties.get(nucleus)[1] * 10 / 2 / _np.pi)
            )
            print(
                "Natural Abundance      [%%] : %5.2f" % (gmrProperties.get(nucleus)[3])
            )
            print("")

    elif len(args) > 2:
        print("Too many input arguments")
