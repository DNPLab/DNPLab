'''dnpTools
Collection of tools and functions useful to process DNP-NMR data
'''

import inspect
import math


def mrProperties(nucleus, *args):
# def mrProperty(nucleus, field):
    '''Return magnetic resonance property of specified isotope.
    This function is model after the Matlab function gmr written by Mirko Hrovat
    https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties

    Reference: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818.
    Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants .
       or  http://physics.nist.gov/PhysRefData/codata86/codata86.html
       or  http://www.isis.rl.ac.uk/neutronSites/constants.htm
    Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

    Args:

        nucleus: String defining the nucleus e.g. 1H, 13C, etc.

        args: If numerical value is given, it is interpreted as the B0 value in Tesla and Larmor frequency is returned. As string the following values are valid:

        gamma: Return Gyromagnetic Ration [radians/T/s]
        spin: Spin number of selected nucleus [1]
        qmom: Quadrupole moment [fm^2} (100 barns)
        natAbundance: Natural abundance [%]
        relSensitivity: Relative sensitiviy with respect to 1H at constant B0

    Returns:

        .. code-block:: python

            dnp.dnpTools.mrProperties('1H')
            26.7522128                          # 1H Gyromagnetic Ratio

            dnp.dnpTools.mrProperties('1H', 0.35)
            14902114.17018196                   # 1H Larmor Frequency at 0.35 T

            dnp.dnpTools.mrProperties('2H', 'qmom')
            0.286                               # Nuclear Quadrupole Moment

            value = dnp.dnpTools.mrProperties('6Li', 'natAbundance')
            7.59                                # Natural Abundance

            value = dnp.dnpTools.mrProperties('6Li', 'relSensitivity')
            0.000645                            # Relative sensitivity


    '''

    isotopeList = ['0e', '1H', '2H', '3H', '3He', '6Li', '7Li', '9Be', '10B', '11B', '13C', '14N', '15N', '17O', '19F', '21Ne', '23Na', '25Mg', '27Al', '29Si', '31P', '33S', '35Cl', '37C1', '39K', '40K', '41K', '43Ca', '45Sc', '47Ti', '49Ti', '50V', '51V'	, '53Cr', '55Mn', '57Fe', '59Co', '61Ni', '63Cu', '65Cu', '67Zn', '69Ga', '71Ga', '73Ge', '75As', '77Se', '79Br', '81Br', '83Kr', '85Rb', '87Rb', '87Sr', '89Y', '91Zr', '93Nb', '95Mo', '97Mo', '99Ru', '99Tc', '101Ru', '103Rh', '105Pd', '107Ag', '109Ag', '111Cd', '113Cd', '113In', '115Sn', '115In', '117Sn', '119Sn', '121Sb', '123Te', '123Sb', '125Te', '127I', '129Xe', '131Xe', '133Cs', '135Ba', '137Ba', '138La', '139La', '141Pr', '143Nd', '145Nd', '147Sm', '149Sm', '151Eu', '153Eu', '155Gd', '157Gd', '159Tb', '161Dy', '163Dy', '165Ho', '167Er', '169Tm', '171Yb', '173Yb', '175Lu', '176Lu', '177Hf', '179Hf', '181Ta', '183W', '185Re', '1870s', '187Re', '1890s', '191Ir', '193Ir', '195Pt', '197Au', '199Hg', '201Hg', '203Tl', '205Tl', '207Pb', '209Bi', '235U']

    spinList = [0.5, 0.5, 1, 0.5, 0.5, 1, 1.5, 1.5, 3, 1.5, 0.5, 1, 0.5, 2.5, 0.5, 1.5, 1.5, 2.5, 2.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 4, 1.5, 3.5, 3.5, 2.5, 3.5, 6, 3.5, 1.5, 2.5, 0.5, 3.5, 1.5, 1.5, 1.5, 2.5, 1.5, 1.5, 4.5, 1.5, 0.5, 1.5, 1.5, 4.5, 2.5, 1.5, 4.5, 0.5, 2.5, 4.5, 2.5, 2.5, 2.5, 4.5, 2.5, 0.5, 2.5, 0.5, 0.5, 0.5, 0.5, 4.5, 0.5, 4.5, 0.5, 0.5, 2.5, 0.5, 3.5, 0.5, 2.5, 0.5, 1.5, 3.5, 1.5, 1.5, 5, 3.5, 2.50, 3.50, 3.50, 3.50, 3.50, 2.50, 2.50, 1.50, 1.50, 1.50, 2.50, 2.50, 3.50, 3.50, 0.50, 0.50, 2.50, 3.50, 7.00, 3.5, 4.5, 3.5, 0.5, 2.5, 0.5, 2.5, 1.5, 1.5, 1.5, 0.5, 1.5, 0.5, 1.5, 0.5, 0.5, 0.5, 4.5, 3.50]

    # Gyromagnetic ratios are givin in 10^7*radians/T/s
    gmrList = [17608.59794, 26.7522128, 4.10662791, 28.5349779, -20.3801587, 3.9371709, 10.3977013, -3.759666, 2.8746786, 8.5847044, 6.728284, 1.9337792, -2.7126180, -3.62808, 25.18148, -2.11308, 7.0808493, -1.63887, 6.9762715, -5.319, 10.8394, 2.055685, 2.624198, 2.184368, 1.2500608, -1.5542854, 0.68606808, -1.803069, 6.5087973, -1.5105, -1.51095, 2.670649, 7.0455117, -1.5152, 6.6452546, 0.8680624, 6.332, -2.3948, 7.111789, 7.60435, 1.676688, 6.438855, 8.181171, -0.9360303, 4.596163, 5.1253857, 6.725616, 7.249776, -1.0331, 2.592705, 8.7864, -1.1639376, -1.3162791, -2.49743, 6.5674, -1.751, -1.788, -1.229, 6.046, -1.377, -0.8468, -1.23, -1.0889181, -1.2518634, -5.6983131, -5.9609155, 5.8845, -8.8013, 5.8972, -9.58879, -10.0317, 6.4435, -7.059098, 3.4892, -8.5108404, 5.389573, -7.3999, 2.209076, 3.5332539, 2.6755, 2.99295, 3.557239, 3.8083318, 8.1907, -1.457, -0.898, -1.115, -0.9192, 6.651, 2.9369, -0.82132, -1.0769, 6.431, -0.9201, 1.289, 5.71, -0.77157, -2.218, 4.7288, -1.3025, 3.0552, 2.1684, 1.086, -0.6821, 3.2438, 1.1282403, 6.1057, 0.6192895, 6.1682, 2.10713, 0.4812, 0.5227, 5.8385, 0.47306, 4.8457916, -1.788769, 15.5393338, 15.6921808, 5.58046, 4.375, -0.52] 

    Q = [0, 0, 0.286, 0, 0, -0.0808, -4.01, 5.288, 8.459, 4.059, 0, 2.044, 0, -2.55, 0, 10.15, 10.4, 19.94, 14.66, 0, 0, -6.78, -8.16, -6.43, 5.85, -7.3, 7.11, -4.08, -22, 30.2, 24.7, 21, -5.2, -15, 33, 0, 42, 16.2, -22, -20.4, 15, 17.1, 10.7, -19.6, 31.4, 0, 31.3, 26.2, 5.9, 27.6, 13.35, 33.5, 0, -17.6, -32, -2.2, 5.5, 7.9, -12.9, 45.7, 0, 66, 0, 0, 0, 0, 79.9, 0, 81, 0, 0, -36, 0, -49, 0, -71, 0, -11.4, -0.34, 16, 24.5, 45, 20, -5.89, -63, -33, -25.9, 7.4, 90.3, 241.2, 127, 135, 143.2, 250.7, 264.8, 358, 356.5, 0, 0, 280, 349, 497, 336.5, 379.3, 317, 0, 218, 0, 207, 85.6, 81.6, 75.1, 0, 54.7, 0, 38.6, 0, 0, 0, -51.6, 493.6]


    natAbundantList = [0, 99.9885, 0.0115, 0, 0.000137, 7.59, 92.41, 100, 19.9, 80.1, 1.07, 99.632, 0.368, 0.038, 100, 0.27, 100, 10, 100, 4.6832, 100, 0.76, 75.78, 24.22, 93.2581, 0.0117, 6.7302, 0.135, 100, 7.44, 5.41, 0.25, 99.75, 9.501, 100, 2.119, 100, 1.1399, 69.17, 30.83, 4.1, 60.108, 39.892, 7.73, 100, 7.63, 50.69, 49.31, 11.49, 72.17, 27.83, 7, 100, 11.22, 100, 15.92, 9.55, 12.76, [], 17.06, 100, 22.33, 51.839, 48.161, 12.8, 12.22, 4.29, 0.34, 95.71, 7.68, 8.59, 57.21, 0.89, 42.79, 7.07, 100, 26.44, 21.18, 100, 6.592, 11.232, 0.09, 99.91, 100, 12.2, 8.3, 14.99, 13.82, 47.81, 52.19, 14.8, 15.65, 100, 18.91, 24.9, 100, 22.93, 100, 14.28, 16.13, 97.41, 2.59, 18.6, 13.62, 99.988, 14.31, 37.4, 1.96, 62.6, 16.15, 37.3, 62.7, 33.832, 100, 16.87, 13.18, 29.524, 70.476, 22.1, 100, 0.72]   

    relativeSensitivityList = [0, 1, 1.11E-06, 0, 6.06E-07, 6.45E-04, 0.271, 1.39E-02, 3.95E-03, 0.132, 1.70E-04, 1.00E+03, 3.84E-06, 1.11E-05, 0.834, 6.65E-06, 9.27E-02, 2.68E-04, 0.207, 3.68E-04, 6.65E-02, 1.72E-05, 3.58E+03, 6.59E-04, 4.76E-04, 6.12E-07, 5.68E-06, 8.68E-06, 0.302, 1.56E-04, 2.05E-04, 1.39E-04, 0.383, 8.63E-05, 0.179, 7.24E-07, 0.278, 4.09E-05, 6.50E-02, 3.54E-02, 1.18E-04, 4.19E-02, 5.71E-02, 1.09E-04, 2.54E-02, 5.37E-04, 4.03E-02, 4.91E-02, 2.18E-04, 7.67E-03, 4.93E-02, 1.90E-04, 1.19E-04, 1.07E+03, 0.488, 5.21E-04, 3.33E-04, 1.44E-04, [], 2.71E-04, 3.17E-05, 2.53E-04, 3.50E-05, 4.94E-05, 1.24E-03, 1.35E-03, 1.51E-02, 1.21E-04, 0.338, 3.54E-03, 4.53E-03, 9.33E-02, 1.64E-04, 1.99E-02, 2.28E-03, 9.54E-02, 5.72E-03, 5.96E-04, 4.84E-02, 3.30E-04, 7.87E-04, 8.46E-05, 6.05E-02, [], [], [], [], [], [], [], [], [], [], [], [], [], [], 5.70E-04, 7.89E-04, [], [], [], 2.61E-04, 7.45E-05, 3.74E-02, 1.07E-05, 5.19E-02, 2.43E-07, 8.95E-02, 3.95E-04, 1.09E-05, 2.34E-05, 3.51E-03, 2.77E-05, 1.00E-03, 1.97E-04, 5.79E-02, 0.142, 2.01E-03, 0.144]

    if isinstance(nucleus, str):
        if nucleus in isotopeList:
            index = isotopeList.index(nucleus)
            gmr = gmrList[index]
        else:
            print("Isotope doesn't exist in list")
            return
    else:            
        print("ERROR: String expected")            

    if len(args) == 0:
        return gmr

    elif len(args) == 1:

        if isinstance(args[0], str):
            if args[0] == 'gamma':
                return gmrList[index] * 1e7
            if args[0] == 'spin':
                return spinList[index]
            elif args[0] == 'qmom':
                return Q[index]
            elif args[0] == 'natAbundance':
                return natAbundantList[index]
            elif args[0] == 'relSensitivity':
                return relativeSensitivityList[index]
            else:
                print("Keyword not recognize")

        else:
            vLarmor = args[0] * gmr * 1e7 / 2 / math.pi
            return vLarmor

    elif len(args) > 1:
        print("Too many input arguments")