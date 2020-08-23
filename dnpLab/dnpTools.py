'''Collection of tools and functions useful to process DNP-NMR data'''


def mrProperty(nucleus, field):
    '''Return magnetic resonance property of specified isotope.
    This function is model after the Matlab function gmr.

           gamma   = Gyromagnetic ratio of selected nucleus.
           spin    = Spin number of selected nucleus
           qmom    = Quadrupole moment in fm^2 (100 barns).
           natabun = Natural abundance in mole percent.
           relsens = Relative sensitiviy with respect to 1H at constant Bo.
           label   = String label for the specified nucleus.
           moment  = Magnetic dipole moment in terms of the nuclear magneton.
                       value of nuclear magneton = 5.0507866e-27  J/T
           qlw     = Quadrupolar line-width factor which combined with 
                       "relsens" is an estimator for the ease of detection.

    Args:

    Returns:

    '''



    # '0e'    0.5     17608.59794     0           0           0           1838.28200037   0
    # '1H'    0.5     26.7522128      0           99.9885     1           4.83735357      0
    # '2H'	1    	4.10662791      0.286       0.0115      1.11E-06    1.21260077      0.41
    # '3H'	0.5 	28.5349779      0           0           0           5.159714367     0
    # '3He'   0.5     -20.3801587     0           0.000137    6.06E-07    -3.685154336	0
    # '6Li'   1    	3.9371709       -0.0808     7.59        6.45E-04    1.1625637       0.033
 

    return 12345

