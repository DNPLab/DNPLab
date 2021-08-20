================
DNPLab Utilities
================

DNPLab comes with many tools, scripts and functions to assist your day-to-day routines.

.. _mr_properties:

MR Properties
=============

The magnetic properties of nuclei of the periodic table can be accessed through the function mr_properties, a function of the dnpTools module. The function returns the  magnetic resonance properties of the specified isotope.  

The function is modeled after the Matlab function gmr written by Mirko Hrovat: https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties  

The properties are compiled using the following references: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818. Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants, http://physics.nist.gov/PhysRefData/codata86/codata86.html, or http://www.isis.rl.ac.uk/neutronSites/constants.htm. The Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.


Examples of MR Properties  
-------------------------

.. code-block:: python

    >>> dnp.dnpTools.mr_properties('1H')
    26.7522128

This returns the 1H gyromagnetic ratio of  = 26.7522128. The unit is radians per Tesla second.


.. code-block:: python

    >>> dnp.dnpTools.mr_properties('1H',0.35)
    14902114.17018196

This returns the 1H Larmor freqeuncy for a proton (1H) at a magnetic field strength of 0.35 T. The unit is in Hz.


.. code-block:: python

    >>> dnp.dnpTools.mr_properties('e',0.35)*1e-9
    9.80873391074068

This returns the (free) electron Larmor frequency at a magnetic field strength of 0.35 T. The returned value is in GHz.


.. code-block:: python

    >>> dnp.dnpTools.mr_properties('6Li','natAbundance')
    7.59

This returns the natural abundance of the selected isotope. The unit of the returned value is %.


.. code-block:: python

    >>> dnp.dnpTools.mr_properties('55Mn','spin')
    2.5

This returns the nuclear spin quantum number of the isotope.


The function can return more properties such as the quadrupole moment of nuclei with a spin quantum number of I > 1/2, or the relative sensitivity compared to proton (1H) detection. For more information take a look a the documentation of the function in the dnpTools section of this help.

