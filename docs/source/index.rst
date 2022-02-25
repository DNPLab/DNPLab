=================
Welcome to DNPLab
=================

Welcome to the DNPLab documentation. DNPLab is an Open Source Python package for importing and processing Dynamic Nuclear Polarization (DNP) data. The aim of the project is to provide a free, turn-key processing package for easy processing and analysis of DNP-NMR data.

DNPLab is a collaborative project created by

* |B12TLink|
* The |HanLabLink| at University of California, Santa Barbara
* The |FranckLabLink| at Syracuse University
* License: `MIT License <https://en.wikipedia.org/wiki/MIT_License>`_

Authors: Timothy Keller, Thomas Casey, Yanxian Lin, John Franck, Thorsten Maly, Songi Han
    
The source code for the project is published here: |dnpLabGitLink|

Please report all issues related to dnpLab using the: |dnpLabGitIssueTrackerLink|

.. list-table::
   :widths: 60 40

   * - Current Release
     - |release|
   * - Documentation Build Date
     - |date|
   * - Author(s)
     - |author|

To check the version of DNPLab that is currently installed type:

.. code-block:: bash

    $ pip show dnplab


Features
========
* Import NMR spectra in various formats (Bruker - TopSpin, Varian - (Open) VnmrJ, Magritek - Kea) 
* Process NMR data
* Extract Hydration Dynamics Information
* Create Publication Quality Figures

.. toctree::
   :maxdepth: 2

   introduction
   install

.. quick-start
.. auto_examples/index
.. dnpData
.. dnpImport
.. dnpIO
.. dnpSave
.. dnpMath
.. dnpTools
.. dnpNMR
.. dnpFit
.. dnpHydration
.. dnpResults

Index
=====

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

Acknowledgements
================

Development of DNPLap is sponsored by grants by the United States and Germany. In particular:

* |B12TLink| acknowledges support from the U.S. National Institutes of Health (NIH), grant |GM116612|.
* The |HanLabLink| acknowledges support from the U.S. National Science Foundation (NSF), Award No. CHE-1800596 and the German Science Foundation (DFG) as part of the Excellence Initiative under RESOLVDFG-EXC-2033 Project No. 390677874.
* The |FranckLabLink| acknowledges the received Startup Funds by Syracuse University.
