.. dnpLab documentation master file, created by
   sphinx-quickstart on Fri Jan 10 17:26:18 2020.

.. .. image:: b12logo.png
..     :alt: Bridge12 Technologies, Inc
..     :target: http://www.bridge12.com

.. .. image:: hanlogo.png
..     :alt: Han Lab UCSB
..     :target: https://han.chem.ucsb.edu/

.. |GM116612| raw:: html

   <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=9896838&icde=51233599" target="_blank"> GM116612</a>

.. |B12TLink| raw:: html

   <a href="http://www.bridge12.com" target="_blank"> Bridge12 Technologies, Inc.</a>

.. |HanLabLink| raw:: html

   <a href="https://han.chem.ucsb.edu/" target="_blank"> Han Lab</a>

.. |FranckLabLink| raw:: html

   <a href="https://jmfrancklab.github.io/" target="_blank"> Franck Lab</a>


.. |dnpLabGitLink| raw:: html

   <a href="https://github.com/DNPLab/dnpLab" target="_blank"> DNPLab on Git Lab</a>

.. |dnpLabGitIssueTrackerLink| raw:: html

   <a href="https://github.com/DNPLab/dnpLab/issues" target="_blank"> DNPLab Git Issue Tracker</a>


=================
Welcome to DNPLab
=================

Welcome to the DNPLab documentation. DNPLab is an Open Source Python package for importing and processing Dynamic Nuclear Polarization (DNP) data. The aim of the project is to provide a free, turn-key processing package for easy processing and analysis of DMP-NMR data.

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
   quick-start
   examples
   dnpData
   dnpImport
   dnpNMR
   dnpFit
   dnpHydration
   dnpResults
   hydrationGUI

Index
=====

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Acknowledgements
================

Development of DNPLap is sponsored by grants by the United States and the Germany. In particular:

* |B12TLink| acknowledges support from the U.S. National Institutes of Health (NIH), grant |GM116612|.
* The |HanLabLink| acknowledges support from the U.S. National Science Foundation (NSF), Award No. CHE-1800596 and the German Science Foundation (DFG) as part of the Excellence Initiative under RESOLVDFG-EXC-2033 Project No. 390677874.
* The |FranckLabLink| acknowledges the received Startup Funds by Syracuse University.
