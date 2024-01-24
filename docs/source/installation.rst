============
Installation
============

Required Packages
=================
The following packages are required to run DNPLab:

.. list-table::
   :widths: 40 60

   * - **Package**
     - **Version**
   * - NumPy
     - 1.19 or higher
   * - SciPy
     - 1.5 or higher
   * - Matplotlib
     - 3.3 or higher
   * - h5py
     - 2.10 or higher


Hydration GUI
-------------
Due to cross-platfrom issues of PyQt5, the Hydration GUI is no longer part of DNPLab and needs to be installed separately. The GUI is part of the |HanLab| python package. Instructions on how to install and use the GUI can be found on the website of the |HanLab| python package.

Installing DNPLab 
=================

Installing using pip
--------------------
The easiest and most convenient way to install DNPLab is by using |pip|. In a terminal simply type the following command:

.. code-block:: bash

   $ python -m pip install dnplab

or simply just:

.. code-block:: bash

   $ pip install dnplab


If you prefer to install DNPLab from the source code, check out our GitHub repository: |dnpLabGitLink|. The newest developments are merged into the *Development* branch.

Confirm Successful Installation
-------------------------------
To confirm that your installation of DNPLab was successful type the following command:

.. code-block:: bash

    $ pip show dnplab

The output will look similar to this (note, the actual version and path to location depends on the local installation):

.. code-block:: bash

    Name: dnplab
    Version: 1.0.3
    Summary: dnpLab - Bringing the Power of Python to DNP-NMR Spectroscopy
    Home-page: http://dnpLab.net
    Author: DNPLab Team
    Author-email: None
    License: MIT
    Location: /Path/to/Package
    Requires: numpy, scipy, matplotlib, h5py
    Required-by: 


Specify DNPLab Version to install
---------------------------------
If you wish to install a specific version of DNPLab typ the following comman in a terminal window:

.. code-block:: bash
    
    $ pip install dnplab==1.0.11

Install Preliminary Release
---------------------------
If you wish to use a pre-release version of DNPLab (downloaded from the GitHub repository) we recommend first uninstalling the current DNPLab version. Clone (or download or fork ...) the desired branch from the GitHub website. In a terminal window navigate into the directory that contains the setup.py file and type the following command into the terminal window:

.. code-block:: bash
    
    $ python setup.py develop

Once you ran the above command, check the path and version of the package by running :code:`pip show dnplab`. If the version does not match the version of the checked-out branch, you may have to first uninstall DNPLab (:code:`pip uninstall dnplab`), then re-install the version you would like to use (:code:`pip install dnplab`) and then running (:code:`python setup.py develop`) if you would like to make your own changes to the code.

Upgrading DNPLab
================
To upgrade your currently installed version of DNPLab type the following command:

.. code-block:: bash

    $ pip install dnplab --upgrade


Uninstalling DNPLab
===================
The safest method to uninstall DNPLab is to use pip. Type the following command in a terminal window:
    
.. code-block:: bash
    
    $ python -m pip uninstall dnplab
