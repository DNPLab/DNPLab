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
     - 2.0.0 or higher
   * - SciPy
     - 1.14.0 or higher
   * - Matplotlib
     - 3.9.1 or higher
   * - h5py
     - 3.11.0 or higher


Hydration GUI
-------------
Due to cross-platfrom issues of PyQt5, the Hydration GUI is no longer part of DNPLab and needs to be installed separately. The GUI is part of the |HanLab| python package. Instructions on how to install and use the GUI can be found on the website of the |HanLab| python package.

Installing DNPLab 
=================

Installing using pip
--------------------
The easiest and most convenient way to install DNPLab is by using |pip|. In a terminal simply type the following command:

.. code-block:: bash

   $ pip install dnplab

or simply just:

.. code-block:: bash

   $ pip install dnplab


If you prefer to install DNPLab from the source code, check out our GitHub repository: |dnpLabGitLink|. The newest developments are merged into the *Development* branch.

Installing with a virtual enviroment
------------------------------------
Starting from Ubuntu 23.10 pip3 will issue a warning when trying to install dnplab from pypy.
It is recommended to not do a global install but use a virtual enviroment (venv).
If you do not have already have a virtual enviroment you can create a folder at a convenient location where the enviroment will be located.

In this example this will be in our home folder and the folder will be named DNPLab.
To create this enviroment use the command

.. code-block:: bash

   $ python3 -m venv ~/DNPLab

Note that you need to activate the venv when you want to use it and install packages via pip3.
You can activate the venv by sourcing the activate script that should now be created under ~/DNPLab/bin

.. code-block:: bash

  $ source ~/DNPLab/bin/activate

you have to do this everytime you start this enviroment, to ease this you can create an alias "dnplab" and add it to your .bash_aliases file

.. code-block:: bash

  $ echo "dnplab = 'source ~/DNPLab/bin/activate'" >> ~/.bash_aliases

to deactivate the enviroment just enter

.. code-block:: bash

  $ deactivate

into your terminal.

Confirm Successful Installation
-------------------------------
To confirm that your installation of DNPLab was successful type the following command:

.. code-block:: bash

    $ pip show dnplab

The output will look similar to this (note, the actual version and path to location depends on the local installation):

.. code-block:: bash

    Name: dnplab
    Version: 2.1.25
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
If you wish to install a specific version of DNPLab typ the following command in a terminal window:

.. code-block:: bash
    
    $ pip install dnplab==2.1.25

Install Preliminary Release
---------------------------
If you wish to use a pre-release version of DNPLab (downloaded from the GitHub repository) we recommend first uninstalling the current DNPLab version. Clone (or download or fork ...) the desired branch from the GitHub website. In a terminal window navigate into the directory that contains the dnplab folder and type the following command into the terminal window:

.. code-block:: bash
    
    $ pip install -e dnplab

Once you ran the above command, check the path and version of the package by running :code:`pip show dnplab`. If the version does not match the version of the checked-out branch, you may have to first uninstall DNPLab (:code:`pip uninstall dnplab`), then re-install the version you would like to use (:code:`pip install dnplab`) and then running (:code:`pip install -e dnplab`) if you would like to make your own changes to the code.

Upgrading DNPLab
================
To upgrade your currently installed version of DNPLab type the following command:

.. code-block:: bash

    $ pip install --upgrade dnplab


Uninstalling DNPLab
===================
The safest method to uninstall DNPLab is to use pip. Type the following command in a terminal window:
    
.. code-block:: bash
    
    $ pip uninstall dnplab
