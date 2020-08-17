.. install:

==================
Quick-Start Guide
==================

Importing the Package
=====================

Once the package has been :ref:`installed via pip <installing>`, you should be able to import dnpLab in a python terminal/script.

.. code-block:: python

   import dnpLab as dnp


Importing data
==============
.. code-block:: python

   # Import module
   import dnpLab as dnp

   # import Topspin Data
   path = 'path/to/data'
   data = dnpLab.dnpImport.topspin.import_topspin(path)

   # create workspace for processing data
   workspace = dnp.create_workspace('raw', data)

Processing NMR Data
===================
.. code-block:: python

   # Remove DC offset from FID
   workspace = dnpLab.dnpNMR.remove_offset(workspace, {})
   # Apply Exponential Apodization to data
   workspace = dnpLab.dnpNMR.window(workspace, {})
   # Apply Fourier Transform to direct dimension by default (t2)
   workspace = dnpLab.dnpNMR.fourier_transform(workspace, {})

Example Script
==============

