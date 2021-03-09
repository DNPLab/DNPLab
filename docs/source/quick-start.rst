.. install:

==================
Quick-Start Guide
==================

Importing the Package
=====================

Once the package has been :ref:`installed via pip <installing>`, you should be able to import DNPLab in a python terminal/script.

.. code-block:: python

   import dnplab as dnp


Importing data
==============
.. code-block:: python

   # Import module
   import dnplab as dnp

   # import Topspin Data
   path = 'path/to/data'
   data = dnp.dnpImport.load(path)

   # create workspace for processing data
   workspace = dnp.create_workspace('raw', data)

Processing NMR Data
===================
.. code-block:: python

   # Remove DC offset from FID
   dnp.dnpNMR.remove_offset(workspace)
   # Apply Exponential Apodization to data
   dnp.dnpNMR.window(workspace)
   # Apply Fourier Transform to direct dimension by default (t2)
   dnp.dnpNMR.fourier_transform(workspace)
