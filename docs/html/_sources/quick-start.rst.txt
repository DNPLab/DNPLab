.. install:

==================
Quick-Start Guide
==================

Importing data
==============
.. code-block:: python

   import odnpLab
   odnpLab.odnpImport.bruker.importBruker()

Processing NMR Data
===================
.. code-block:: python

   paramDict = {'linewidth' : 20
                }
   dataDict = odnpLab.odnpNMR.removeOffset(dataDict,paramDict)
   dataDict = odnpLab.odnpNMR.window(dataDict,paramDict)
   dataDict = odnpLab.odnpNMR.fourierTransform(dataDict,paramDict)

Example Script
==============

