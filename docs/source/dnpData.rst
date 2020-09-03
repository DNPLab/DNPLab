.. _dnpdata:

=======
dnpData
=======

dnpdata Class Overview
======================

The dnpdata class is a flexible data container for N-dimensional data. The dnpdata class stores data, axes, parameters and processing information in a single object. 

The dnpdata class integrates some concepts from, and designed for ongoing compatibility with, `pySpecData <https://jmfrancklab.github.io/pyspecdata/>`_, an object-oriented spectral data processing package developed by the Franck Lab at Syracuse University, with ongoing collaborations between the development teams at Bridge12 Technologies, Inc. and the Han Lab at University of California, Santa Barabara.

dnpdata Attributes
==================

The attributes in the dnpdata object are named according to the convention by Pandas and xarray. 

+---------------+----------------------------+----------------------------------------------+
| **attribute** | **type**                   | **description**                              |
+---------------+----------------------------+----------------------------------------------+
| values        | numpy.ndarray              | Numpy array of data values                   |
+---------------+----------------------------+----------------------------------------------+
| dims          | list of str                | Names for each of the N-dimensions in values |
+---------------+----------------------------+----------------------------------------------+
| coords        | list of numpy.ndarray      | Axes values for each of the N-dimensions     |
+---------------+----------------------------+----------------------------------------------+
| attrs         | dict                       | Dictionary of miscellaneous                  |
+---------------+----------------------------+----------------------------------------------+
| proc_attrs    | list of tuples (str, dict) | List which stores each processing step       |
+---------------+----------------------------+----------------------------------------------+

dnpdata Examples
================

A dnpdata object can be defined as follows:

.. code-block:: python

   import dnplab as dnp
   import numpy as np

   x = np.r_[-10:10:100j]
   values = x**2.
   coords = [x]
   dims = ['x']

   data = dnp.dnpdata(values, coords, dims)

   print(data)


The dnpdata class has a number of methods for manipulating the data.

The dimensions can be renamed:

.. code-block:: python
   
   data.rename('x', 't')



A N-dimensional data set can be reshaped by it's dimension labels. If a data set has 

.. code-block:: python
   
   data.reorder(['z', 'x', 'y'])

Indexing
========

A number of methods can be used to index the data based on the coordinates.

.. code-block:: python
   
   data_slice = data['t', 0:10] # return first 10 points of data down 't' dimension

   data_slice = data['t', (0.1, 0.5)] # return data in range from 0.1 to 0.5 

   data_slice = data['t', 0.1] # for single float, return data at index nearest to 0.1 in time coordinates

   data_slice = data['t', 10] # for single int, return data at index 10 


dnpdata Methods
===============

.. autoclass:: dnplab.dnpdata
   :members:
   :show-inheritance:
   :member-order: bysource
   :inherited-members:


dnpdata_collection (Workspace)
==============================

To store multiple data objects, the user can create a workspace which is a dict like object for storing dnpdata objects.

.. code-block:: python

   import dnplab as dnp
   
   ws = dnp.create_workspace()
   ws['raw'] = data

Processing Buffer
=================

The workspace has an attribute called processing_buffer. The processing buffer indicates for functions which operate on the workspace, which dnpdata object should be operated on. By default, the processing_buffer is called "proc".

.. code-block:: python

   ws.copy('raw', 'proc') # copy some data into the default processing buffer

At any time, the processing buffer can be changed, however you need to make sure to move data into the processing buffer before any processing steps.

.. code-block:: python

   ws.processing_buffer = 'new_proc'

Saving the Workspace
====================

The workspace can be saved in h5 format with the saveh5 function:

.. code-block:: python

   dnplab.dnpImport.saveh5('test.h5', ws)


Loading the Workspace
=====================

A workspace can also be loaded with the loadh5 function.

.. code-block:: python

   dnplab.dnpImport.loadh5('test.h5')


dnpdata_collection Methods
==========================

.. autoclass:: dnplab.dnpdata_collection
   :members:
   :show-inheritance:
   :member-order: bysource

