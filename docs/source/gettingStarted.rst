=================
Getting Started
=================

This Getting Started page introduces some of the basic features of DNPLab. DNPLab is a python package for processing EPR, NMR, and DNP data.

.. code-block:: Python


    import dnplab as dnp # this will import the dnplab package
    import numpy as np # we will also use numpy


The DNPData object is central to the DNPLab package. The DNPData object stores the data, coordinates, dimensions, and other meta-data in a single object.

Let's create a DNPData object from scratch. We need to define a data array (values), dimension (dims), and coordinates (coords) for the data. 

.. code-block:: Python


    dims = ['t'] # create a list for the dimensions
    coords = np.r_[0:10:1000j] # create array from 0 to 1

    values = np.exp(-1*coords) # create a array of values for the data

    data = dnp.DNPData(values, dims, coords) # create a DNPData object

    print(data)



