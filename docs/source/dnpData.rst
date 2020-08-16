=======
dnpData
=======

dnpdata Class Overview
======================

The dnpdata class is a flexible data container for N-dimensional data. The dnpdata class stores data, axes, parameters and processing information in a single object. 

Attributes
==========

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

.. automodule:: dnpLab.dnpdata
   :members:
   :show-inheritance:
   :member-order: bysource
   :inherited-members:


.. .. automodule:: dnpLab.dnpdata
..    :members:
..    :show-inheritance:
..    :member-order: bysource
..    :inherited-members:
