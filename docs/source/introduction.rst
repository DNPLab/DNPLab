======================
Introduction to dnpLab
======================

In this section, we introduce the intended workflow for processing data with dnpLab. The dnpLab python package designed to import data from standard spectrometer formats.

The general workflow follows:

1. Import Data
2. Create Workspace
3. Process Data
4. Save data in h5 format

Workflow
========

.. figure:: _static/images/dnpLab_workflow.png
    :width: 400
    :alt: dnpLab Workflow
    :align: center

    dnpLab Workflow

Importing Data
--------------
The data is imported by the :ref:`dnpImport <import>` sub-package. In this sub-package there are modules for importing various spectrometer formats (e.g. :ref:`Topspin <topspin>`, :ref:`VnmrJ <vnmrj>`, :ref:`Prospa <prospa>`).

The data is imported as a :ref:`dnpdata <dnpdata>` object. The dnpdata object is a container for data (values), coordinates for each dimension (coords), dimension labels (dims), and experimental parameters (attrs). In addition, each processing step applied to the data is saved in the dnpdata object (stored as proc_attrs).

The dnpdata object is a flexible data format which can handle N-dimensional data and coordinates together.

Creating a workspace
--------------------
The workspace can be created with the "create_workspace" function in dnpLab. Once the data is imported, it is added to a workspace which is a python dict-like class that stores multiple dnpdata objects. A workspace is a collection of dnpdata objects and allows for raw and processed data to be saved in the same h5 file.


Processing Data
---------------
To process the data, we move or copy data into the processing buffer. The dnpLab workspace has the concept of a "processing_buffer". The processing buffer specifies the data which is meant for processing. Typically one will add raw data to the workspace and copy the data to the processing buffer.
dnpLab is primarily designed for processing DNP NMR data. For processing NMR data, this is performed in the :ref:`dnpNMR <dnpNMR>` module. 

Saving Data in h5 format
------------------------
Once the data is processed, the data can be saved in h5 format. This is done using the :ref:`h5 <h5>` module. The workspace can then be loaded, subsequent processing can be performed and the data can be saved again.

