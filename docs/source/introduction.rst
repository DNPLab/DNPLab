======================
Introduction to dnpLab
======================

In this introduction, we introduce the intended workflow for processing data with dnpLab. The dnpLab python package is intended to import data in standard spectrometer formats and convert this into a standard dnpLab format called "dnpdata".

Once the data is imported, it is added to a workspace which is a python dict-like class that stores multiple dnpdata objects. 

Workflow
========

Data container, workspace

.. figure:: _static/images/dnpLab_workflow.png
    :width: 400
    :alt: dnpLab Workflow
    :align: center

    dnpLab Workflow

Generalized data container (dnpdata)
====================================

Workspace Concept (dnpdata_collection)
======================================

