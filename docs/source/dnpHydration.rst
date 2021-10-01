.. _dnpHydration:

============
dnpHydration
============

This module implements the process of calculating parameters that describe hydration dynamics using ODNP data as described in various studies from Songi Han and John Franck. The calculations follow:

J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
http://dx.doi.org/10.1016/j.pnmrs.2013.06.001

J.M. Franck, S. Han; Methods in Enzymology, Chapter 5, Volume 615, (2019) 131-175
https://doi.org/10.1016/bs.mie.2018.09.024

To use the dnpHydration module first create a dictionary with the necessary inputs and add it to a workspace as **'hydration_inputs'**. For example, start by defining the inputs dictionary,

.. code-block:: python

    import dnplab
    
    Enhancements = # list of signal enhancements
    Enhancement_powers = # list of powers in Watts corresponding to Enhancements
    T1s = # list of T1 values in seconds
    T1_powers = # list of powers in Watts corresponding to T1s
    
    inputs = {
              'E_array' : np.array(Enhancements),
              'E_powers' : np.array(Enhancement_powers),
              'T1_array' : np.array(T1s),
              'T1_powers' : np.array(T1_powers),
              'T10': 2.0, # T1 measured with power=0
              'T100': 2.5, # T1 measured with SL=0 and power=0
              'spin_C': 100, # spin concentration in micromolar
              'field': 350, # magnetic field in mT
              'smax_model': 'tethered', # choice of smax model
              'interpolate_method': 'second_order' # choice of interpolation method
              }
    

Now you can either create a workspace and add the dictionary under the key **'hydration_inputs'**,

.. code-block:: python

    workspace = dnplab.create_workspace('hydration_inputs', inputs)


Or add to an existing workspace,

.. code-block:: python

    workspace.add('hydration_inputs', inputs)


In rare cases the bulk water or second order T1 interpolation constants may need to be altered. This is not necessary for the dnpHydration module to operate, but if needed this can be done by adding the dictionary **'hydration_constants'** to the workspace. For example,

.. code-block:: python
    
    constants = {
                 'ksigma_bulk': 95.4, # bulk ksigma value
                 'krho_bulk': 353.4, # bulk krho value
                 'klow_bulk': 366, # bulk klow value
                 'tcorr_bulk': 54, # bulk tcorr value
                 'D_H2O': 2.3e-9, # bulk water diffusivity
                 'D_SL': 4.1e-10, # diffusivity of spin probe in bulk water
                 'delta_T1_water': 1 # change in water proton T1 due to microwaves
                 'T1_water': 2.5, # T1 of bulk water protons
                 'macro_C': 100, # concentration of macromolecule in uM
                 }

    workspace.add('hydration_constants', constants)


Next, pass the workspace to dnpHydration to perform calculations using,

.. code-block:: python

    hydration_results = dnplab.dnpHydration.hydration(workspace)


Or for in-place operation simply use,

.. code-block:: python

    dnplab.dnpHydration.hydration(workspace)


If returned, **hydration_results** is a dictionary that has the elements listed in the table below. If only in-place operation, the workspace will now contain a **'hydration_results'** dictionary. *Note: even if the dictionary is returned, the 'hydration_results' dictionary is still added to the workspace*

+-------------------+-------------+------------------------------------------------------------------------------------------+
| key               | type        | description                                                                              |
+===================+=============+==========================================================================================+
| uncorrected_Ep    | numpy array | fit to Equation 12 by varying coupling factor and p\ :sub:`1/2`                          |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| uncorrected_xi    | float       | coupling factor from uncorrected_Ep fit (unitless)                                       |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| interpolated_T1   | numpy array | interpolation of T1 measurements                                                         | 
+-------------------+-------------+------------------------------------------------------------------------------------------+
| ksigma_array      | numpy array | left side of Equation 42                                                                 |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| ksigma_fit        | numpy array | fit to Equation 42 by varying κ\ :sub:`σ` and p\ :sub:`1/2`                              |          
+-------------------+-------------+------------------------------------------------------------------------------------------+
| ksigma            | float       | cross-relaxivity, κ\ :sub:`σ`, (s\ :sup:`-1` M\ :sup:`-1`)                               |   
+-------------------+-------------+------------------------------------------------------------------------------------------+
| ksigma_stdd       | float       | standard deviation in κ\ :sub:`σ` (s\ :sup:`-1` M\ :sup:`-1`)                            |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| ksigma_bulk_ratio | float       | ratio of κ\ :sub:`σ` to bulk value (κ\ :sub:`σ,bulk` = 95.4 s\ :sup:`-1` M\ :sup:`-1`).  |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| krho              | float       | self-relaxivity, κ\ :sub:`ρ`, (s\ :sup:`-1` M\ :sup:`-1`)                                | 
+-------------------+-------------+------------------------------------------------------------------------------------------+
| krho_bulk_ratio   | float       | ratio of κ\ :sub:`ρ` to bulk value (κ\ :sub:`ρ,bulk` = 353.4 s\ :sup:`-1` M\ :sup:`-1`)  |          
+-------------------+-------------+------------------------------------------------------------------------------------------+
| klow              | float       | [(5/3)κ\ :sub:`ρ` - (7/3)κ\ :sub:`σ`]   (s\ :sup:`-1` M\ :sup:`-1`)                      |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| klow_bulk_ratio   | float       | ratio of κ\ :sub:`low` to bulk value (κ\ :sub:`low,bulk` = 366 s\ :sup:`-1` M\ :sup:`-1`)|          
+-------------------+-------------+------------------------------------------------------------------------------------------+
| coupling_factor   | float       | κ\ :sub:`σ` / κ\ :sub:`ρ` (unitless)                                                     |   
+-------------------+-------------+------------------------------------------------------------------------------------------+
| tcorr             | float       | translational correlation time, τ\ :sub:`corr` (ps), see Equations. 21-23                |
+-------------------+-------------+------------------------------------------------------------------------------------------+
| tcorr_bulk_ratio  | float       | ratio of τ\ :sub:`corr` to bulk value (τ\ :sub:`corr,bulk` = 54 ps)                      |          
+-------------------+-------------+------------------------------------------------------------------------------------------+
| Dlocal            | float       | local diffusivity, D\ :sub:`local`, (m\ :sup:`2`/s), see Equations 18-20                 |   
+-------------------+-------------+------------------------------------------------------------------------------------------+

If needed, access the results individually as follows,

.. code-block:: python
     
     interpolated_t1 = hydration_results['interpolated_T1']
     ksigma_array = hydration_results['ksigma_array']     
     ksigma = hydration_results['ksigma']
     coupling_factor = hydration_results['coupling_factor']
     etc.

Or,

.. code-block:: python
     
     interpolated_t1 = workspace['hydration_results']['interpolated_T1']
     ksigma_array = workspace['hydration_results']['ksigma_array']     
     ksigma = workspace['hydration_results']['ksigma']
     coupling_factor = workspace['hydration_results']['coupling_factor']
     etc.

For explanation of 'smax_model' see https://doi.org/10.1039/c0cp02126a. For explanations of 'interpolate_method' options or any of the equations used to calculate the hydration parameters refer to http://dx.doi.org/10.1016/j.pnmrs.2013.06.001 and https://doi.org/10.1016/bs.mie.2018.09.024.


Detailed Descriptions of Functions
==================================

.. automodule:: dnplab.dnpHydration
   :members:
   :show-inheritance:
   :member-order: bysource
