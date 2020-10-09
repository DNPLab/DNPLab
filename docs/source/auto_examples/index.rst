:orphan:



.. _sphx_glr_auto_examples:

DNPLab Examples
===============

Below is a gallery of examples


.. |ExamplesLink| raw:: html

   <a href="https://github.com/DNPLab/dnpLab/tree/master/examples" target="_blank"> Link to Examples</a>

.. |DataLink| raw:: html

   <a href="https://github.com/DNPLab/dnpLab/tree/master/data" target="_blank"> Link to Data</a>


DNPLab comes with many example scripts to demonstrate how the package can be used to import data from different spectrometer platforms, process NMR data, and extract enhancement data or hydration information. The example scripts are located in the *examples* folder using sample data located in the *data* folder.

If you installed DNPlab using pip you can download the example scripts and data from the GitHub repository:

.. list-table::
   :widths: 50 50

   * - Example Scripts:
     - |ExamplesLink|
   * - Example Data:
     - |DataLink|





.. Import Data and Process FID (Bruker Format)
.. ===========================================
.. This example uses the example script: *example_process_1Dbruker.py*. The script demonstrates the following features of DNPLab:

.. #. Load a single FID (Bruker format)
.. #. Perform an offset correction
.. #. Apply apodization to the FID
.. #. Perform a Fourier transformation
.. #. Phase correct the resulting spectrum


.. If you installed DNPlab using pip, you only need to import NumPy and DNPLab using:

.. .. code-block:: python

..    import numpy as np
..    import dnplab


.. If you downloaded DNPLab via GitHub and have not installed, you must add the directory for DNPLab to the system path before importing DNPLab. Instead of the above, use:

.. .. code-block:: python

..    import sys
..    sys.path.append('path/to/dnplab/package')

..    import numpy as np
..    import dnplab


.. In the next step load a single FID in Bruker format:

.. .. code-block:: python

..     path = 'path/to/data/topspin/'
..     folder = 20

..     data = dnplab.dnpImport.topspin.import_topspin(path,folder)

.. The topspin import module requires the path and the folder number.

.. In the next step the workspace is set up and the imported data is added to the *raw* workspace and the same data is copied to the *proc* workspace.

.. .. code-block:: python

..     ws = dnplab.create_workspace()
..     ws.add('raw', data)
..     ws.copy('raw', 'proc')

.. .. note::

..     When working with DNPLab one of the first steps is to copy the imported data to the *raw* workspace. That way the raw data and all its attributes will be always accessible to the user. When saving data with DNPLab the raw data is saved together with the processed data. DNPLab uses the h5 format to store data. 

.. In the following steps, the FID is processed and the spectrum is plotted.

.. .. code-block:: python

..     dnplab.dnpNMR.remove_offset(ws,{})
..     dnplab.dnpNMR.window(ws,{'linewidth' : 10})
..     dnplab.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
..     dnplab.dnpNMR.autophase(ws,{})


.. In this example, a baseline correction is performed (dnpNMR.remove_offset), apodization is applied to the FID (dnpNMR.window) and a line broadening of 10 Hz is applied. The next step is to Fourier transform the FID (dnpNMR.fourier_transform) and phase the spectrum (dnpNMR.autophase).

.. To plot the NMR spectrum: 

.. .. code-block:: python

..     dnplab.dnpResults.figure()
..     dnplab.dnpResults.plot(ws['proc'].real)
..     dnplab.dnpResults.xlim([-35,50])
..     dnplab.dnpResults.plt.xlabel('Chemical Shift [ppm]')
..     dnplab.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
..     dnplab.dnpResults.show()

.. .. _Index_1DBrukerReal:
.. .. figure:: _static/images/example_process_1dbruker_real.png
..     :width: 400
..     :alt: 1D NMR Spectrum (Bruker Format)
..     :align: center

..     1D NMR Spectrum Imported in Bruker Format

.. Here only the real part of the spectrum is displayed (dnpResults.plot(ws['proc'].real)). The imaginary part of the spectrum can be displayed by changing the second line to

.. .. code-block:: python

..     dnplab.dnpResults.plot(ws['proc'].imag)

.. To display the unprocessed raw FID:

.. .. code-block::

..     dnplab.dnpResults.figure()
..     dnplab.dnpResults.plot(ws['raw'].real)
..     dnplab.dnpResults.plt.xlabel('t2 [s]')
..     dnplab.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
..     dnplab.dnpResults.show()

.. .. _Index_1DFIDBrukerReal:
.. .. figure:: _static/images/example_FID_1dbruker_real.png
..     :width: 400
..     :alt: Raw FID (Bruker Format)
..     :align: center

..     1D FID from raw data (Bruker Format)


.. Determine T1 from an Inversion Recovery Experiment
.. ==================================================

.. In this example, the data from an inversion recovery experiment is analyzed to extract the longitudinal relaxation time T1. This example uses the example script: *example_process_IRbruker.py*.

.. Import DNPLab, load data, and create a workspace in the same manner as demonstrated above in the first example.

.. For the 2D dataset, add the align function to the processing:

.. .. code-block:: python

..     dnplab.dnpNMR.remove_offset(ws,{})
..     dnplab.dnpNMR.window(ws,{'linewidth' : 10})
..     dnplab.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
..     dnplab.dnpNMR.align(ws, {})
..     dnplab.dnpNMR.autophase(ws,{})

.. To plot the processed NMR spectra:

.. .. code-block:: python

..     dnplab.dnpResults.plot(ws['ft'].real)
..     dnplab.dnpResults.xlim([-30,50])
..     dnplab.dnpResults.plt.xlabel('Chemical Shift [ppm]')
..     dnplab.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
..     dnplab.dnpResults.figure()

.. .. _Index_IRBruker:
.. .. figure:: _static/images/example_process_IRbruker.png
..     :width: 400
..     :alt: Processed IR spectra
..     :align: center

..     Processed inversion recovery spectra (Bruker Format)

.. Next, the processed NMR spectra are copied to *ft* within the workspace, the signal amplitude for each NMR spectrum is integrated and the data is fitted to a function, describing inversion recovery polarization build-up.

.. .. code-block:: python

..     ws.copy('proc', 'ft')
..     dnplab.dnpNMR.integrate(ws, {'integrate_width' : 100, 'integrate_center' : 0})
..     dnplab.dnpFit.t1Fit(ws)

.. The T1 value can be displayed using:

.. .. code-block:: python

..     print('T1 value (sec) = ' + str(ws['fit'].attrs['t1']))
..     T1 value (sec) = 2.045498109768188


.. To plot the inversion-recovery build-up curve (experimental and fitted data):

.. .. code-block:: python

..     dnplab.dnpResults.plot(ws['proc'].real, 'o')
..     dnplab.dnpResults.plot(ws['fit'])
..     dnplab.dnpResults.show()

.. .. _Index_IRBuildUp:
.. .. figure:: _static/images/example_process_IRbuildup.png
..     :width: 400
..     :alt: Inversion Recovery Build-up
..     :align: center

..     Inversion recovery build-up (experimental and fit)


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This first example demonstrates how to import two NMR spectra, one recorded with a microwave po...">

.. only:: html

 .. figure:: /auto_examples/images/thumb/sphx_glr_plot_exampleDNPOnOff_thumb.png
     :alt: 01 - Microwave On/Off Signal

     :ref:`sphx_glr_auto_examples_plot_exampleDNPOnOff.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/plot_exampleDNPOnOff
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-gallery


  .. container:: sphx-glr-download sphx-glr-download-python

    :download:`Download all examples in Python source code: auto_examples_python.zip </auto_examples/auto_examples_python.zip>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

    :download:`Download all examples in Jupyter notebooks: auto_examples_jupyter.zip </auto_examples/auto_examples_jupyter.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
