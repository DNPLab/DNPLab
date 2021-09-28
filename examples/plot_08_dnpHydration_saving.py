# %%
"""
.. _08_dnpHydration_saving:

========================================
08 - dnpHydration calculating and saving
========================================

This example demonstrates how to calculate hydration parameters using the dnpHydration module and save the results in h5 format.

"""
# %%
# Import DNPLab, numpy for creating arrays, and pyplot for plotting,

import dnplab as dnp
import numpy as np
import matplotlib.pyplot as plt

# # %%
# # Define the enhancement and T1 arrays along with their corresponding power measurements,

# TESTSET = {
#     "T1_array": np.array(
#         [
#             2.020153734009,
#             2.276836030132750,
#             2.3708172489377400,
#             2.4428968088189100,
#             2.5709096032675700,
#         ]
#     ),
#     "T1_powers": np.array(
#         [
#             0.000589495934876689,
#             0.024242327290569100,
#             0.054429505156431400,
#             0.0862844940360515,
#             0.11617812912435900,
#         ]
#     ),
#     "E_array": np.array(
#         [
#             0.57794113752189,
#             -0.4688718613022250,
#             -0.5464528159680670,
#             -1.0725090541762200,
#             -1.4141203961920700,
#             -1.695789643686440,
#             -1.771840068080760,
#             -1.8420812985152700,
#             -1.97571340381877,
#             -2.091405209753480,
#             -2.1860546327712800,
#             -2.280712535872610,
#             -2.4709892163826400,
#             -2.5184316153191200,
#             -2.556110148443770,
#             -2.576413132701720,
#             -2.675593912859120,
#             -2.8153300703866400,
#             -2.897475156648710,
#             -3.0042154567120800,
#             -3.087886507216510,
#         ]
#     ),
#     "E_powers": np.array(
#         [
#             0.0006454923080882520,
#             0.004277023425898170,
#             0.004719543572446050,
#             0.00909714298712173,
#             0.01344187403986090,
#             0.01896059941058610,
#             0.02101937603827090,
#             0.022335737104727900,
#             0.026029715703921800,
#             0.02917012237740640,
#             0.0338523245243911,
#             0.03820738749745440,
#             0.04733370907740660,
#             0.05269608016472140,
#             0.053790874615060400,
#             0.05697639350179900,
#             0.06435487925718170,
#             0.07909179437004270,
#             0.08958910066880800,
#             0.1051813598911370,
#             0.11617812912435900,
#         ]
#     ),
# }

# # %%
# # Next, create a 'hydration_inputs' dictionary using your ODNP data. See dnpHydration documentation for detailed descriptions of the inputs,

# data = {
#             "E_array": TESTSET["E_array"],
#             "E_powers": TESTSET["E_powers"],
#             "T1_array": TESTSET["T1_array"],
#             "T1_powers": TESTSET["T1_powers"],
#             "field": 348.5,
#             "spin_C": 100,
#             "T10": 2.0,
#             "T100": 2.5,
#             "smax_model": "free",
#             "interpolate_method": "linear"
#         }
# # %%
# # If you would like to redefine the constants used in the calculations you may do so with a constants dictionary. This is optional, if this is not defined by you the defaults will be automatically used,

# constants = {
#             "ksigma_bulk": 95.4,
#             "krho_bulk": 353.4,
#             "klow_bulk": 366,
#             "tcorr_bulk": 54,
#             "D_H2O": 2.3e-9,
#             "D_SL": 4.1e-10,
#             "delta_T1_water": False,
#             "T1_water": False,
#             "macro_C": False,
#         }

# # %%
# # Create a workspace with the hydration inputs, or add to an existing workspace. In this case a new workspace is created,

# ws = dnp.create_workspace("hydration_inputs", data)

# # %%
# # If you redefined the constants make sure to add them to the workspace.

# ws.add("hydration_constants", constants)

# # %%
# # Pass the workspace to the dnpHydration module. If you return the results as is shown below, the results will be a dictionary. The workspace will now contain a 'hydration_results' dictionary.

# results = dnp.dnpHydration.hydration(ws)

# # %%
# # Save the inputs and results workspace in h5 format,

# dnp.dnpSave.save(ws, filename="../dnpHydration_example.h5")

# # %%
# # Plot your results if you wish. You can either use the 'hydration_results' dictionary from the workspace,

# plt.figure()
# plt.plot(ws["hydration_inputs"]["E_powers"], ws["hydration_results"]["ksigma_array"], label="Data")
# plt.plot(ws["hydration_inputs"]["E_powers"], ws["hydration_results"]["ksigma_fit"], label="Fit")
# plt.xlabel("Microwave Power")
# plt.ylabel("ksigma")
# plt.legend()
# plt.grid(True)
# plt.show()

# # %%
# # Or you may use the returned results dictionary,

# plt.figure()
# plt.plot(ws["hydration_inputs"]["E_powers"], results["ksigma_array"], label="Data")
# plt.plot(ws["hydration_inputs"]["E_powers"], results["ksigma_fit"], label="Fit")
# plt.xlabel("Microwave Power")
# plt.ylabel("ksigma")
# plt.legend()
# plt.grid(True)
# plt.show()
