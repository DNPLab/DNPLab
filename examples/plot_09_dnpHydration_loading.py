# %%
"""
.. _09_dnpHydration_loading:

===========================================
09 - dnpHydration importing and calculating
===========================================

This example demonstrates how to load ODNP results from an h5 file and calculate hydration parameters.

"""
# # %%
# # Import DNPLab and also pyplot for plotting,

# import dnplab as dnp
# import matplotlib.pyplot as plt

# # %%
# # Load the saved h5 file containing the inputs needed for dnpHydration,

# ws = dnp.dnpImport.load("../example_hydration_inputs.h5")

# # %%
# # Pass the loaded workspace to the dnpHydration module. If you return the results as is shown below, the results will be a dictionary. The workspace will now contain a 'hydration_results' dictionary.

# results = dnp.dnpHydration.hydration(ws)

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
