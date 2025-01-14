# %%
"""
.. _plot_06_fitting:

=============================
Using the DNPLab Fit function
=============================

This example demonstrates how to use the DNPLab fit function on a dnpdata object.

"""
# %%
# The following example shows how the autophase function can be conveniently used.
# To get started, first, setup the python environment:

import dnplab as dnp
import numpy as np

# %%
# Let's generate a test data set with some noise
# the test data is a lorentzian distribution with some noise added

pts = 1024
x = np.r_[-50:50:1j*pts]

np.random.seed(101)

values = dnp.math.lineshape.lorentzian(x, 0, 0.5, 1.5)
values += np.random.randn(pts)*0.05

data = dnp.DNPData(values, ['f2'], [x])


# %%
# Now we guess the initial parameters of the fit
# we create a spectrum with the initial guess to compare to our test data

init_guess = [0, 0.5, 1.0]

guess_values = dnp.math.lineshape.lorentzian(x, *init_guess)

guess = dnp.DNPData(guess_values, ['f2'], [x])

# %%
# now we perform the fit
# the output is a dictionary of DNPData objects containting the "fit" and optimal parameters "popt"

out = dnp.fit(dnp.math.lineshape.lorentzian, data, 'f2', init_guess)
fit = out['fit']
popt = out['popt']

print('Optimal Fit Values')
print(popt.values) # print optimal fitting values


# %%
# Now we plot the data, initial guess and fit

figure('data')
dnp.plot(data, label = 'data')
dnp.plot(guess, label = 'guess')
dnp.plot(fit, label = 'fit')
legend()
show()
