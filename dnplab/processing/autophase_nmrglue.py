import nmrglue as ng
import dnplab as dnp
import numpy as np

import matplotlib.pyplot as plt

if __name__ == "__main__":

    # %%
    # Let's load some example data
    # the data consists of 4 fid that are phase cycled (0-90-180-270)

    data = dnp.load("../../data/prospa/water_phase_cycled/data.2d")

    # %%
    # data can be autophased, but beware that the data is not necessarily good to be autophased

    data = dnp.fourier_transform(data)
    print(data.dims,data.dtype,data.size,data.shape,len(data.coords['f2']),len(data['f2',:]),np.array(data['f2',:]).shape )

    # do some tests for derivatives
    values=data.values[:,0]

    phased=ng.proc_autophase.autops(values,'acme')

    plt.figure()
    data_real = np.real(phased)
    plt.plot(data_real,label='autophased')
    plt.plot(values,label='orig')
    plt.gca().set_title('nmrglue')

    plt.show()
