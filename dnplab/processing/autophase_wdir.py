# bsd3 license
# from nmrglue
import numpy as np
import numpy as _np
import scipy.optimize
import math

# from dnplab
import phase

# nonsymetric stencils and nonuniform stencils would be possible but might only come at a later point
# https://en.wikipedia.org/wiki/Finite_difference_coefficient
cStencil = {
    1: [1 / 12, -2 / 3, 0, 2 / 3, -1 / 12],
    2: [-1 / 12, 4 / 3, -5 / 2, 4 / 3, -1 / 12],
    3: [1 / 8, -1, 13 / 8, 0, -13 / 8, 1, -1 / 8],
    4: [-1 / 6, 2, -13 / 2, 28 / 3, -13 / 2, 2, -1 / 6],
}
# forward and backward parts for endings
fStencil = {
    1: [-25 / 12, 4, -3, 4 / 3, -1 / 4],
    2: [15 / 4, -77 / 6, 107 / 6, -13, 61 / 12, -5 / 6],
    3: [-49 / 8, 29, -461 / 8, 62, -307 / 8, 13, -15 / 8],
    4: [28 / 3, -111 / 2, 142, -1219 / 6, 176, -185 / 2, 82 / 3, -7 / 2],
}
bStencil = {}
for key, value in fStencil.items():
    value = _np.array(value)
    if key % 2 == 0:
        bStencil[key] = np.array(list(reversed(value)))
    else:
        bStencil[key] = -1 * np.array(list(reversed(value)))


def deriveF(f, dx, deriv=4):
    df = np.zeros(f.size)  # f is 1d!

    global cStencil
    global fStencil
    global bStencil
    # central part
    try:
        cStencil_loc = cStencil[deriv]
        fStencil_loc = fStencil[deriv]
        bStencil_loc = bStencil[deriv]
    except KeyError:
        raise NotImplementedError(
            "only the following derivatives orders are implemented: {0}, but {1} was selected".format(
                cStencil.keys(), deriv
            )
        )

    n = len(cStencil_loc)

    startInd = n // 2
    stopInd = -(n // 2)
    for k in range(n):
        if k != n - 1:
            df[startInd:stopInd] = df[startInd:stopInd] + (
                cStencil_loc[k] * f[k : -(n - k - 1)]
            )
        else:
            df[startInd:stopInd] = df[startInd:stopInd] + (cStencil_loc[k] * f[k:])

    # forward & backward part
    fwLen = len(fStencil_loc)
    for k in range(fwLen):
        df[:startInd] = df[:startInd] + fStencil_loc[k] * f[k : startInd + k]

    bwLen = len(bStencil_loc)  # should be same as fwLen
    for k in range(bwLen):
        if k != (bwLen - 1):
            df[stopInd:] = (
                df[stopInd:]
                + bStencil_loc[k] * f[-(bwLen - k - 1) + stopInd : -(bwLen - k - 1)]
            )
        else:
            df[stopInd:] = df[stopInd:] + bStencil_loc[k] * f[stopInd:]

    return df / dx**deriv


def autophase(data, dim="f2", pivot=False, deriv=1, gamma=1e-5):
    """
    Phase correction according to
    An efficient algorithm for automatic phase correction of
    NMR spectra based on entropy minimization
    By Li Chen et Al
    JMR 158 (2002) 164-168

    - autophases all spectra along dim
    - if pivot is not False it is a ??? that defines the pivot element
    """
    if pivot == False:
        raw_data = np.array(data[dim, :])
        size = len(data.coords[dim])
        raw_data = raw_data.flat[
            :size
        ]  # now this is not necessarily true if the dimension is not the first ?
    else:
        # pivot is now a element along an axis, e.g. ('Average',5)
        raise NotImplementedError

    dx = _np.diff(data.coords[dim])[0]

    def optimfun(phaseList, data, deriv, gamma, dx, *args):
        ph0, ph1 = phaseList
        phasingFaktor = np.exp(
            1.0j * (ph0 + (ph1 * np.arange(data.size) / data.size))
        ).astype(data.dtype)

        # phase and take real part
        R = _np.real(data * phasingFaktor)
        dR = deriveF(R, dx, deriv=deriv)

        h = _np.abs(dR) / _np.sum(_np.abs(dR))
        h = h + _np.max(h) * 1e-12  # ensure positivity

        P = gamma * _np.sum((R - _np.abs(R)) ** 2 / 4)

        E = _np.sum(-_np.log(h) * h + P)

        return E

    # simple and easy for now

    # "simplex" method: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin.html#scipy.optimize.fmin

    # fopt, iterations, funcalls, warnflag, allvecs
    xopt, fopt, iter, funcalls, warnflags = scipy.optimize.fmin(
        optimfun,
        [0, 0],
        args=(raw_data, deriv, gamma, dx),
        disp=False,
        full_output=True,
    )
    if warnflags == 1:
        warnings.warn(
            "Optimization reached maximum number of function evaluations, solution might not be optimal!"
        )
    if warnflags == 2:
        warnings.warn(
            "Optimization reached maximum number of iterations, solution might not be optimal!"
        )
    ph0, ph1 = xopt

    data.add_proc_attrs("autophase", {"pivot": pivot, "deriv": deriv, "dim": dim})
    phasedSpectrum = phase.phase(
        data, dim=dim, p0=ph0 / 2 / math.pi * 360, p1=ph1 / 2 / math.pi * 360
    )
    """
    phasedSpectrum = (np.atleast_2d(
        np.exp(1.0j * (ph0 + (ph1 * np.arange(raw_data.size) / raw_data.size))).astype(
            data.dtype
        ) ).T
        * data
    )
    """
    return phasedSpectrum


import matplotlib.pyplot as plt

if __name__ == "__main__":
    import dnplab as dnp
    import nmrglue as ng

    # %%
    # Let's load some example data
    # the data consists of 4 fid that are phase cycled (0-90-180-270)

    data = dnp.load("../../data/prospa/water_phase_cycled/data.2d")

    # %%
    # data can be autophased, but beware that the data is not necessarily good to be autophased

    data = dnp.fourier_transform(data)
    print(
        data.dims,
        data.dtype,
        data.size,
        data.shape,
        len(data.coords["f2"]),
        len(data["f2", :]),
        np.array(data["f2", :]).shape,
    )

    # do some tests for derivatives
    n_avg = 1
    values = data.values[:, n_avg]

    data = data["Average", n_avg]

    phased = ng.proc_autophase.autops(values, "acme")
    phasedData = autophase(data, deriv=1)
    phasedData4 = autophase(data, deriv=4)

    data_real = np.real(phased)
    dnp4 = np.real(phasedData4).values[:, 0]
    dnp1 = np.real(phasedData).values[:, 0]

    print(phasedData4.proc_info())

    plt.figure()

    plt.plot(data_real, label="nmrglue")
    plt.plot(np.real(values), label="orig")
    plt.plot(np.real(phasedData).values[:, 0], linestyle="--", label="dnplab - order 1")
    plt.plot(
        np.real(phasedData4).values[:, 0], linestyle="--", label="dnplab - order 4"
    )
    plt.gca().set_title("nmrglue & dnplab")

    plt.legend()

    plt.figure()
    plt.plot(
        (data_real - dnp1) / dnp1,
        linestyle="--",
        label="(nmrglue - dnplab order 1)/dnplab order 1",
    )
    plt.plot(
        (data_real - dnp4) / dnp4,
        linestyle="--",
        label="(nmrglue - dnplab order 4)/dnplab order 1",
    )
    plt.plot((data_real - dnp1) / 1, label="(nmrglue - dnplab order 1)")
    plt.plot((data_real - dnp4) / 1, label="(nmrglue - dnplab order 4)")
    plt.gca().set_title("nmrglue & dnplab difference")

    plt.legend()

    plt.show()

    """
    # do some tests for derivatives

    x = np.linspace(-1, 1, 100)
    f1 = lambda x: x**3 + x**2 + x
    f1prime_1 = lambda x: 3 * x**2 + 2 * x + 1
    f1prime_4 = lambda x: np.zeros(x.size)
    f2 = lambda x: 2 * x**2 + 3 * x
    f2prime = lambda x: 4 * x + 3

    # plt.plot(x,f1(x))
    # plt.plot(x,f1prime(x),linestyle="dotted")
    dx = np.diff(x)[0]

    df1_4 = deriveF(f1(x), dx, deriv=4)
    df1_1 = deriveF(f1(x), dx, deriv=1)
    plt.plot(x, df1_4 - f1prime_4(x), linestyle="dotted", label="fourth deriv")
    plt.plot(x, df1_1 - f1prime_1(x), label="first deriv")
    plt.legend()

    plt.figure()
    plt.plot(x, deriveF(f2(x), dx, deriv=1) - f2prime(x))

    phasedData = autophase (data,deriv=1)
    phasedData4 = autophase (data,deriv=4)

    plt.figure()
    data_real = np.real(data)
    dnp.fancy_plot(data_real,title='orig')

    plt.figure()
    phasedReal = np.real(phasedData)
    dnp.fancy_plot(phasedReal,title='deriv = 1')

    plt.figure()
    phasedReal = np.real(phasedData4)
    dnp.fancy_plot(phasedReal,title='deriv = 4')

    """
    plt.show()
