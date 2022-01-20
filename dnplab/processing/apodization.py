import numpy as np


def apodize(
    data,
    dim="t2",
    type="exponential",
    linewidth=1,
    gaussian_max=0,
    inverse=False,
):
    """Apply Apodization to data down given dimension

    .. math::

        \mathrm{exponential}    &=  \exp(-2t * \mathrm{linewidth}) &

        \mathrm{gaussian}       &=  \exp((\mathrm{linewidth[0]} * t) - (\mathrm{linewidth[1]} * t^{2})) &

        \mathrm{hamming}        &=  0.53836 + 0.46164\cos(\pi * n/(N-1)) &

        \mathrm{han}            &=  0.5 + 0.5\cos(\pi * n/(N-1)) &

        \mathrm{sin2}           &=  \cos((-0.5\pi * n/(N - 1)) + \pi)^{2} &

        \mathrm{lorentz\_gauss} &=  \exp(L -  G^{2}) &

               L(t)    &=  \pi * \mathrm{linewidth[0]} * t &

               G(t)    &=  0.6\pi * \mathrm{linewidth[1]} * (\mathrm{gaussian\_max} * (N - 1) - t) &

        \mathrm{traf}           &=  (f1 * (f1 + f2)) / (f1^{2} + f2^{2}) &

               f1(t)   &=  \exp(-t * \pi * \mathrm{linewidth[0]}) &

               f2(t)   &=  \exp((t - T) * \pi * \mathrm{linewidth[1]}) &


    Args:
        all_data (dnpdata, dict): data container

    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | parameter       | type                    | default       | description                                       |
    +=================+=========================+===============+===================================================+
    | dim             | str                     | 't2'          | Dimension to apply exponential apodization        |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | type            | str                     | 'exponential' | type of apodization                               |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | linewidth       | float, list, or ndarray | 1             | linewidths  in Hz                                 |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | gaussian_max    | float                   | 0             | Location of gaussian component maximum            |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | inverse         | boolean                 | False         | invert the window function                        |
    +-----------------+-------------------------+---------------+---------------------------------------------------+

    Returns:
        dnpdata: data object with window function applied, including attr "window"
    """

    dim_size = data.coords[dim].size
    shape_data = np.shape(data.values)
    index = data.index(dim)

    if (isinstance(linewidth, np.ndarray) or isinstance(linewidth, list)) and len(
        linewidth
    ) == 2:
        exp_lw = linewidth[0]
        gauss_lw = linewidth[1]
    elif isinstance(linewidth, (int, float)):
        exp_lw = linewidth
        gauss_lw = linewidth
    else:
        raise ValueError("linewidth must be int/float, or list/ndarray with len==2")

    if type == "exponential":
        apwin = dnpMath.exponential_window(all_data, dim, linewidth)
    elif type == "gaussian":
        apwin = dnpMath.gaussian_window(all_data, dim, [exp_lw, gauss_lw])
    elif type == "hamming":
        apwin = dnpMath.hamming_window(dim_size)
    elif type == "hann":
        apwin = dnpMath.hann_window(dim_size)
    elif type == "lorentz_gauss":
        apwin = dnpMath.lorentz_gauss_window(
            all_data, dim, exp_lw, gauss_lw, gaussian_max=gaussian_max
        )
    elif type == "sin2":
        apwin = dnpMath.sin2_window(dim_size)
    elif type == "traf":
        apwin = dnpMath.traf_window(all_data, dim, exp_lw)
    else:
        raise ValueError("Invalid window type")

    apwin.reshape(dim_size)

    if inverse:
        apwin = 1 / apwin

    new_shape = [1 if ix != index else shape_data[index] for ix in range(data.ndim)]
    apwin = np.reshape(apwin, new_shape)

    data.values *= apwin

    proc_parameters = {
        "type": type,
        "linewidth": linewidth,
        "dim": dim,
        "gaussian_max": gaussian_max,
        "inverse": inverse,
    }
    proc_attr_name = "window"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data