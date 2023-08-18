from warnings import warn

import numpy as _np
from scipy.constants import *


def autophase(
    data,
    dim="f2",
    method="search",
    reference_range=None,
    pts_lim=None,
    order="zero",
    pivot=0,
    delta=0,
    phase=None,
    reference_slice=None,
    force_positive=False,
):
    """Automatically phase correct data, or apply manual phase correction

    Args:
        data (DNPData): Data object to autophase
        dim (str): Dimension to autophase
        method (str): Autophase method, "search" by default
        reference_range:
        pts_lim:
        order:
        pivot:
        delta:
        phase:
        reference_slice:
        force_positive:

    Returns:
        DNPData: Autophased data, including attrs "phase0" for order="zero", and "phase1" if order="first"

    .. math::

        \mathrm{data}         &= \exp(-1j * \mathrm{phase}) &

        \mathrm{phase(arctan)} &= \mathrm{arctan}(\mathrm{sum}(\mathrm{data.imag}) / \mathrm{sum}(\mathrm{data.real})) &

        \mathrm{phase(search)} &= \mathrm{argmax}(\mathrm{sum}(phased\_real^{2}) / \mathrm{sum}(phased\_imag^{2})) &

        phased\_real          &= \mathrm{data.real} * \exp(-1j * \mathrm{phase}) &

        phased\_imag          &= \mathrm{data.imag} * \exp(-1j * \mathrm{phase}) &

    """

    out = data.copy()
    if reference_slice == 0:
        raise ValueError(
            "please use indices from 1 to number of slices, i.e. use 1 instead of 0"
        )

    shape_data = _np.shape(out.values)
    index = out.dims.index(dim)

    if phase is not None:
        method = "manual"

    if method == "manual":
        if order == "zero" and isinstance(phase, (int, float)):
            out.attrs["phase0"] = phase
        elif order == "zero" and not isinstance(phase, (int, float)):
            raise ValueError(
                "for a zero order phase correction you must supply a single phase"
            )
        elif order == "first" and isinstance(phase, (int, float)):
            out.attrs["phase0"] = phase
            order = "zero"
            warn(
                "method=manual and order=first but only a single phase was given, switching to order=zero"
            )
        elif (
            order == "first"
            and isinstance(phase, (list, _np.ndarray))
            and len(phase) == shape_data[index]
        ):
            out.attrs["phase1"] = _np.array(phase)
        else:
            raise ValueError(
                "Invalid combination of phase order and phase value(s). Supply float for zero order, array or list for first order"
            )
    else:
        if isinstance(reference_range, (list, tuple)) and len(reference_range) == 2:
            check_data = out[dim, (reference_range[0], reference_range[1])]
        else:
            check_data = data.copy()
            if reference_range is not None:
                warn("reference_range must be None or list/tuple length=2")

        if reference_slice is not None:
            if len(shape_data) == 1:
                reference_slice = None
                temp_data = check_data.values
                warn("ignoring reference_slice, this is 1D data")
            else:
                temp_data = check_data.values[:, reference_slice - 1]
        else:
            temp_data = check_data.values

        if method == "arctan":
            out.attrs["phase0"] = _np.arctan(
                _np.sum(_np.imag(temp_data.reshape(-1, 1)))
                / _np.sum(_np.real(temp_data.reshape(-1, 1)))
            )
        elif method == "search":
            if pts_lim is not None:
                if len(check_data.coords[dim]) > pts_lim:
                    phasing_x = _np.linspace(
                        min(check_data.coords[dim]),
                        max(check_data.coords[dim]),
                        int(pts_lim),
                    ).reshape(-1)
                    if len(check_data.dims) == 1:
                        temp_data = _np.interp(
                            phasing_x, check_data.coords[dim], check_data.values
                        ).reshape(-1)
                    else:
                        ind_dim = list(set(out.dims) - set([dim]))[0]
                        ind_shape = out.shape[out.index(ind_dim)]
                        temp_data = _np.array(
                            [
                                _np.interp(
                                    phasing_x,
                                    check_data.coords[dim],
                                    check_data[dim, :].values[:, indx],
                                ).reshape(-1)
                                for indx in range(ind_shape)
                            ]
                        ).reshape(pts_lim, ind_shape)
            phases_0 = _np.linspace(-pi / 2, pi / 2, 180).reshape(-1)
            rotated_data = (temp_data.reshape(-1, 1)) * _np.exp(-1j * phases_0)
            real_imag_ratio = (_np.real(rotated_data) ** 2).sum(axis=0) / (
                (_np.imag(rotated_data) ** 2).sum(axis=0)
            )
            out.attrs["phase0"] = phases_0[_np.argmax(real_imag_ratio)]
        else:
            raise TypeError("Invalid autophase method")

    if order == "zero":
        out.values *= _np.exp(-1j * out.attrs["phase0"])
    elif order == "first":
        if method == "manual":
            out.attrs["phase1"] = phase
        else:
            pivot_ratio = pivot / len(out.values)
            out.attrs["phase1"] = _np.linspace(
                out.attrs["phase0"] - delta * pivot_ratio,
                out.attrs["phase0"] + delta * (1 - pivot_ratio),
                len(out.values),
            )
        out.values.T.dot(_np.exp(-1j * out.attrs["phase1"]))

    else:
        raise TypeError("Invalid order or order & phase pair")

    if force_positive:
        if _np.sum(out.values) < 0:
            out.values *= -1

    proc_parameters = {
        "method": method,
        "reference_range": reference_range,
        "reference_slice": reference_slice,
        "force_positive": force_positive,
        "order": order,
        "pivot": pivot,
        "delta": delta,
    }
    proc_attr_name = "autophase"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def phase_cycle(data, dim, receiver_phase):
    """Apply phase cycle to data

    Args:
        all_data (dnpdata_collection, dnpdata): data to process
        dim (str): dimension to perform phase cycle
        receiver_phase (numpy.array, list): Receiver Phase 0 (x), 1 (y), 2 (-x), 3 (-y)

    Returns:
        dnpdata: data object after phase cycle applied
    """

    out = data.copy()

    if dim not in out.dims:
        raise ValueError("dim not in dims")

    coord = out.coords[dim]
    receiver_phase = _np.array(receiver_phase).ravel()

    proc_parameters = {"dim": dim, "receiver_phase": receiver_phase}

    receiver_phase = _np.tile(receiver_phase, int(coord.size / receiver_phase.size))

    index = out.dims.index(dim)

    reshape_size = [1 for k in out.dims]
    reshape_size[index] = len(out.coords[dim])

    out *= _np.exp(-1j * (pi / 2.0) * receiver_phase.reshape(reshape_size))

    proc_attr_name = "phasecycle"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def phase(data, dim="f2", p0=0.0, p1=0.0, pivot=None):
    """Apply phase correction to DNPData object

    Args:
        data (DNPData): Data object to phase
        dim (str): Dimension to phase, default is "f2"
        p0 (float, array): Zero order phase correction (degree, 0 - 360)
        p1 (float, array): First order phase correction (degree, 0 - 360)
        picot (float): Pivot point for first order phase correction

    Returns:
        data (DNPData): Phased data, including new attributes "p0", "p1", and "pivot"

    Examples:

        0th-order phase correction of 1D or 2D DNPData object. If the DNPData object has multiple 1D spectra the same phase p0 is applied to all spectra.

        >>> data = dnp.phase(data,p0)

        0th-order phase correction of all spectra of a 2D DNPData object using a (numpy) array p0 of phases:

        >>> p0 = np.array([15, 15, 5, -5, 0])
        >>> data = dnp.phase(data, p0)

    .. Note::
        A 2D DNPData object can either be phase using a single p0 (p1) value, or using an array of phases. When using an array, the size of the phase array has to be equal to the number of spectra to be phased.

    """

    p0 = _np.mod(p0, 360)
    p1 = _np.mod(p1, 360)

    p0 = _np.array(p0 * pi / 180.0)  # p0 in radians
    p1 = _np.array(p1 * pi / 180.0)  # p1 in radians

    out = data.copy()
    out.unfold(dim)
    coord = out.coords[dim]

    phase = _np.exp(
        1.0j
        * (
            p0.reshape(1, -1)
            + (p1.reshape(1, -1) * _np.arange(coord.size).reshape(-1, 1) / coord.size)
        )
    )

    out *= phase

    out.fold()

    proc_parameters = {
        "p0": p0,
        "p1": p1,
        "pivot": pivot,
    }

    proc_attr_name = "phase_correction"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out
