"""This module contains Matplotlib widgets for manually processing data. Examples include manually phasing and manually aligning NMR spectra
"""
import numpy as _np

from . import return_data, dnpdata, dnpdata_collection

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


def manual_align(all_data, dim):
    """Manually align spectra"""

    data, isDict = return_data(all_data)

    coord = data.coords[dim]
    max_index = int(data.size / (coord.size ** 2.0))

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)
    init_index = 0
    delta_index = 1
    l = plt.plot(data.coords["f2"], _np.real(data.values))
    ax.margins(x=0)

    axcolor = "lightgoldenrodyellow"
    axindex = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

    sindex = Slider(
        axindex,
        "index",
        -1 * max_index,
        max_index,
        valinit=init_index,
        valstep=delta_index,
    )

    def update(val):
        index = sindex.val
        ix = 0
        for line in l:
            line.set_ydata(_np.roll(data[dim, ix].values.ravel(), index * ix))
            ix += 1
        fig.canvas.draw_idle()

    sindex.on_changed(update)

    reset_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
    reset_button = Button(reset_ax, "Reset", color=axcolor, hovercolor="0.975")

    inc_ax = plt.axes([0.6, 0.025, 0.1, 0.04])
    inc_button = Button(inc_ax, "+", color=axcolor, hovercolor="0.975")

    dec_ax = plt.axes([0.4, 0.025, 0.1, 0.04])
    dec_button = Button(dec_ax, "-", color=axcolor, hovercolor="0.975")

    def reset(event):
        sindex.reset()

    def inc(event):
        sindex.set_val(sindex.val + 1)

    def dec(event):
        sindex.set_val(sindex.val - 1)

    reset_button.on_clicked(reset)
    inc_button.on_clicked(inc)
    dec_button.on_clicked(dec)

    plt.show()
    manual_index = sindex.val

    for ix, x in enumerate(data.coords[dim]):
        data[dim, ix] = _np.roll(data[dim, ix].values, manual_index * ix)
        ix += 1

    proc_parameters = {
        "dim": dim,
    }
    proc_attr_name = "manualalign"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def manualphase(all_data):
    """Manually Phase NMR Spectra"""

    data, isDict = return_data(all_data)

    fig, ax = plt.subplots()

    plt.title("Manual Phase")
    plt.subplots_adjust(bottom=0.25)
    init_phase = 0
    delta_phase = 0.1
    (l0,) = plt.plot(data.coords["f2"], _np.real(data.values), alpha=0.5)
    (l,) = plt.plot(data.coords["f2"], _np.real(data.values))
    ax.margins(x=0)

    axcolor = "lightgoldenrodyellow"
    axphase = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)

    sphase = Slider(
        axphase, "phase (deg)", -180, 180, valinit=init_phase, valstep=delta_phase
    )

    def update(val):
        phase = sphase.val
        y_data = _np.real(_np.exp(-1j * _np.pi * phase / 180.0) * data.values)
        margin = 0.1 * (_np.max(y_data) - _np.min(y_data))
        l.set_ydata(y_data)
        ax.set_ylim(_np.min(y_data) - margin, _np.max(y_data) + margin)
        fig.canvas.draw_idle()

    sphase.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, "Reset", color=axcolor, hovercolor="0.975")

    def reset(event):
        sphase.reset()

    button.on_clicked(reset)

    plt.show()
    manual_phase = sphase.val

    data *= _np.exp(-1j * _np.pi * manual_phase / 180.0)

    proc_attr_name = "manualphase"
    proc_parameters = {"phase": manual_phase}
    data.add_proc_attrs(proc_attr_name, proc_parameters)
    data.attrs["phase0"] = manual_phase

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data
