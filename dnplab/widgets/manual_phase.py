import numpy as _np

from . import return_data, dnpdata, dnpdata_collection

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


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
