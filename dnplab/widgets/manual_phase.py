import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


def manual_phase(data, dim="f2"):
    """Manually Phase NMR Spectra"""

    fig, ax = plt.subplots()

    plt.title("Manual Phase")
    plt.subplots_adjust(bottom=0.25)
    init_phase = 0
    delta_phase = 1.0
    coord = data.coords[dim]
    size = coord.size
    values = data.values
    values = values.reshape(size, -1)
    values = values[:, 0].reshape(-1, 1)

    l0 = plt.plot(coord, values, alpha=0.5)
    l = plt.plot(coord, values)
    ax.margins(x=0)

    axcolor = "lightgoldenrodyellow"
    axphase = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)

    sphase = Slider(
        axphase, "phase (deg)", -180, 180, valinit=init_phase, valstep=delta_phase
    )

    def update(val):
        phase = sphase.val
        margin = 0
        min_y_data = 0
        max_y_data = 0
        for ix, line in enumerate(l):
            y_data = np.real(np.exp(-1j * np.pi * phase / 180.0) * values[:, ix])
            margin = max(0.1 * (np.max(y_data) - np.min(y_data)), margin)
            line.set_ydata(y_data)
            min_y_data = min(min_y_data, np.min(y_data))
            max_y_data = max(max_y_data, np.max(y_data))
        ax.set_ylim(min_y_data - margin, max_y_data + margin)
        fig.canvas.draw_idle()

    sphase.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    reset_button = Button(resetax, "Reset", color=axcolor, hovercolor="0.975")

    plus_90 = plt.axes([0.6, 0.025, 0.1, 0.04])
    plus_90_button = Button(plus_90, "+90", color=axcolor, hovercolor="0.975")

    minus_90 = plt.axes([0.4, 0.025, 0.1, 0.04])
    minus_90_button = Button(minus_90, "-90", color=axcolor, hovercolor="0.975")

    def reset(event):
        sphase.reset()

    def plus_90_phase(event):
        new_phase = ((sphase.val + 270) % 360) - 180
        sphase.set_val(new_phase)

    def minus_90_phase(event):
        new_phase = ((sphase.val + 90) % 360) - 180
        sphase.set_val(new_phase)

    reset_button.on_clicked(reset)
    plus_90_button.on_clicked(plus_90_phase)
    minus_90_button.on_clicked(minus_90_phase)

    plt.show()
    manual_phase = sphase.val

    data *= np.exp(-1j * np.pi * manual_phase / 180.0)

    proc_attr_name = "manualphase"
    proc_parameters = {"phase": manual_phase}
    data.add_proc_attrs(proc_attr_name, proc_parameters)
    data.attrs["phase0"] = manual_phase

    return data
