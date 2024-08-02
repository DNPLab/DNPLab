import numpy as _np

import matplotlib.pyplot as _plt
from matplotlib.widgets import Slider, Button, RadioButtons


def align_widget(data, dim):
    """Manually align spectra"""

    coord = data.coords[dim]
    max_index = int(data.size / (coord.size**2.0))

    fig, ax = _plt.subplots()
    _plt.subplots_adjust(left=0.25, bottom=0.25)
    init_index = 0
    delta_index = 1
    l = _plt.plot(data.coords["f2"], _np.real(data.values))
    ax.margins(x=0)

    axcolor = "lightgoldenrodyellow"
    axindex = _plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

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

    reset_ax = _plt.axes([0.8, 0.025, 0.1, 0.04])
    reset_button = Button(reset_ax, "Reset", color=axcolor, hovercolor="0.975")

    inc_ax = _plt.axes([0.6, 0.025, 0.1, 0.04])
    inc_button = Button(inc_ax, "+", color=axcolor, hovercolor="0.975")

    dec_ax = _plt.axes([0.4, 0.025, 0.1, 0.04])
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

    _plt.show()
    manual_index = sindex.val

    for ix, x in enumerate(data.coords[dim]):
        data[dim, ix] = _np.roll(data[dim, ix].values, manual_index * ix)
        ix += 1

    proc_parameters = {
        "dim": dim,
    }
    proc_attr_name = "manualalign"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data
