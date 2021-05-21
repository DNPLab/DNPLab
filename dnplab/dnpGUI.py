""" dnpGUI

A graphical user interface for using DNPLab

"""

import sys
import os

from PyQt5 import QtWidgets, QtGui, QtCore

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QSizePolicy,
    QWidget,
    QPushButton,
    QLineEdit,
    QSlider,
    QLabel,
    QCheckBox,
    QFileDialog,
    QComboBox,
    QSpinBox,
)
from PyQt5.QtCore import Qt

if hasattr(QtCore.Qt, "AA_EnableHighDpiScaling"):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, "AA_UseHighDpiPixmaps"):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np
from scipy.io import loadmat, savemat
import copy

sys.path.append("/Users/thomascasey/dnplab")
import dnplab


class dnpGUI(QMainWindow):
    def __init__(self):

        super().__init__()

        self.setWindowTitle("DNPLab GUI")
        self.setGeometry(10, 10, 1200, 700)
        self.setContentsMargins(0, 0, 0, 0)

        self.dataplt = PlotCanvas(self, width=9, height=6.2)
        self.dataplt.move(5, 40)

        self.loadfilebutton = QPushButton("Load File", self)
        self.loadfilebutton.setStyleSheet(
            "font : bold ; color : rgb(70, 129, 43) ; background-color : rgb(167, 169, 172)"
        )
        self.loadfilebutton.move(5, 5)
        self.loadfilebutton.resize(90, 30)
        self.loadfilebutton.clicked.connect(self.import_data_file)

        self.loaddirbutton = QPushButton("Load Directory", self)
        self.loaddirbutton.setStyleSheet(
            "font : bold ; color : rgb(70, 129, 43) ; background-color : rgb(167, 169, 172)"
        )
        self.loaddirbutton.move(100, 5)
        self.loaddirbutton.resize(110, 30)
        self.loaddirbutton.clicked.connect(self.import_data_dir)

        self.pathLabel = QLabel(self)
        self.pathLabel.setStyleSheet("font : bold 14px; color : rgb(0, 0, 0)")
        self.pathLabel.move(215, 13)
        self.pathLabel.resize(690, 20)
        self.pathLabel.setText("Data folder path")

        self.realCheckbox = QCheckBox(self)
        self.realCheckbox.setStyleSheet("font : bold 14px")
        self.realCheckbox.move(20, 670)
        self.realCheckbox.resize(100, 20)
        self.realCheckbox.setText("Real")
        self.realCheckbox.setChecked(True)

        self.imagCheckbox = QCheckBox(self)
        self.imagCheckbox.setStyleSheet("font : bold 14px")
        self.imagCheckbox.move(90, 670)
        self.imagCheckbox.resize(100, 20)
        self.imagCheckbox.setText("Imag")
        self.imagCheckbox.setChecked(True)

        self.wskeySelect = QComboBox(self)
        self.wskeySelect.move(760, 670)
        self.wskeySelect.resize(150, 25)

        self.rmoffsetCheckbox = QCheckBox(self)
        self.rmoffsetCheckbox.setStyleSheet("font : bold 14px")
        self.rmoffsetCheckbox.move(920, 40)
        self.rmoffsetCheckbox.resize(200, 20)
        self.rmoffsetCheckbox.setText("Remove Offset")
        self.rmoffsetCheckbox.setChecked(True)

        self.leftshiftLabel = QLabel(self)
        self.leftshiftLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.leftshiftLabel.move(920, 65)  # 123, 590
        self.leftshiftLabel.resize(490, 30)
        self.leftshiftLabel.setText("Left Shift")

        self.leftshiftSpinbox = QSpinBox(self)
        # self.leftshiftSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.leftshiftSpinbox.move(1000, 70)  # 123, 590
        self.leftshiftSpinbox.resize(100, 20)

        self.windowLabel = QLabel(self)
        self.windowLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.windowLabel.move(920, 95)  # 123, 590
        self.windowLabel.resize(490, 30)
        self.windowLabel.setText("Window")

        self.windowSelect = QComboBox(self)
        self.windowSelect.move(915, 125)
        self.windowSelect.resize(150, 25)
        self.windowSelect.addItems(
            [
                "exponential",
                "gaussian",
                "hamming",
                "hann",
                "lorentz_gauss",
                "sin2",
                "traf",
            ]
        )

        self.explwLabel = QLabel(self)
        self.explwLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.explwLabel.move(925, 150)  # 123, 590
        self.explwLabel.resize(290, 30)
        self.explwLabel.setText("lw (exp)")

        self.explwSpinbox = QSpinBox(self)
        # self.explwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.explwSpinbox.move(990, 155)  # 123, 590
        self.explwSpinbox.resize(75, 20)

        self.gaussLabel = QLabel(self)
        self.gaussLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gaussLabel.move(925, 180)  # 123, 590
        self.gaussLabel.resize(290, 30)
        self.gaussLabel.setText("gauss")

        self.gausslwLabel = QLabel(self)
        self.gausslwLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gausslwLabel.move(978, 180)  # 123, 590
        self.gausslwLabel.resize(290, 30)
        self.gausslwLabel.setText("lw")

        self.gausslwSpinbox = QSpinBox(self)
        # self.gausslwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.gausslwSpinbox.move(995, 185)  # 123, 590
        self.gausslwSpinbox.resize(70, 20)

        self.gaussmaxLabel = QLabel(self)
        self.gaussmaxLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gaussmaxLabel.move(1080, 180)  # 123, 590
        self.gaussmaxLabel.resize(290, 30)
        self.gaussmaxLabel.setText("max")

        self.gaussmaxSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.gaussmaxSpinbox.move(1110, 185)  # 123, 590
        self.gaussmaxSpinbox.resize(70, 20)

        self.zerofillLabel = QLabel(self)
        self.zerofillLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.zerofillLabel.move(920, 215)  # 123, 590
        self.zerofillLabel.resize(550, 30)
        self.zerofillLabel.setText("Zero-fill Factor")

        self.zerofillSpinbox = QSpinBox(self)
        # self.leftshiftSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.zerofillSpinbox.move(1035, 220)  # 123, 590
        self.zerofillSpinbox.resize(100, 20)
        self.zerofillSpinbox.setValue(2)

        self.phaseLabel = QLabel(self)
        self.phaseLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.phaseLabel.move(920, 245)  # 123, 590
        self.phaseLabel.resize(490, 30)
        self.phaseLabel.setText("Phase")

        self.optphaseCheckbox = QCheckBox(self)
        self.optphaseCheckbox.setStyleSheet("font : 14px")
        self.optphaseCheckbox.move(970, 250)
        self.optphaseCheckbox.resize(100, 20)
        self.optphaseCheckbox.setText("Optimize")

        self.zerothLabel = QLabel(self)
        self.zerothLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.zerothLabel.move(930, 270)  # 123, 590
        self.zerothLabel.resize(40, 30)
        self.zerothLabel.setText("0th")

        self.zerothSlider = QSlider(Qt.Horizontal, self)
        # self.zerothSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.zerothSlider.setGeometry(960, 275, 225, 20)

        self.firstLabel = QLabel(self)
        self.firstLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.firstLabel.move(930, 300)  # 123, 590
        self.firstLabel.resize(40, 30)
        self.firstLabel.setText("1st")

        self.firstSlider = QSlider(Qt.Horizontal, self)
        # self.firstSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.firstSlider.setGeometry(960, 305, 225, 20)

        self.pivotLabel = QLabel(self)
        self.pivotLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.pivotLabel.move(920, 330)  # 123, 590
        self.pivotLabel.resize(40, 30)
        self.pivotLabel.setText("Pivot")

        self.pivotSlider = QSlider(Qt.Horizontal, self)
        # self.pivotSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.pivotSlider.setGeometry(960, 340, 225, 20)

        self.alignCheckbox = QCheckBox(self)
        self.alignCheckbox.setStyleSheet("font : bold 14px")
        self.alignCheckbox.move(920, 370)
        self.alignCheckbox.resize(100, 20)
        self.alignCheckbox.setText("Align")

        self.baselineLabel = QLabel(self)
        self.baselineLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.baselineLabel.move(920, 390)  # 123, 590
        self.baselineLabel.resize(490, 30)
        self.baselineLabel.setText("Baseline")

        self.baselineSelect = QComboBox(self)
        self.baselineSelect.move(915, 420)
        self.baselineSelect.resize(150, 25)
        self.baselineSelect.addItems(["none", "polynomial", "exponential"])

        self.orderLabel = QLabel(self)
        self.orderLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.orderLabel.move(1070, 415)  # 123, 590
        self.orderLabel.resize(290, 30)
        self.orderLabel.setText("order")

        self.baselineSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.baselineSpinbox.move(1110, 420)  # 123, 590
        self.baselineSpinbox.resize(70, 20)

        self.integrateCheckbox = QCheckBox(self)
        self.integrateCheckbox.setStyleSheet("font : bold 14px")
        self.integrateCheckbox.move(920, 455)
        self.integrateCheckbox.resize(100, 20)
        self.integrateCheckbox.setText("Integrate")

        self.centerLabel = QLabel(self)
        self.centerLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.centerLabel.move(935, 475)  # 123, 590
        self.centerLabel.resize(290, 30)
        self.centerLabel.setText("center")

        self.centerSpinbox = QSpinBox(self)
        # self.gausslwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.centerSpinbox.move(980, 480)  # 123, 590
        self.centerSpinbox.resize(70, 20)

        self.widthLabel = QLabel(self)
        self.widthLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.widthLabel.move(1070, 475)  # 123, 590
        self.widthLabel.resize(290, 30)
        self.widthLabel.setText("width")

        self.widthSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.widthSpinbox.move(1110, 480)  # 123, 590
        self.widthSpinbox.resize(70, 20)

        self.enhancementCheckbox = QCheckBox(self)
        self.enhancementCheckbox.setStyleSheet("font : bold 14px")
        self.enhancementCheckbox.move(920, 520)
        self.enhancementCheckbox.resize(200, 20)
        self.enhancementCheckbox.setText("Enhancements")

        self.s2nCheckbox = QCheckBox(self)
        self.s2nCheckbox.setStyleSheet("font : bold 14px")
        self.s2nCheckbox.move(920, 545)
        self.s2nCheckbox.resize(200, 20)
        self.s2nCheckbox.setText("Signal-to-Noise")

        self.fitLabel = QLabel(self)
        self.fitLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.fitLabel.move(920, 580)  # 123, 590
        self.fitLabel.resize(490, 30)
        self.fitLabel.setText("Fitting")

        self.fitSelect = QComboBox(self)
        self.fitSelect.move(915, 610)
        self.fitSelect.resize(150, 25)
        self.fitSelect.addItems(
            ["none", "T1 inv recov", "buildup", "T2", "mono-exp", "bi-exp"]
        )

        self.sliceLabel = QLabel(self)
        self.sliceLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.sliceLabel.move(235, 663)  # 123, 590
        self.sliceLabel.resize(100, 30)
        self.sliceLabel.setText("Slice")

        self.sliceSlider = QSlider(Qt.Horizontal, self)
        # self.sliceSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.sliceSlider.setGeometry(275, 670, 330, 20)

        self.sliceSpinbox = QSpinBox(self)
        # self.sliceSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.sliceSpinbox.move(615, 670)  # 123, 590
        self.sliceSpinbox.resize(70, 20)
        self.sliceSpinbox.setValue(1)

        self.savebutton = QPushButton("Save", self)
        self.savebutton.setStyleSheet(
            "font : bold 16px ; color : rgb(70, 129, 43) ; background-color : rgb(167, 169, 172)"
        )
        self.savebutton.move(1115, 650)
        self.savebutton.resize(75, 40)

        self.gui = {"cur_data": {}}
        self.ws = dnplab.create_workspace()

        self.initUI()

    def initUI(self):

        self.gui["rm_offset"] = False
        self.show()

    def connect_widgets(self):

        self.realCheckbox.clicked.connect(self.plot_data)
        self.imagCheckbox.clicked.connect(self.plot_data)

        self.wskeySelect.currentIndexChanged.connect(self.ws_key_select)
        self.sliceSlider.valueChanged[int].connect(self.plot_data)
        self.sliceSpinbox.valueChanged.connect(self.plot_data)

        self.rmoffsetCheckbox.clicked.connect(self.process_data)

        self.leftshiftSpinbox.valueChanged.connect(self.process_data)

        self.windowSelect.currentIndexChanged.connect(self.process_data)
        self.explwSpinbox.valueChanged.connect(self.process_data)
        self.gausslwSpinbox.valueChanged.connect(self.process_data)
        self.gaussmaxSpinbox.valueChanged.connect(self.process_data)

        self.zerofillSpinbox.valueChanged.connect(self.process_data)

        self.optphaseCheckbox.clicked.connect(self.phase_data)
        self.zerothSlider.valueChanged[int].connect(self.zeroth_phase_data)
        self.firstSlider.valueChanged[int].connect(self.first_phase_data)
        self.pivotSlider.valueChanged[int].connect(self.pivot_phase_data)

        self.alignCheckbox.clicked.connect(self.align_data)

        self.baselineSelect.currentIndexChanged.connect(self.baseline_correct_data)
        self.baselineSpinbox.valueChanged.connect(self.baseline_correct_data)

        self.integrateCheckbox.clicked.connect(self.integrate_data)
        self.centerSpinbox.valueChanged.connect(self.integrate_data)
        self.widthSpinbox.valueChanged.connect(self.integrate_data)

        self.enhancementCheckbox.clicked.connect(self.calc_enhancement)
        self.s2nCheckbox.clicked.connect(self.calc_s2n)

        self.savebutton.clicked.connect(self.save_data)

    def import_data_file(self):

        dirname = QFileDialog.getOpenFileName(self)
        if dirname[0]:
            flname = dirname[0]
        else:
            return
        self.flname = os.path.normpath(flname)

        self.import_data()

    def import_data_dir(self):

        # dirname = QFileDialog.getExistingDirectory(self)
        dirname = "/Users/thomascasey/dnplab/data/topspin/5"
        if dirname:
            flname = dirname
        else:
            return
        self.flname = os.path.normpath(flname)

        self.import_data()

    def import_data(self):

        self.wskeySelect.clear()
        self.wskeySelect.addItems(["raw", "proc"])

        data = dnplab.dnpImport.load(self.flname)

        self.ws.add("raw", data)

        self.process_data()

    def ws_key_select(self):

        self.gui["cur_data"] = copy.deepcopy(self.ws[self.wskeySelect.currentText()])

        self.gui["cur_data"].rename(self.gui["cur_data"].dims[0], "x2")

        if self.gui["cur_data"].ndim == 2:
            self.gui["cur_data"].rename(self.gui["cur_data"].dims[1], "x1")

        self.plot_data()

    def process_data(self):

        self.ws.copy("raw", "proc")

        if self.rmoffsetCheckbox.isChecked():
            dnplab.dnpNMR.remove_offset(self.ws)

        if self.leftshiftSpinbox.value() != 0:
            dnplab.dnpNMR.left_shift(
                self.ws,
            )

        if self.windowSelect.currentText() == "exponential":
            dnplab.dnpNMR.window(
                self.ws, type="exponential", linewidth=self.explwSpinbox.value()
            )
        elif self.windowSelect.currentText() == "gaussian":
            dnplab.dnpNMR.window(self.ws, type="gaussian")
        elif self.windowSelect.currentText() == "hamming":
            dnplab.dnpNMR.window(self.ws, type="hamming")
        elif self.windowSelect.currentText() == "hann":
            dnplab.dnpNMR.window(self.ws, type="hann")
        elif self.windowSelect.currentText() == "lorentz_gauss":
            dnplab.dnpNMR.window(self.ws, type="lorentz_gauss")
        elif self.windowSelect.currentText() == "sin2":
            dnplab.dnpNMR.window(self.ws, type="sin2")
        elif self.windowSelect.currentText() == "traf":
            dnplab.dnpNMR.window(self.ws, type="traf")

        dnplab.dnpNMR.fourier_transform(
            self.ws, zero_fill_factor=self.zerofillSpinbox.value()
        )

        if self.optphaseCheckbox.isChecked():
            self.phase_data()

        if self.ws["proc"].ndim > 1 and self.alignCheckbox.isChecked():
            self.align_data()

        if self.baselineSelect.currentText() in ["polynomial", "exponential"]:
            self.baseline_correct_data()

        if self.integrateCheckbox.isChecked():
            self.integrate_data()

        if self.enhancementCheckbox.isChecked():
            self.calc_enhancement()

        if self.s2nCheckbox.isChecked():
            self.calc_s2n()

        self.ws.copy("proc", "basic_proc")

        self.wskeySelect.setCurrentText("proc")

        self.ws_key_select()

    def phase_data(self):
        dnplab.dnpNMR.autophase(self.ws)

    def align_data(self):
        dnplab.dnpNMR.align(self.ws)

    def baseline_correct_data(self):
        dnplab.dnpTools.baseline(self.ws)

    def integrate_data(self):
        dnplab.dnpNMR.integrate(self.ws)

        if "integrals" not in [
            self.wskeySelect.itemText(i) for i in range(self.wskeySelect.count())
        ]:
            self.wskeySelect.addItem("integrals")

    def calc_enhancement(self):
        dnplab.dnpNMR.calculate_enhancement(self.ws)

    def calc_s2n(self):
        dnplab.dnpTools.signal_to_noise(self.ws)

    def save_data(self):
        pass

    # --Plot Colors--#
    # dark_green = '#46812B'
    # light_green = '#67AE3E'
    # dark_grey = '#4D4D4F'
    # light_grey = '#A7A9AC'
    # orange = '#F37021'
    # ucsb navy = '#003660'
    # ucsb yellow = '#FEBC11'

    def plot_data(self):

        self.dataplt.axes.cla()

        if self.realCheckbox.isChecked():
            self.dataplt.axes.plot(
                self.gui["cur_data"].coords["x2"], self.gui["cur_data"].values.real
            )

        if self.imagCheckbox.isChecked():
            self.dataplt.axes.plot(
                self.gui["cur_data"].coords["x2"], self.gui["cur_data"].values.imag
            )

        """
        if self.integrateCheckbox.isChecked():
            self.dataplt.axes.plot(self.gui["cur_data"].coords[self.int_x_lower], self.int_line)
            self.dataplt.axes.plot(self.gui["cur_data"].coords[self.int_x_upper], self.int_line)
            self.dataplt.axes.plot(self.gui["cur_data"].coords[self.int_x_center], self.gui["cur_data"].real[self.int_x_center],marker="o",
                    linestyle="none")

        if self.pivot_position:
            self.dataplt.axes.plot(self.gui["cur_data"].coords[self.pivot_position], self.gui["cur_data"].real[self.pivot_position],marker="x",
                    linestyle="none")
        """

        self.dataplt.draw()


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=2, height=1, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


def main_func():
    app = QApplication(sys.argv)
    ex = dnpGUI()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main_func()
