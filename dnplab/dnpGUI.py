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

        self.loadbutton = QPushButton("Load", self)
        self.loadbutton.setStyleSheet(
            "font : bold ; color : rgb(70, 129, 43) ; background-color : rgb(167, 169, 172)"
        )
        self.loadbutton.move(5, 5)
        self.loadbutton.resize(65, 30)
        self.loadbutton.clicked.connect(self.import_data)

        self.pathLabel = QLabel(self)
        self.pathLabel.setStyleSheet("font : bold 14px; color : rgb(0, 0, 0)")
        self.pathLabel.move(75, 13)
        self.pathLabel.resize(825, 20)
        self.pathLabel.setText("Data folder path")

        self.realCheckbox = QCheckBox(self)
        self.realCheckbox.setStyleSheet("font : bold 14px")
        self.realCheckbox.move(20, 670)
        self.realCheckbox.resize(100, 20)
        self.realCheckbox.setText("Real")

        self.imagCheckbox = QCheckBox(self)
        self.imagCheckbox.setStyleSheet("font : bold 14px")
        self.imagCheckbox.move(90, 670)
        self.imagCheckbox.resize(100, 20)
        self.imagCheckbox.setText("Imag")

        self.wskeySelect = QComboBox(self)
        self.wskeySelect.move(760, 670)
        self.wskeySelect.resize(150, 25)

        self.rmoffsetCheckbox = QCheckBox(self)
        self.rmoffsetCheckbox.setStyleSheet("font : bold 14px")
        self.rmoffsetCheckbox.move(920, 40)
        self.rmoffsetCheckbox.resize(200, 20)
        self.rmoffsetCheckbox.setText("Remove Offset")

        self.alignCheckbox = QCheckBox(self)
        self.alignCheckbox.setStyleSheet("font : bold 14px")
        self.alignCheckbox.move(1070, 40)
        self.alignCheckbox.resize(100, 20)
        self.alignCheckbox.setText("Align")

        self.leftshiftLabel = QLabel(self)
        self.leftshiftLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.leftshiftLabel.move(920, 70)  # 123, 590
        self.leftshiftLabel.resize(490, 30)
        self.leftshiftLabel.setText("Left Shift")

        self.leftshiftSpinbox = QSpinBox(self)
        # self.leftshiftSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.leftshiftSpinbox.move(1000, 75)  # 123, 590
        self.leftshiftSpinbox.resize(100, 20)

        self.windowLabel = QLabel(self)
        self.windowLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.windowLabel.move(920, 105)  # 123, 590
        self.windowLabel.resize(490, 30)
        self.windowLabel.setText("Window")

        self.windowSelect = QComboBox(self)
        self.windowSelect.move(915, 135)
        self.windowSelect.resize(150, 25)

        self.explwLabel = QLabel(self)
        self.explwLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.explwLabel.move(925, 160)  # 123, 590
        self.explwLabel.resize(290, 30)
        self.explwLabel.setText("lw (exp)")

        self.explwSpinbox = QSpinBox(self)
        # self.explwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.explwSpinbox.move(990, 165)  # 123, 590
        self.explwSpinbox.resize(75, 20)

        self.gaussLabel = QLabel(self)
        self.gaussLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gaussLabel.move(925, 190)  # 123, 590
        self.gaussLabel.resize(290, 30)
        self.gaussLabel.setText("gauss")

        self.gausslwLabel = QLabel(self)
        self.gausslwLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gausslwLabel.move(978, 190)  # 123, 590
        self.gausslwLabel.resize(290, 30)
        self.gausslwLabel.setText("lw")

        self.gausslwSpinbox = QSpinBox(self)
        # self.gausslwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.gausslwSpinbox.move(995, 195)  # 123, 590
        self.gausslwSpinbox.resize(70, 20)

        self.gaussmaxLabel = QLabel(self)
        self.gaussmaxLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.gaussmaxLabel.move(1080, 190)  # 123, 590
        self.gaussmaxLabel.resize(290, 30)
        self.gaussmaxLabel.setText("max")

        self.gaussmaxSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.gaussmaxSpinbox.move(1110, 195)  # 123, 590
        self.gaussmaxSpinbox.resize(70, 20)

        self.phaseLabel = QLabel(self)
        self.phaseLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.phaseLabel.move(920, 230)  # 123, 590
        self.phaseLabel.resize(490, 30)
        self.phaseLabel.setText("Phase")

        self.optphaseCheckbox = QCheckBox(self)
        self.optphaseCheckbox.setStyleSheet("font : 14px")
        self.optphaseCheckbox.move(970, 235)
        self.optphaseCheckbox.resize(100, 20)
        self.optphaseCheckbox.setText("Optimize")

        self.zerothLabel = QLabel(self)
        self.zerothLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.zerothLabel.move(930, 255)  # 123, 590
        self.zerothLabel.resize(40, 30)
        self.zerothLabel.setText("0th")

        self.zerothSlider = QSlider(Qt.Horizontal, self)
        # self.zerothSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.zerothSlider.setGeometry(960, 260, 225, 20)

        self.firstLabel = QLabel(self)
        self.firstLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.firstLabel.move(930, 285)  # 123, 590
        self.firstLabel.resize(40, 30)
        self.firstLabel.setText("1st")

        self.firstSlider = QSlider(Qt.Horizontal, self)
        # self.firstSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.firstSlider.setGeometry(960, 290, 225, 20)

        self.pivotLabel = QLabel(self)
        self.pivotLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.pivotLabel.move(920, 315)  # 123, 590
        self.pivotLabel.resize(40, 30)
        self.pivotLabel.setText("Pivot")

        self.pivotSlider = QSlider(Qt.Horizontal, self)
        # self.pivotSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.pivotSlider.setGeometry(960, 325, 225, 20)

        self.baselineLabel = QLabel(self)
        self.baselineLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.baselineLabel.move(920, 350)  # 123, 590
        self.baselineLabel.resize(490, 30)
        self.baselineLabel.setText("Baseline")

        self.baselineSelect = QComboBox(self)
        self.baselineSelect.move(915, 380)
        self.baselineSelect.resize(150, 25)

        self.orderLabel = QLabel(self)
        self.orderLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.orderLabel.move(1070, 375)  # 123, 590
        self.orderLabel.resize(290, 30)
        self.orderLabel.setText("order")

        self.baselineSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.baselineSpinbox.move(1110, 380)  # 123, 590
        self.baselineSpinbox.resize(70, 20)

        self.integrateCheckbox = QCheckBox(self)
        self.integrateCheckbox.setStyleSheet("font : bold 14px")
        self.integrateCheckbox.move(920, 415)
        self.integrateCheckbox.resize(100, 20)
        self.integrateCheckbox.setText("Integrate")

        self.centerLabel = QLabel(self)
        self.centerLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.centerLabel.move(935, 435)  # 123, 590
        self.centerLabel.resize(290, 30)
        self.centerLabel.setText("center")

        self.centerSpinbox = QSpinBox(self)
        # self.gausslwSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.centerSpinbox.move(980, 440)  # 123, 590
        self.centerSpinbox.resize(70, 20)

        self.widthLabel = QLabel(self)
        self.widthLabel.setStyleSheet(
            "font : 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.widthLabel.move(1070, 435)  # 123, 590
        self.widthLabel.resize(290, 30)
        self.widthLabel.setText("width")

        self.widthSpinbox = QSpinBox(self)
        # self.gaussmaxSpinbox.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.widthSpinbox.move(1110, 440)  # 123, 590
        self.widthSpinbox.resize(70, 20)

        self.enhancementCheckbox = QCheckBox(self)
        self.enhancementCheckbox.setStyleSheet("font : bold 14px")
        self.enhancementCheckbox.move(920, 480)
        self.enhancementCheckbox.resize(200, 20)
        self.enhancementCheckbox.setText("Enhancements")

        self.s2nCheckbox = QCheckBox(self)
        self.s2nCheckbox.setStyleSheet("font : bold 14px")
        self.s2nCheckbox.move(920, 505)
        self.s2nCheckbox.resize(200, 20)
        self.s2nCheckbox.setText("Signal-to-Noise")

        self.fitLabel = QLabel(self)
        self.fitLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.fitLabel.move(920, 540)  # 123, 590
        self.fitLabel.resize(490, 30)
        self.fitLabel.setText("Fitting")

        self.fitLabel = QComboBox(self)
        self.fitLabel.move(915, 570)
        self.fitLabel.resize(150, 25)

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

        self.savebutton = QPushButton("Save", self)
        self.savebutton.setStyleSheet(
            "font : bold 16px ; color : rgb(70, 129, 43) ; background-color : rgb(167, 169, 172)"
        )
        self.savebutton.move(1115, 650)
        self.savebutton.resize(75, 40)

        self.initUI()

    def initUI(self):
        self.show()

    def import_data(self):
        pass

    # --Plot Colors--#
    # dark_green = '#46812B'
    # light_green = '#67AE3E'
    # dark_grey = '#4D4D4F'
    # light_grey = '#A7A9AC'
    # orange = '#F37021'
    # ucsb navy = '#003660'
    # ucsb yellow = '#FEBC11'

    def dataplt(self):
        self.dataplt.axes.cla()


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
