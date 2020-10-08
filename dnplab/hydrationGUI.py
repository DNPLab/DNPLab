""" hydrationGUI

A graphical user interface for using DNPLab to process Han Lab format ODNP data and calculating hydration parameters
using the dnpHydration module.

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
import re
import datetime
import time

import dnplab


class hydrationGUI(QMainWindow):
    def __init__(self):

        super().__init__()
        self.testmode = False  # set to True for testing, False for normal use
        self.testpath = ".."  # same as sys path to dnplab

        # self.setStyleSheet('background-color : #A7A9AC')

        self.setWindowTitle("ODNP Processing")
        self.setGeometry(10, 10, 1050, 625)
        self.setContentsMargins(0, 0, 0, 0)

        self.dataplt = PlotCanvas(self, width=7.2, height=4.8)
        self.dataplt.move(5, 40)

        self.enhplt = PlotCanvas(self, width=3.15, height=2)
        self.enhplt.move(730, 40)

        self.t1plt = PlotCanvas(self, width=3.15, height=2)
        self.t1plt.move(730, 260)

        self.guiresultButton = QPushButton("GUI Result", self)
        self.guiresultButton.setStyleSheet(
            "font : bold ; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.guiresultButton.move(5, 5)
        self.guiresultButton.resize(80, 30)
        self.guiresultButton.clicked.connect(self.GUI_Result_Button)

        self.workupButton = QPushButton("Workup", self)
        self.workupButton.setStyleSheet(
            "font : bold ; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.workupButton.move(90, 5)
        self.workupButton.resize(80, 30)
        self.workupButton.clicked.connect(self.Workup_Button)

        self.singleButton = QPushButton("Bruker", self)
        self.singleButton.setStyleSheet(
            "font : bold ; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.singleButton.move(175, 5)
        self.singleButton.resize(80, 30)
        self.singleButton.clicked.connect(self.Bruker_Button)

        self.hanlabButton = QPushButton("Han Lab", self)
        self.hanlabButton.setStyleSheet(
            "font : bold ; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.hanlabButton.move(260, 5)
        self.hanlabButton.resize(80, 30)
        self.hanlabButton.clicked.connect(self.Han_Lab_Button)

        self.pathLabel = QLabel(self)
        self.pathLabel.setStyleSheet("font : bold 14px; color : rgb(0, 54, 96)")
        self.pathLabel.move(345, 13)
        self.pathLabel.resize(700, 20)
        self.pathLabel.setText("Data folder path")

        self.phaseLabel = QLabel(self)
        self.phaseLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.phaseLabel.move(10, 520)  # 123, 590
        self.phaseLabel.resize(490, 30)
        self.phaseLabel.setText("   Phase Adjust:")

        self.phaseSlider = QSlider(Qt.Horizontal, self)
        # self.phaseSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.phaseSlider.setGeometry(120, 526, 365, 20)

        self.intcenterLabel = QLabel(self)
        self.intcenterLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.intcenterLabel.move(5, 551)
        self.intcenterLabel.resize(490, 30)
        self.intcenterLabel.setText(" Window center:")

        self.intcenterSlider = QSlider(Qt.Horizontal, self)
        # self.intcenterSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.intcenterSlider.setGeometry(120, 557, 365, 20)

        self.intwindowLabel = QLabel(self)
        self.intwindowLabel.setStyleSheet(
            "font : bold 14px"
        )  # ; color : rgb(254, 188, 17)') ; background-color : rgb(0, 54, 96)')
        self.intwindowLabel.move(6, 582)
        self.intwindowLabel.resize(490, 30)
        self.intwindowLabel.setText("  Window width:")

        self.intwindowSlider = QSlider(Qt.Horizontal, self)
        # self.intwindowSlider.setStyleSheet('background-color : rgb(0, 54, 96)')
        self.intwindowSlider.setGeometry(195, 588, 290, 20)

        self.intwindowEdit = QLineEdit(self)
        self.intwindowEdit.move(120, 588)
        self.intwindowEdit.resize(35, 25)
        self.intwindowEdit.setText("10")

        self.inteditLabel = QLabel(self)
        self.inteditLabel.setStyleSheet("font : bold 14px")
        self.inteditLabel.move(160, 585)
        self.inteditLabel.resize(50, 30)
        self.inteditLabel.setText("ppm")

        self.optphsCheckbox = QCheckBox(self)
        self.optphsCheckbox.setStyleSheet("font : bold 14px")
        self.optphsCheckbox.move(490, 526)
        self.optphsCheckbox.resize(100, 20)
        self.optphsCheckbox.setText("Optimize")

        self.optcentCheckbox = QCheckBox(self)
        self.optcentCheckbox.setStyleSheet("font : bold 14px")
        self.optcentCheckbox.move(490, 557)
        self.optcentCheckbox.resize(100, 20)
        self.optcentCheckbox.setText("Optimize")

        self.optwidthCheckbox = QCheckBox(self)
        self.optwidthCheckbox.setStyleSheet("font : bold 14px")
        self.optwidthCheckbox.move(490, 588)
        self.optwidthCheckbox.resize(100, 20)
        self.optwidthCheckbox.setText("Optimize")

        self.onlyT1pCheckbox = QCheckBox(self)
        self.onlyT1pCheckbox.setStyleSheet("font : bold 14px")
        self.onlyT1pCheckbox.move(727, 575)
        self.onlyT1pCheckbox.resize(100, 20)
        self.onlyT1pCheckbox.setText("Only T1(p)")
        self.onlyT1pCheckbox.setChecked(False)

        self.onlyT10Checkbox = QCheckBox(self)
        self.onlyT10Checkbox.setStyleSheet("font : bold 14px")
        self.onlyT10Checkbox.move(727, 595)
        self.onlyT10Checkbox.resize(100, 20)
        self.onlyT10Checkbox.setText("Only T1(0)")

        self.nextButton = QPushButton("Next", self)
        self.nextButton.setStyleSheet(
            "font : bold 14px; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.nextButton.move(625, 525)
        self.nextButton.resize(100, 40)

        self.autoButton = QPushButton("Auto Process", self)
        self.autoButton.setStyleSheet(
            "font : bold 14px; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.autoButton.move(740, 525)
        self.autoButton.resize(100, 40)

        self.backButton = QPushButton("Back", self)
        self.backButton.setStyleSheet(
            "font : bold 14px; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.backButton.move(625, 575)
        self.backButton.resize(100, 40)

        self.dnpLab_errorLabel = QLabel(self)
        self.dnpLab_errorLabel.setStyleSheet("font : bold 14px")
        self.dnpLab_errorLabel.move(615, 545)
        self.dnpLab_errorLabel.resize(500, 20)
        self.dnpLab_errorLabel.setText("DNPLab fit Error")

        self.workup_errorLabel = QLabel(self)
        self.workup_errorLabel.setStyleSheet("font : bold 14px")
        self.workup_errorLabel.move(615, 525)
        self.workup_errorLabel.resize(500, 20)
        self.workup_errorLabel.setText("Workup fit Error")

        self.t1fitLabel = QLabel(self)
        self.t1fitLabel.setStyleSheet("font : bold 14px")
        self.t1fitLabel.move(750, 470)
        self.t1fitLabel.resize(230, 20)
        self.t1fitLabel.setText("T1 interpolation:")

        self.linearfitCheckbox = QCheckBox(self)
        self.linearfitCheckbox.setStyleSheet("font : bold 14px")
        self.linearfitCheckbox.move(865, 470)
        self.linearfitCheckbox.resize(100, 20)
        self.linearfitCheckbox.setText("Linear")

        self.order2fitCheckbox = QCheckBox(self)
        self.order2fitCheckbox.setStyleSheet("font : bold 14px")
        self.order2fitCheckbox.move(930, 470)
        self.order2fitCheckbox.resize(100, 20)
        self.order2fitCheckbox.setText("2nd Order")

        self.exclude1T1Checkbox = QCheckBox(self)
        self.exclude1T1Checkbox.setStyleSheet("font : bold 14px")
        self.exclude1T1Checkbox.move(865, 500)
        self.exclude1T1Checkbox.resize(150, 20)
        self.exclude1T1Checkbox.setText("Exclude first T1(p)")

        self.t10Label = QLabel(self)
        self.t10Label.setStyleSheet("font : bold 14px")
        self.t10Label.move(73, 525)
        self.t10Label.resize(80, 20)
        self.t10Label.setText("T1(0) (s):")

        self.t10Edit = QLineEdit(self)
        self.t10Edit.move(140, 525)
        self.t10Edit.resize(65, 25)
        self.t10Edit.setText("2.5")

        self.workupt10Label = QLabel(self)
        self.workupt10Label.setStyleSheet("font : bold 14px")
        self.workupt10Label.move(420, 525)
        self.workupt10Label.resize(150, 20)
        self.workupt10Label.setText("workup T1(0) (s):")

        self.workupt10Edit = QLineEdit(self)
        self.workupt10Edit.move(545, 525)
        self.workupt10Edit.resize(65, 25)
        self.workupt10Edit.setText("2.5")

        self.show_wrkupCheckbox = QCheckBox(self)
        self.show_wrkupCheckbox.setStyleSheet("font : bold 14px")
        self.show_wrkupCheckbox.move(420, 550)
        self.show_wrkupCheckbox.resize(130, 20)
        self.show_wrkupCheckbox.setText("Show workup")

        self.fit_wrkupCheckbox = QCheckBox(self)
        self.fit_wrkupCheckbox.setStyleSheet("font : bold 14px")
        self.fit_wrkupCheckbox.move(420, 570)
        self.fit_wrkupCheckbox.resize(130, 20)
        self.fit_wrkupCheckbox.setText("Fit workup")

        self.t100Label = QLabel(self)
        self.t100Label.setStyleSheet("font : bold 14px")
        self.t100Label.move(68, 560)
        self.t100Label.resize(80, 20)
        self.t100Label.setText("T1<sub>0</sub>(0) (s):")

        self.t100Edit = QLineEdit(self)
        self.t100Edit.move(140, 560)
        self.t100Edit.resize(65, 25)
        self.t100Edit.setText("2.5")

        self.slcLabel = QLabel(self)
        self.slcLabel.setStyleSheet("font : bold 14px")
        self.slcLabel.move(43, 595)
        self.slcLabel.resize(150, 20)
        self.slcLabel.setText("Spin [C] (uM):")

        self.slcEdit = QLineEdit(self)
        self.slcEdit.move(140, 595)
        self.slcEdit.resize(65, 25)
        self.slcEdit.setText("1")

        self.fieldLabel = QLabel(self)
        self.fieldLabel.setStyleSheet("font : bold 14px")
        self.fieldLabel.move(225, 525)
        self.fieldLabel.resize(150, 20)
        self.fieldLabel.setText("Field (mT):")

        self.fieldEdit = QLineEdit(self)
        self.fieldEdit.move(305, 525)
        self.fieldEdit.resize(65, 25)
        self.fieldEdit.setText("348.5")

        self.smaxLabel = QLabel(self)
        self.smaxLabel.setStyleSheet("font : bold 14px")
        self.smaxLabel.move(225, 560)
        self.smaxLabel.resize(100, 20)
        self.smaxLabel.setText("s<sub>max</sub> model:")

        self.tetheredCheckbox = QCheckBox(self)
        self.tetheredCheckbox.setStyleSheet("font : bold 14px")
        self.tetheredCheckbox.move(305, 560)
        self.tetheredCheckbox.resize(100, 20)
        self.tetheredCheckbox.setText("Tethered")

        self.freeCheckbox = QCheckBox(self)
        self.freeCheckbox.setStyleSheet("font : bold 14px")
        self.freeCheckbox.move(305, 580)
        self.freeCheckbox.resize(100, 20)
        self.freeCheckbox.setText("Free")

        self.matoutButton = QPushButton("Save results", self)
        self.matoutButton.setStyleSheet(
            "font : bold 14px; color : rgb(254, 188, 17) ; background-color : rgb(0, 54, 96)"
        )
        self.matoutButton.move(925, 590)
        self.matoutButton.resize(120, 30)

        self.gui_dict = {
            "gui_function": {},
            "folder_structure": {},
            "rawdata_function": {},
            "processing_spec": {},
            "workup_function": {},
            "dnpLab_function": {},
            "workup_data": {},
            "dnpLab_data": {},
            "hydration_results": {},
            "data_plot": {},
            "enhancement_plot": {},
            "t1_plot": {},
            "t1_fit": {},
        }

        self.initUI()

    def initUI(self):

        self.gui_dict["gui_function"]["buttons"] = False
        self.gui_dict["gui_function"]["sliders"] = False
        self.gui_dict["gui_function"]["T100_process"] = False

        self.intwindowSlider.setMinimum(1)
        self.intwindowSlider.setMaximum(100)
        self.gui_dict["processing_spec"]["integration_width"] = 10
        self.intwindowSlider.setValue(
            self.gui_dict["processing_spec"]["integration_width"]
        )

        self.gui_dict["processing_spec"]["integration_center"] = 0
        self.intcenterSlider.setMinimum(
            self.gui_dict["processing_spec"]["integration_center"] - 50
        )
        self.intcenterSlider.setMaximum(
            self.gui_dict["processing_spec"]["integration_center"] + 50
        )
        self.intcenterSlider.setValue(
            self.gui_dict["processing_spec"]["integration_center"]
        )

        self.gui_dict["processing_spec"]["original_phase"] = 0

        # set blank plots
        self.reset_plots()
        self.connect_widgets()

        self.show()

    def show_hide_components(self):

        if self.gui_dict["gui_function"]["calculating"]:
            self.t100Label.setVisible(True)
            self.t100Edit.setVisible(True)
            self.t1fitLabel.setVisible(True)
            self.linearfitCheckbox.setVisible(True)
            self.order2fitCheckbox.setVisible(True)
            self.exclude1T1Checkbox.setVisible(True)
            self.slcLabel.setVisible(True)
            self.slcEdit.setVisible(True)
            self.fieldLabel.setVisible(True)
            self.fieldEdit.setVisible(True)
            self.smaxLabel.setVisible(True)
            self.tetheredCheckbox.setVisible(True)
            self.freeCheckbox.setVisible(True)
            self.matoutButton.setVisible(True)
            self.backButton.setVisible(True)
            self.onlyT1pCheckbox.setVisible(True)
            self.onlyT10Checkbox.setVisible(True)

            self.nextButton.setVisible(False)
            self.autoButton.setVisible(False)
            self.backButton.setText("Restart")
            self.intcenterLabel.setVisible(False)
            self.intcenterSlider.setVisible(False)
            self.intwindowLabel.setVisible(False)
            self.intwindowSlider.setVisible(False)
            self.intwindowEdit.setVisible(False)
            self.inteditLabel.setVisible(False)
            self.phaseLabel.setVisible(False)
            self.phaseSlider.setVisible(False)
            self.optcentCheckbox.setVisible(False)
            self.optphsCheckbox.setVisible(False)
            self.optwidthCheckbox.setVisible(False)

            if (
                self.gui_dict["gui_function"]["isWorkup"]
                or self.gui_dict["gui_function"]["isLab"]
            ):
                self.backButton.setVisible(False)
                self.onlyT1pCheckbox.setVisible(False)
                self.onlyT10Checkbox.setVisible(False)

            if self.gui_dict["gui_function"]["isWorkup"]:
                self.workupt10Label.setVisible(True)
                self.workupt10Edit.setVisible(True)
            else:
                self.t10Label.setVisible(True)
                self.t10Edit.setVisible(True)

            if self.gui_dict["gui_function"]["addWorkup"]:
                self.show_wrkupCheckbox.setVisible(True)
                self.fit_wrkupCheckbox.setVisible(True)
                self.show_wrkupCheckbox.setChecked(True)
                self.gui_dict["workup_function"]["show"] = True

                self.gui_dict["workup_function"]["fit"] = False
                self.fit_wrkupCheckbox.setChecked(False)

        else:
            self.intcenterLabel.setVisible(True)
            self.intcenterSlider.setVisible(True)
            self.intwindowLabel.setVisible(True)
            self.intwindowSlider.setVisible(True)
            self.intwindowEdit.setVisible(True)
            self.inteditLabel.setVisible(True)
            self.phaseLabel.setVisible(True)
            self.phaseSlider.setVisible(True)
            self.optcentCheckbox.setVisible(True)
            self.optphsCheckbox.setVisible(True)
            self.optwidthCheckbox.setVisible(True)
            self.autoButton.setVisible(True)
            self.backButton.setVisible(True)
            self.nextButton.setVisible(True)

            self.t10Label.setVisible(False)
            self.t10Edit.setVisible(False)
            self.workupt10Label.setVisible(False)
            self.workupt10Edit.setVisible(False)
            self.t100Label.setVisible(False)
            self.t100Edit.setVisible(False)
            self.t1fitLabel.setVisible(False)
            self.linearfitCheckbox.setVisible(False)
            self.order2fitCheckbox.setVisible(False)
            self.exclude1T1Checkbox.setVisible(False)
            self.slcLabel.setVisible(False)
            self.slcEdit.setVisible(False)
            self.fieldLabel.setVisible(False)
            self.fieldEdit.setVisible(False)
            self.smaxLabel.setVisible(False)
            self.tetheredCheckbox.setVisible(False)
            self.freeCheckbox.setVisible(False)
            self.matoutButton.setVisible(False)
            self.show_wrkupCheckbox.setVisible(False)
            self.fit_wrkupCheckbox.setVisible(False)
            self.backButton.setText("Back")
            self.onlyT1pCheckbox.setVisible(False)
            self.onlyT10Checkbox.setVisible(False)

        self.t1plt.setVisible(True)
        self.enhplt.setVisible(True)

        self.dnpLab_errorLabel.setVisible(False)
        self.workup_errorLabel.setVisible(False)

    def connect_widgets(self):

        self.intcenterSlider.valueChanged[int].connect(self.Integration_Center_Slider)
        self.intwindowSlider.valueChanged[int].connect(self.Integration_Window_Slider)
        self.intwindowEdit.editingFinished.connect(self.Integration_Window_Edit)
        self.phaseSlider.valueChanged[int].connect(self.Spectrum_Phase_Slider)
        self.optcentCheckbox.clicked.connect(self.Optimize_Center_Checkbox)
        self.optcentCheckbox.setChecked(True)
        self.optphsCheckbox.clicked.connect(self.Optimize_Phase_Checkbox)
        self.optphsCheckbox.setChecked(True)
        self.optwidthCheckbox.clicked.connect(self.Optimize_Width_Checkbox)
        self.optwidthCheckbox.setChecked(False)
        self.nextButton.clicked.connect(self.Next_Button)
        self.autoButton.clicked.connect(self.Auto_Process_Button)
        self.backButton.clicked.connect(self.Back_Button)
        self.onlyT1pCheckbox.clicked.connect(self.Only_T1p_Checkbox)
        self.onlyT1pCheckbox.setChecked(False)
        self.onlyT10Checkbox.clicked.connect(self.Only_T10_Checkbox)
        self.onlyT10Checkbox.setChecked(False)
        self.show_wrkupCheckbox.clicked.connect(self.Show_Workup_Checkbox)
        self.show_wrkupCheckbox.setChecked(True)
        self.fit_wrkupCheckbox.clicked.connect(self.Fit_Workup_Checkbox)
        self.fit_wrkupCheckbox.setChecked(False)
        self.t100Edit.editingFinished.connect(self.Edit_Hydration_Inputs)
        self.t10Edit.editingFinished.connect(self.Edit_Hydration_Inputs)
        self.workupt10Edit.editingFinished.connect(self.Edit_Hydration_Inputs)
        self.linearfitCheckbox.clicked.connect(self.Linear_Interpolation_Checkbox)
        self.linearfitCheckbox.setChecked(False)
        self.order2fitCheckbox.clicked.connect(self.SecondOrder_Interpolation_Checkbox)
        self.order2fitCheckbox.setChecked(True)
        self.exclude1T1Checkbox.clicked.connect(self.Exclude_FirstT1_Checkbox)
        self.exclude1T1Checkbox.setChecked(False)
        self.slcEdit.editingFinished.connect(self.Edit_Hydration_Inputs)
        self.fieldEdit.editingFinished.connect(self.Edit_Hydration_Inputs)
        self.tetheredCheckbox.clicked.connect(self.Smax_Tethered_Checkbox)
        self.tetheredCheckbox.setChecked(True)
        self.freeCheckbox.clicked.connect(self.Smax_Free_Checkbox)
        self.freeCheckbox.setChecked(False)
        self.matoutButton.clicked.connect(self.Save_Results_Button)

    def reset_plots(self):

        self.gui_dict["gui_function"]["autoProcess"] = False

        self.gui_dict["data_plot"]["xdata"] = []
        self.gui_dict["data_plot"]["ydata"] = []
        self.gui_dict["data_plot"]["xmin"] = -1
        self.gui_dict["data_plot"]["xmax"] = 1
        self.gui_dict["data_plot"]["plotksig"] = False
        self.gui_dict["data_plot"]["title"] = "Spectrum"
        self.plot_data()

        self.gui_dict["enhancement_plot"]["xdata"] = []
        self.gui_dict["enhancement_plot"]["ydata"] = []
        self.gui_dict["enhancement_plot"]["xmin"] = 0
        self.gui_dict["enhancement_plot"]["xmax"] = 1
        self.gui_dict["enhancement_plot"]["title"] = "E[p]"
        self.gui_dict["enhancement_plot"]["xLabel"] = "microwave power"
        self.gui_dict["enhancement_plot"]["yLabel"] = "enhancement"
        self.gui_dict["enhancement_plot"]["ytick"] = [0]
        self.gui_dict["enhancement_plot"]["ytickLabel"] = ["0"]
        self.gui_dict["enhancement_plot"]["tau"] = []
        self.gui_dict["enhancement_plot"]["t1Amps"] = []
        self.gui_dict["enhancement_plot"]["t1Fit"] = []
        self.gui_dict["enhancement_plot"]["plotT1fit"] = False
        self.gui_dict["enhancement_plot"]["plotEpfit"] = False
        self.plot_enh()

        self.gui_dict["t1_plot"]["xdata"] = []
        self.gui_dict["t1_plot"]["ydata"] = []
        self.gui_dict["t1_plot"]["xmin"] = 0
        self.gui_dict["t1_plot"]["xmax"] = 1
        self.gui_dict["t1_plot"]["ymin"] = 0
        self.gui_dict["t1_plot"]["ymax"] = 4
        self.gui_dict["t1_plot"]["title"] = r"$T_1[p]$"
        self.gui_dict["t1_plot"]["xLabel"] = "microwave power"
        self.gui_dict["t1_plot"]["yLabel"] = r"$T_1 (s)$"
        self.gui_dict["t1_plot"]["ytick"] = [1, 3]
        self.gui_dict["t1_plot"]["ytickLabel"] = ["1", "3"]
        self.gui_dict["t1_plot"]["plotT1interp"] = False
        self.plot_t1()

        self.gui_dict["gui_function"]["hydrationEdits"] = False
        self.gui_dict["gui_function"]["calculating"] = False

        self.gui_dict["folder_structure"]["index"] = 0
        self.gui_dict["processing_spec"]["phase_factor"] = 0

        self.show_hide_components()

    def plot_setter(self):

        if self.gui_dict["rawdata_function"]["folder"] == -1:
            self.gui_dict["data_plot"]["title"] = (
                r"$T_1$ Measurement, Folder # " + self.singlefolder
            )
            self.gui_dict["enhancement_plot"]["title"] = r"$T_1$ Fit"
            self.gui_dict["enhancement_plot"]["xLabel"] = r"$\tau$"
            self.gui_dict["enhancement_plot"]["yLabel"] = r"$M_z$"

        elif self.gui_dict["rawdata_function"]["folder"] == -2:
            self.gui_dict["data_plot"]["title"] = (
                "1D Data, Folder # " + self.singlefolder
            )
        else:

            if (
                self.gui_dict["rawdata_function"]["folder"]
                == self.gui_dict["folder_structure"]["p0"]
            ):
                self.gui_dict["data_plot"][
                    "title"
                ] = "Signal with power=0, Folder # " + str(
                    self.gui_dict["rawdata_function"]["folder"]
                )
                self.backButton.setText("Back")
                self.onlyT1pCheckbox.setVisible(False)
                self.onlyT10Checkbox.setVisible(False)

            elif (
                self.gui_dict["rawdata_function"]["folder"]
                == self.gui_dict["folder_structure"]["T10"]
            ):
                self.gui_dict["data_plot"][
                    "title"
                ] = "T1 with power=0, Folder # " + str(
                    self.gui_dict["rawdata_function"]["folder"]
                )

            elif (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["enh"]
            ):
                self.gui_dict["data_plot"][
                    "title"
                ] = "Enhanced Signal, Folder # " + str(
                    self.gui_dict["rawdata_function"]["folder"]
                )

            elif (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["T1"]
            ):
                self.gui_dict["data_plot"][
                    "title"
                ] = r"$T_1$ Measurement, Folder # " + str(
                    self.gui_dict["rawdata_function"]["folder"]
                )

            if (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["T1"]
                or self.gui_dict["rawdata_function"]["folder"]
                == self.gui_dict["folder_structure"]["T10"]
            ):
                self.gui_dict["enhancement_plot"]["title"] = r"$T_1$ Fit"
                self.gui_dict["enhancement_plot"]["xLabel"] = r"$\tau$"
                self.gui_dict["enhancement_plot"]["yLabel"] = r"$M_z$"
                self.gui_dict["enhancement_plot"]["ytick"] = [0]
                self.gui_dict["enhancement_plot"]["ytickLabel"] = ["0"]
                self.gui_dict["enhancement_plot"]["plotT1fit"] = True
                self.gui_dict["enhancement_plot"]["plotEpfit"] = False
                self.gui_dict["enhancement_plot"]["plotT1interp"] = False
            else:
                self.gui_dict["enhancement_plot"]["title"] = "E[p]"
                self.gui_dict["enhancement_plot"]["xLabel"] = "microwave power"
                self.gui_dict["enhancement_plot"]["yLabel"] = "enhancement"
                self.gui_dict["enhancement_plot"]["plotT1fit"] = False
                self.gui_dict["enhancement_plot"]["plotEpfit"] = False

            self.gui_dict["processing_spec"]["phase_factor"] = 0
            self.gui_dict["gui_function"]["sliders"] = False
            self.phaseSlider.setValue(0)
            self.gui_dict["gui_function"]["sliders"] = True

    def optCenter(self, width):

        optcenter_workspace = copy.deepcopy(self.processing_workspace)
        intgrl_array = []
        max_index = np.argmax(optcenter_workspace["proc"].values, axis=0)
        if max_index.size != 1:
            max_index = max_index[0]
        starting_center = round(optcenter_workspace["proc"].coords["t2"][max_index])
        indx = range(starting_center - 50, starting_center + 50)
        optcenter_workspace = self.phs_workspace(
            optcenter_workspace, self.gui_dict["processing_spec"]["original_phase"]
        )
        int_params = {"integrate_width": width}
        for k in indx:
            iterativeopt_workspace = copy.deepcopy(optcenter_workspace)
            int_params["integrate_center"] = k
            iterativeopt_workspace = self.int_workspace(
                iterativeopt_workspace, int_params
            )
            if len(iterativeopt_workspace["proc"].values) > 1:
                intgrl_array.append(sum(abs(iterativeopt_workspace["proc"].values)))
            else:
                intgrl_array.append(abs(iterativeopt_workspace["proc"].values[0]))
        cent = np.argmax(intgrl_array)
        self.gui_dict["processing_spec"]["integration_center"] = indx[cent]

    @staticmethod
    def import_create_workspace(dir, folder):

        data = dnplab.dnpImport.topspin.import_topspin(dir, folder)
        workspace = dnplab.create_workspace("raw", data)
        workspace.copy("raw", "proc")

        return workspace

    @staticmethod
    def int_workspace(workspace, int_params):

        dnplab.dnpNMR.integrate(
            workspace,
            integrate_center=int_params["integrate_center"],
            integrate_width=int_params["integrate_width"],
        )
        workspace["proc"].values = np.real(workspace["proc"].values)

        return workspace

    @staticmethod
    def phs_workspace(workspace, phase):

        workspace["proc"].values = workspace["proc"].values * np.exp(-1j * phase)

        return workspace

    def proc_T100(self):

        data = dnplab.dnpImport.topspin.import_topspin(
            self.gui_dict["rawdata_function"]["directory"], 100
        )
        self.dnpLab_workspace = dnplab.create_workspace("raw", data)
        self.dnpLab_workspace.copy("raw", "proc")

        self.processData()

        int_params = {
            "integrate_center": self.gui_dict["processing_spec"]["integration_center"],
            "integrate_width": self.gui_dict["processing_spec"]["integration_width"],
        }

        t100proc_workspace = self.phs_workspace(
            self.processing_workspace, self.gui_dict["processing_spec"]["phase"]
        )
        t100proc_workspace = self.int_workspace(t100proc_workspace, int_params)

        dnplab.dnpFit.t1Fit(t100proc_workspace)

        T100 = t100proc_workspace["fit"].attrs["t1"]
        self.t100Edit.setText(str(round(T100, 4)))

    def GUI_Result_Button(self):
        """Select either the h5 or the .mat files previously saved using the 'Save results' button."""
        try:
            if self.testmode:
                flname = os.path.join(
                    self.testpath,
                    "data",
                    "topspin",
                    "GUI_results hydrationGUI Results",
                    "GUI_results hydration_parameters.h5",
                )
                # flname = os.path.join(self.testpath, 'data', 'topspin', 'GUI_results hydrationGUI Results', 'GUI_results xODNP.mat')
            else:
                dirname = QFileDialog.getOpenFileName(self)

                if dirname[0]:
                    flname = dirname[0]
                else:
                    return

            print("GUI Results: " + flname)

            x = flname.split(os.sep)
            exten = flname.split(".")

            self.pathLabel.setText(
                "GUI RESULTS DIRECTORY: "
                + x[len(x) - 2]
                + " "
                + os.sep
                + " "
                + x[len(x) - 1]
            )

            self.ksiglabel = "DNPLab"
            self.gui_dict["rawdata_function"]["folder"] = -3
            self.gui_dict["gui_function"]["isLab"] = True
            self.gui_dict["gui_function"]["isWorkup"] = False
            self.gui_dict["gui_function"]["addWorkup"] = False
            self.gui_dict["workup_function"]["show"] = False
            self.gui_dict["workup_function"]["fit"] = False

            self.reset_plots()

            if "mat" in exten:
                matin = loadmat(flname)

                self.t10Edit.setText(str(round(float(matin["odnp"]["T10"]), 4)))

                self.gui_dict["dnpLab_data"]["T10"] = float(matin["odnp"]["T10"])
                self.gui_dict["dnpLab_data"]["T10_stdd"] = float(
                    matin["odnp"]["T10_stdd"]
                )
                epows = matin["odnp"]["Epowers"][0]
                self.gui_dict["dnpLab_data"]["Epowers"] = np.ravel(epows[0])
                ep = matin["odnp"]["Ep"][0]
                self.Ep = np.ravel(ep[0])
                t1pows = matin["odnp"]["T1powers"][0]
                self.gui_dict["dnpLab_data"]["T1powers"] = np.ravel(t1pows[0])
                t1p = matin["odnp"]["T1p"][0]
                self.T1p = np.ravel(t1p[0])
                t1perr = matin["odnp"]["T1p_stdd"][0]
                self.T1p_stdd = np.ravel(t1perr[0])

            elif "h5" in exten:
                h5in = dnplab.dnpImport.h5.loadh5(flname)

                self.t10Edit.setText(
                    str(round(float(h5in["hydration_inputs"]["T10"]), 4))
                )

                self.gui_dict["dnpLab_data"]["T10"] = float(
                    h5in["hydration_inputs"]["T10"]
                )
                self.gui_dict["dnpLab_data"]["T10_stdd"] = float(
                    h5in["hydration_results"]["T10_stdd"]
                )
                self.gui_dict["dnpLab_data"]["Epowers"] = h5in["hydration_inputs"][
                    "E_power"
                ]
                self.Ep = h5in["hydration_inputs"]["E"]
                self.gui_dict["dnpLab_data"]["T1powers"] = h5in["hydration_inputs"][
                    "T1_power"
                ]
                self.T1p = h5in["hydration_inputs"]["T1"]
                self.T1p_stdd = h5in["hydration_results"]["T1_stdd"]

            self.gui_dict["rawdata_function"]["nopowers"] = False

            self.finishProcessing()

            self.gui_dict["gui_function"]["buttons"] = True

        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText("File type error ")
            self.gui_dict["gui_function"]["buttons"] = False

    def Workup_Button(self):
        """Select the "Workup" folder that is the output of workup software used by the Han Lab."""
        try:
            if self.testmode:
                pthnm = os.path.join(self.testpath, "data", "topspin", "Workup")
            else:
                dirname = QFileDialog.getExistingDirectory(self)

                if dirname:
                    pthnm = dirname
                else:
                    return

            self.gui_dict["workup_function"]["directory"] = pthnm + os.sep
            print("Workup: " + pthnm)
            x = pthnm.split(os.sep)
            self.pathLabel.setText(
                "WORKUP DIRECTORY: "
                + x[len(x) - 3]
                + " "
                + os.sep
                + " "
                + x[len(x) - 2]
            )

            self.ksiglabel = "Workup"
            self.gui_dict["rawdata_function"]["folder"] = -3

            self.gui_dict["gui_function"]["isWorkup"] = True
            self.gui_dict["gui_function"]["isLab"] = False
            self.gui_dict["gui_function"]["addWorkup"] = False
            self.gui_dict["workup_function"]["show"] = True
            self.gui_dict["workup_function"]["fit"] = True

            self.reset_plots()
            self.processWorkup()

            self.t10Edit.setText(str(round(self.gui_dict["workup_data"]["T10"], 4)))

            self.gui_dict["rawdata_function"]["nopowers"] = False

            self.finishProcessing()

            self.gui_dict["gui_function"]["buttons"] = True

        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText("File type error")
            self.gui_dict["gui_function"]["buttons"] = False

    def processWorkup(self):

        Etest = np.loadtxt(
            self.gui_dict["workup_function"]["directory"] + "enhancementPowers.csv",
            delimiter=",",
            usecols=range(0, 3),
            max_rows=22,
            skiprows=1,
        )
        if Etest[0, 2] == 0:
            Eraw = Etest
        else:
            try:
                Eraw = np.loadtxt(
                    self.gui_dict["workup_function"]["directory"]
                    + "enhancementPowers.csv",
                    delimiter=",",
                    usecols=range(0, 2),
                    max_rows=30,
                    skiprows=1,
                )
            except IndexError:
                Eraw = np.loadtxt(
                    self.gui_dict["workup_function"]["directory"]
                    + "enhancementPowers.csv",
                    delimiter=",",
                    usecols=range(0, 2),
                    max_rows=22,
                    skiprows=1,
                )

        T1test = np.loadtxt(
            self.gui_dict["workup_function"]["directory"] + "t1Powers.csv",
            delimiter=",",
            usecols=range(0, 4),
            max_rows=6,
            skiprows=1,
        )
        if T1test[5, 3] == 304:
            T1raw = T1test
        else:
            T1test = np.loadtxt(
                self.gui_dict["workup_function"]["directory"] + "t1Powers.csv",
                delimiter=",",
                usecols=range(0, 4),
                max_rows=9,
                skiprows=1,
            )
            if T1test[8, 3] == 36:
                T1raw = T1test
            else:
                T1test = np.loadtxt(
                    self.gui_dict["workup_function"]["directory"] + "t1Powers.csv",
                    delimiter=",",
                    usecols=range(0, 4),
                    max_rows=10,
                    skiprows=1,
                )
                if T1test[9, 3] == 37:
                    T1raw = T1test
                else:
                    T1test = np.loadtxt(
                        self.gui_dict["workup_function"]["directory"] + "t1Powers.csv",
                        delimiter=",",
                        usecols=range(0, 4),
                        max_rows=11,
                        skiprows=1,
                    )
                    if T1test[10, 3] == 304:
                        T1raw = T1test

        ePows = Eraw[:, 0].reshape(-1)
        eP = Eraw[:, 1].reshape(-1)
        self.gui_dict["workup_data"]["Epowers"] = ePows[1 : len(ePows)]
        self.gui_dict["workup_data"]["Ep"] = eP[1 : len(ePows)]

        t1Pows = T1raw[:, 0].reshape(-1)
        t1P = T1raw[:, 1].reshape(-1)
        t1P_stdd = T1raw[:, 2].reshape(-1)
        self.gui_dict["workup_data"]["T1powers"] = t1Pows[0 : len(t1Pows) - 1]
        self.gui_dict["workup_data"]["T1p"] = t1P[0 : len(t1P) - 1]
        self.gui_dict["workup_data"]["T1p_stdd"] = t1P_stdd[0 : len(t1P_stdd) - 1]
        self.gui_dict["workup_data"]["T10"] = t1P[len(t1P) - 1]
        self.gui_dict["workup_data"]["T10_stdd"] = t1P_stdd[len(t1P_stdd) - 1]
        self.workupt10Edit.setText(str(round(self.gui_dict["workup_data"]["T10"], 4)))

        wrkupksig = np.loadtxt(
            self.gui_dict["workup_function"]["directory"] + "kSigma.csv",
            delimiter=",",
            usecols=range(0, 2),
            max_rows=1,
            skiprows=1,
        )

        self.gui_dict["workup_data"]["kSigma"] = wrkupksig[0] * 1e6
        self.gui_dict["workup_data"]["kSigma_stdd"] = wrkupksig[1] * 1e6

        # wrkupksig_array = np.loadtxt(self.gui_dict['workup_function']['directory'] + 'kSigma.csv', delimiter = ',', usecols = range(0,1), max_rows=21, skiprows = 6)
        # self.gui_dict['workup_data']['ksigma_powers'] = wrkupksig_array[:,1].reshape(-1)])
        # self.gui_dict['workup_data']['ksigma_array'] = wrkupksig_array[:,0].reshape(-1)])
        # ksig = np.transpose(ksig)
        # ksig = ksig[ksig[:,0].argsort()]
        # self.workup_ksig_array = ksig[:,1]

    @staticmethod
    def get_powers(fullPath, powerFile, bufferVal, exps):

        expTime = []
        absTime = []
        for exp in exps:
            opened = open(os.path.join(fullPath + str(exp), "audita.txt"))
            lines = opened.readlines()
            absStart = lines[8].split(" ")[2] + " " + lines[8].split(" ")[3]
            splitup = re.findall(r"[\w']+", absStart)
            absStart = datetime.datetime(
                int(splitup[0]),
                int(splitup[1]),
                int(splitup[2]),
                int(splitup[3]),
                int(splitup[4]),
                int(splitup[5]),
                int(splitup[6]),
            )
            absStart = time.mktime(absStart.utctimetuple())
            start = lines[8].split(" ")[3]
            start = start.split(":")
            hour = int(start[0], 10) * 3600
            minute = int(start[1], 10) * 60
            second = int(start[2].split(".")[0], 10)
            start = second + minute + hour
            absStop = lines[6].split("<")[1].split(">")[0].split(" ")
            absStop = absStop[0] + " " + absStop[1]
            splitup = re.findall(r"[\w']+", absStop)
            absStop = datetime.datetime(
                int(splitup[0]),
                int(splitup[1]),
                int(splitup[2]),
                int(splitup[3]),
                int(splitup[4]),
                int(splitup[5]),
                int(splitup[6]),
            )
            absStop = time.mktime(absStop.utctimetuple())
            stop = lines[6].split(" ")[4]
            stop = stop.split(":")
            hour = int(stop[0], 10) * 3600
            minute = int(stop[1], 10) * 60
            second = int(stop[2].split(".")[0], 10)
            stop = second + minute + hour
            expTime.append(stop - start)
            absTime.append((absStart, absStop))

        threshold = 20

        if os.path.isfile(fullPath + powerFile + ".mat"):
            print("Extracted powers from " + powerFile + ".mat file")
            openfile = loadmat(os.path.join(fullPath, powerFile + ".mat"))
            power = openfile.pop("powerlist")
            power = np.array([x for i in power for x in i])
            exptime = openfile.pop("timelist")
            exptime = np.array([x for i in exptime for x in i])
        elif os.path.isfile(fullPath + powerFile + ".csv"):
            print("Extracted powers from " + powerFile + ".csv file")
            openfile = open(os.path.join(fullPath, powerFile + ".csv", "r"))
            lines = openfile.readlines()
            if len(lines) == 1:
                lines = lines[0].split("\r")
            lines.pop(0)
            timeList = []
            powerList = []
            for line in lines:
                exptime, power = line.split("\r")[0].split(",")
                timeList.append(float(exptime))
                powerList.append(float(power))
            exptime = np.array(timeList)
            power = np.array(powerList)

        step = exptime[1] - exptime[0]
        dp = []
        for i in range(len(power) - 1):
            dp.append((power[i + 1] - power[i]) / step)
        dp = abs(np.array(dp))

        timeBreak = []
        for i in range(len(dp)):
            if dp[i] >= threshold:
                timeBreak.append(exptime[i])

        timeBreak.sort()
        absTime.sort(key=lambda tup: tup[0])
        offSet = absTime[-1][1] - timeBreak[-1] + bufferVal

        power_List = []
        for timeVals in absTime:
            start = int(timeVals[0] - offSet + bufferVal)
            stop = int(timeVals[1] - offSet - bufferVal)
            cutPower = []
            for k in range(0, len(exptime) - 1):
                if start <= exptime[k] <= stop:
                    cutPower.append(power[k])
            powers = round(np.average(cutPower), 3)
            power_List.append(float(powers))

        return power_List

    def Bruker_Button(self):
        """Select any numbered folder of a topspin dataset that contains 1D or 2D data."""
        try:
            if self.testmode:
                pthnm = os.path.join(self.testpath, "data", "topspin", "304")
            else:
                dirname = QFileDialog.getExistingDirectory(self)

                if dirname:
                    pthnm = dirname
                else:
                    return

            pthnm = pthnm + os.sep
            x = pthnm.split(os.sep)
            self.pathLabel.setText(
                "DATA DIRECTORY: " + x[len(x) - 3] + " " + os.sep + " " + x[len(x) - 2]
            )

            self.singlefolder = x[len(x) - 2]
            path = pthnm.replace(os.sep + str(self.singlefolder) + os.sep, os.sep)

            data = dnplab.dnpImport.topspin.import_topspin(path, self.singlefolder)
            self.dnpLab_workspace = dnplab.create_workspace("raw", data)
            self.dnpLab_workspace.copy("raw", "proc")

            if self.dnpLab_workspace["proc"].ndim == 2:
                print("T1 Measurement: " + pthnm)
                self.gui_dict["rawdata_function"]["folder"] = -1
            elif self.dnpLab_workspace["proc"].ndim == 1:
                print("1D Data: " + pthnm)
                self.gui_dict["rawdata_function"]["folder"] = -2

            self.reset_plots()
            self.plot_setter()

            self.gui_dict["gui_function"]["buttons"] = False
            self.gui_dict["gui_function"]["sliders"] = True
            self.gui_dict["gui_function"]["T100_process"] = False
            self.optcentCheckbox.setChecked(True)
            self.optphsCheckbox.setChecked(True)
            self.gui_dict["gui_function"]["isWorkup"] = False
            self.gui_dict["gui_function"]["addWorkup"] = False
            self.gui_dict["gui_function"]["isLab"] = False
            self.gui_dict["workup_function"]["fit"] = False
            self.gui_dict["workup_function"]["show"] = False
            self.gui_dict["enhancement_plot"]["plotT1fit"] = True
            self.backButton.setVisible(False)
            self.onlyT1pCheckbox.setVisible(False)
            self.onlyT10Checkbox.setVisible(False)
            self.nextButton.setVisible(False)
            self.autoButton.setVisible(False)
            self.t1plt.setVisible(False)

            self.processData()

            if self.gui_dict["rawdata_function"]["folder"] == -2:
                self.enhplt.setVisible(False)
            elif self.gui_dict["rawdata_function"]["folder"] == -1:
                self.enhplt.setVisible(True)
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText("Bruker data error")
            self.gui_dict["gui_function"]["sliders"] = False

    def Han_Lab_Button(self):
        """Select the base folder of a dataset generated using the 'rb_dnp1' command in topspin at UCSB.

        Required data:
            Folder 5: 1D spectrum that is collected without microwave power
            Folders 6-26: 1D spectra that are collected at different microwave powers specified in the power.mat file
            Folders 28-32: 2D inversion recovery experiments at different microwave powers specified in the t1powers.mat file
            Folder 304: 2D inversion recovery experiment collected without microwave power
            Additional required files: power.mat and t1powers.mat OR power.csv and t1power.csv files that are the measurements of applied microwave powers

        """
        try:
            if self.testmode:
                pthnm = os.path.join(self.testpath, "data", "topspin")
            else:
                dirname = QFileDialog.getExistingDirectory(self)

                if dirname:
                    pthnm = dirname
                else:
                    return

            pthnm = pthnm + os.sep
            self.gui_dict["rawdata_function"]["directory"] = pthnm
            print("Data: " + pthnm)
            x = pthnm.split(os.sep)
            self.pathLabel.setText(
                "DATA DIRECTORY: " + x[len(x) - 3] + " " + os.sep + " " + x[len(x) - 2]
            )

            self.gui_dict["folder_structure"] = {}
            self.gui_dict["gui_function"]["T100_process"] = False
            self.gui_dict["gui_function"]["isWorkup"] = False
            self.gui_dict["gui_function"]["isLab"] = False
            self.gui_dict["gui_function"]["addWorkup"] = False
            self.gui_dict["workup_function"]["show"] = False
            self.gui_dict["workup_function"]["fit"] = False
            self.nextButton.setText("Next")

            if os.path.exists(self.gui_dict["rawdata_function"]["directory"] + "40"):

                self.gui_dict["folder_structure"]["p0"] = 5
                self.gui_dict["folder_structure"]["enh"] = list(range(6, 30))
                self.gui_dict["folder_structure"]["T1"] = range(31, 41)
                self.gui_dict["folder_structure"]["T10"] = 304

            else:

                self.gui_dict["folder_structure"]["p0"] = 5
                self.gui_dict["folder_structure"]["enh"] = range(6, 27)
                self.gui_dict["folder_structure"]["T1"] = range(28, 33)
                self.gui_dict["folder_structure"]["T10"] = 304

            self.gui_dict["rawdata_function"]["nopowers"] = True

            if (
                os.path.exists(pthnm + "Workup" + os.sep)
                and os.path.isfile(
                    os.path.join(pthnm + "Workup", "enhancementPowers.csv")
                )
                and os.path.isfile(os.path.join(pthnm + "Workup", "kSigma.csv"))
                and os.path.isfile(os.path.join(pthnm + "Workup", "t1Powers.csv"))
            ):

                self.gui_dict["gui_function"]["addWorkup"] = True
                self.gui_dict["workup_function"]["show"] = True

                self.gui_dict["workup_function"]["directory"] = (
                    pthnm + "Workup" + os.sep
                )

                self.processWorkup()

                if len(self.gui_dict["workup_data"]["Epowers"]) == len(
                    self.gui_dict["folder_structure"]["enh"]
                ) and len(self.gui_dict["workup_data"]["T1powers"]) == len(
                    self.gui_dict["folder_structure"]["T1"]
                ):
                    Epowers = self.gui_dict["workup_data"]["Epowers"]
                    T1powers = self.gui_dict["workup_data"]["T1powers"]

                    print("Found Workup output, using power values from Workup.")
                    self.gui_dict["rawdata_function"]["nopowers"] = False

            if self.gui_dict["rawdata_function"]["nopowers"]:
                if os.path.isfile(os.path.join(pthnm, "power.mat")) or os.path.isfile(
                    os.path.join(pthnm, "power.csv")
                ):

                    if os.path.isfile(pthnm + "t1_powers.mat") or os.path.isfile(
                        pthnm + "t1_powers.csv"
                    ):
                        print("No Workup output found, using power readings files.")

                        E_power_List = self.get_powers(
                            self.gui_dict["rawdata_function"]["directory"],
                            "power",
                            2.5,
                            self.gui_dict["folder_structure"]["enh"],
                        )
                        # {{ These corrections to the power values are here to bring the powers to roughly the same magnitude as the results of the workup processing but should not be considered to be the actual correction. This can only be known by measuring the degree of attenuation difference between the path to the power meter and the path to the resonator
                        Epowers = np.add(E_power_List, 21.9992)
                        Epowers = np.divide(Epowers, 10)
                        Epowers = np.power(10, Epowers)
                        Epowers = np.multiply(1e-3, Epowers)
                        # }}

                        T1_power_List = self.get_powers(
                            self.gui_dict["rawdata_function"]["directory"],
                            "t1_powers",
                            20 * 2.5,
                            self.gui_dict["folder_structure"]["T1"],
                        )
                        # {{ These corrections to the power values are here to bring the powers to roughly the same magnitude as the results of the workup processing but should not be considered to be the actual correction. This can only be known by measuring the degree of attenuation difference between the path to the power meter and the path to the resonator
                        T1powers = np.add(T1_power_List, 21.9992)
                        T1powers = np.divide(T1powers, 10)
                        T1powers = np.power(10, T1powers)
                        T1powers = np.multiply(1e-3, T1powers)
                        # }}

                        self.gui_dict["rawdata_function"]["nopowers"] = False

            if self.gui_dict["rawdata_function"]["nopowers"]:
                print("No power readings found.")
                print("Trying to find power settings in experiment titles...")
                try:
                    Eplist = []
                    for k in self.gui_dict["folder_structure"]["enh"]:
                        title = dnplab.dnpImport.topspin.load_title(pthnm, expNum=k)
                        splitTitle = title.split(" ")
                        Eplist.append(float(splitTitle[-1]))

                    T1plist = []
                    for k in self.gui_dict["folder_structure"]["T1"]:
                        title = dnplab.dnpImport.topspin.load_title(pthnm, expNum=k)
                        splitTitle = title.split(" ")
                        T1plist.append(float(splitTitle[-1]))

                    # {{ These corrections to the power values are here to bring the powers to roughly the same magnitude as the results of the workup processing but should not be considered to be the actual correction. This can only be known by measuring the relationship between the attenuation setting, the power meter reading, and the power delivered to the resonator.
                    Epowers = np.multiply(-1, Eplist)
                    Epowers = np.add(Epowers, 29.01525)
                    Epowers = np.divide(Epowers, 10)
                    Epowers = np.power(10, Epowers)
                    Epowers = np.multiply(1e-3, Epowers)

                    T1powers = np.multiply(-1, T1plist)
                    T1powers = np.add(T1powers, 29.01525)
                    T1powers = np.divide(T1powers, 10)
                    T1powers = np.power(10, T1powers)
                    T1powers = np.multiply(1e-3, T1powers)
                    # }}

                    print(
                        "Powers taken from experiment titles. *WARNING: this is not accurate!"
                    )
                    self.gui_dict["rawdata_function"]["nopowers"] = False

                except:
                    print(
                        "No power readings available. E[p] and T1[p] are indexed by folder #. *WARNING: this is not accurate!"
                    )
                    Epowers = self.gui_dict["folder_structure"]["enh"]
                    T1powers = self.gui_dict["folder_structure"]["T1"]

            self.gui_dict["folder_structure"]["all"] = []
            self.gui_dict["folder_structure"]["all"].append(
                self.gui_dict["folder_structure"]["p0"]
            )
            for k in self.gui_dict["folder_structure"]["enh"]:
                self.gui_dict["folder_structure"]["all"].append(k)
            for k in self.gui_dict["folder_structure"]["T1"]:
                self.gui_dict["folder_structure"]["all"].append(k)
            self.gui_dict["folder_structure"]["all"].append(
                self.gui_dict["folder_structure"]["T10"]
            )

            self.Ep = []
            self.T1p = []
            self.T1p_stdd = []
            self.gui_dict["dnpLab_data"]["Epowers"] = Epowers
            self.gui_dict["dnpLab_data"]["T1powers"] = T1powers
            self.originalEPowers = self.gui_dict["dnpLab_data"]["Epowers"]
            self.originalT1Powers = self.gui_dict["dnpLab_data"]["T1powers"]
            self.gui_dict["gui_function"]["buttons"] = True
            self.gui_dict["gui_function"]["sliders"] = True
            self.optcentCheckbox.setChecked(True)
            self.optphsCheckbox.setChecked(True)

            self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                "folder_structure"
            ]["p0"]
            self.ksiglabel = "DNPLab"

            self.reset_plots()
            self.plot_setter()

            self.dnpLab_workspace = self.import_create_workspace(
                self.gui_dict["rawdata_function"]["directory"],
                self.gui_dict["rawdata_function"]["folder"],
            )

            self.processData()
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText("Han Lab data error ")
            self.gui_dict["gui_function"]["buttons"] = False
            self.gui_dict["gui_function"]["sliders"] = False

    def Next_Button(self):
        """Use the Next button to step through the data folders."""
        if self.gui_dict["gui_function"]["buttons"]:

            nextproc_workspace = self.phs_workspace(
                self.processing_workspace, self.gui_dict["processing_spec"]["phase"]
            )

            int_params = {
                "integrate_center": self.gui_dict["processing_spec"][
                    "integration_center"
                ],
                "integrate_width": self.gui_dict["processing_spec"][
                    "integration_width"
                ],
            }
            nextproc_workspace = self.int_workspace(nextproc_workspace, int_params)

            if (
                self.gui_dict["rawdata_function"]["folder"]
                == self.gui_dict["folder_structure"]["p0"]
            ):
                self.gui_dict["dnpLab_data"]["p0"] = nextproc_workspace["proc"].values[
                    0
                ]
            elif (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["enh"]
            ):
                Ep = (
                    nextproc_workspace["proc"].values[0]
                    / self.gui_dict["dnpLab_data"]["p0"]
                )
                self.Ep.append(np.real(Ep))
                if self.gui_dict["gui_function"]["autoProcess"]:
                    pass
                else:
                    self.gui_dict["enhancement_plot"]["xdata"] = self.gui_dict[
                        "dnpLab_data"
                    ]["Epowers"][0 : len(self.Ep)]
                    self.gui_dict["enhancement_plot"]["ydata"] = self.Ep
                    self.gui_dict["enhancement_plot"]["ytick"] = [0, min(self.Ep)]
                    if min(self.Ep) <= -10:
                        self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                            "0",
                            str(int(min(self.Ep))),
                        ]
                    else:
                        self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                            "0",
                            str(round(min(self.Ep), 1)),
                        ]
                    self.plot_enh()

            elif (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["T1"]
                or self.gui_dict["rawdata_function"]["folder"]
                == self.gui_dict["folder_structure"]["T10"]
            ):

                try:
                    dnplab.dnpFit.t1Fit(nextproc_workspace)

                    if (
                        self.gui_dict["rawdata_function"]["folder"]
                        in self.gui_dict["folder_structure"]["T1"]
                    ):
                        self.T1p.append(nextproc_workspace["fit"].attrs["t1"])
                        self.T1p_stdd.append(nextproc_workspace["fit"].attrs["t1_stdd"])
                    elif (
                        self.gui_dict["rawdata_function"]["folder"]
                        == self.gui_dict["folder_structure"]["T10"]
                    ):
                        self.gui_dict["dnpLab_data"]["T10"] = nextproc_workspace[
                            "fit"
                        ].attrs["t1"]
                        self.gui_dict["dnpLab_data"]["T10_stdd"] = nextproc_workspace[
                            "fit"
                        ].attrs["t1_stdd"]
                        self.t10Edit.setText(
                            str(round(self.gui_dict["dnpLab_data"]["T10"], 4))
                        )

                    if self.gui_dict["gui_function"]["autoProcess"]:
                        pass
                    else:

                        self.gui_dict["t1_fit"]["tau"] = np.reshape(
                            nextproc_workspace["proc"].coords, -1
                        )
                        self.gui_dict["t1_fit"]["t1Amps"] = nextproc_workspace[
                            "proc"
                        ].values

                        self.gui_dict["t1_fit"]["xaxis"] = nextproc_workspace[
                            "fit"
                        ].coords[0]
                        self.gui_dict["t1_fit"]["t1Fit"] = nextproc_workspace[
                            "fit"
                        ].values
                        self.gui_dict["t1_fit"]["t1Val"] = nextproc_workspace[
                            "fit"
                        ].attrs["t1"]

                        self.gui_dict["t1_plot"]["xdata"] = self.gui_dict[
                            "dnpLab_data"
                        ]["T1powers"][0 : len(self.T1p)]
                        self.gui_dict["t1_plot"]["ydata"] = self.T1p
                        self.gui_dict["t1_plot"]["ymin"] = (
                            min(self.gui_dict["t1_plot"]["ydata"]) * 0.9
                        )
                        self.gui_dict["t1_plot"]["ymax"] = (
                            max(self.gui_dict["t1_plot"]["ydata"]) * 1.1
                        )
                        self.gui_dict["t1_plot"]["ytick"] = [max(self.T1p)]
                        self.gui_dict["t1_plot"]["ytickLabel"] = [
                            str(round(max(self.T1p), 1))
                        ]

                        self.plot_t1()
                        self.plot_enh()

                except:
                    if (
                        self.gui_dict["folder_structure"]["all"][
                            self.gui_dict["folder_structure"]["index"]
                        ]
                        == 28
                    ):
                        print(
                            "WARNING: Error in first T1(p) fit, setting to ~0 and excluding from dnpHydration"
                        )
                        self.exclude1T1Checkbox.setChecked(True)
                        self.T1p.append(0.001)
                        self.T1p_stdd.append(0)
                    elif (
                        self.gui_dict["folder_structure"]["all"][
                            self.gui_dict["folder_structure"]["index"]
                        ]
                        == self.gui_dict["folder_structure"]["T10"]
                    ):
                        print(
                            "WARNING: Error in T1(0) fit, arbitrarily setting T1(0) = 2s"
                        )
                        self.gui_dict["dnpLab_data"]["T10"] = 2.0
                        self.gui_dict["dnpLab_data"]["T10_stdd"] = 0
                        self.t10Edit.setText(
                            str(round(self.gui_dict["dnpLab_data"]["T10"], 4))
                        )
                    else:
                        print(
                            "WARNING: Error in T1(p) fit for folder "
                            + str(
                                self.gui_dict["folder_structure"]["all"][
                                    self.gui_dict["folder_structure"]["index"]
                                ]
                            )
                            + ", setting equal to previous T1(p)"
                        )
                        self.T1p.append(self.T1p[-1])
                        self.T1p_stdd.append(0)

            self.gui_dict["folder_structure"]["index"] += 1
            if self.gui_dict["gui_function"]["autoProcess"]:
                print(
                    "Finished with Folder #"
                    + str(self.gui_dict["folder_structure"]["index"])
                    + " of "
                    + str(len(self.gui_dict["folder_structure"]["all"]))
                )

            if self.gui_dict["folder_structure"]["index"] >= len(
                self.gui_dict["folder_structure"]["all"]
            ):
                if os.path.exists(
                    self.gui_dict["rawdata_function"]["directory"] + "100" + os.sep
                ):
                    self.gui_dict["gui_function"]["T100_process"] = True
                    self.proc_T100()

                self.finishProcessing()

            else:
                self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                    "folder_structure"
                ]["all"][self.gui_dict["folder_structure"]["index"]]

                if self.gui_dict["gui_function"]["autoProcess"]:
                    pass
                else:
                    if (
                        self.gui_dict["folder_structure"]["index"]
                        == len(self.gui_dict["folder_structure"]["all"]) - 1
                    ):
                        self.nextButton.setText("Finish")
                    self.plot_setter()

                self.dnpLab_workspace = self.import_create_workspace(
                    self.gui_dict["rawdata_function"]["directory"],
                    self.gui_dict["rawdata_function"]["folder"],
                )

                self.processData()

        else:
            pass

    def Back_Button(self):
        """Use the Back button to return to the previous data folder."""
        if self.gui_dict["gui_function"]["buttons"]:

            self.gui_dict["folder_structure"]["index"] -= 1

            if self.gui_dict["folder_structure"]["index"] <= 0:
                self.gui_dict["folder_structure"]["index"] = 0
                self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                    "folder_structure"
                ]["p0"]

            if (
                self.gui_dict["folder_structure"]["index"]
                >= len(self.gui_dict["folder_structure"]["all"]) - 1
            ):
                self.reset_plots()
                self.gui_dict["dnpLab_data"]["Epowers"] = self.originalEPowers
                self.gui_dict["dnpLab_data"]["T1powers"] = self.originalT1Powers
                if self.onlyT10Checkbox.isChecked():
                    self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                        "folder_structure"
                    ]["T10"]
                    self.gui_dict["folder_structure"]["index"] = (
                        len(self.gui_dict["folder_structure"]["all"]) - 1
                    )
                    self.nextButton.setText("Finish")
                else:
                    self.nextButton.setText("Next")
                    if self.onlyT1pCheckbox.isChecked():
                        self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                            "folder_structure"
                        ]["T1"][0]
                        self.gui_dict["folder_structure"]["index"] = (
                            len(self.gui_dict["folder_structure"]["all"])
                            - 1
                            - len(self.gui_dict["folder_structure"]["T1"])
                        )
                        self.T1p = []
                        self.T1p_stdd = []
                    else:
                        self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                            "folder_structure"
                        ]["p0"]
                        self.Ep = []
                        self.T1p = []
                        self.T1p_stdd = []
            else:
                self.gui_dict["rawdata_function"]["folder"] = self.gui_dict[
                    "folder_structure"
                ]["all"][self.gui_dict["folder_structure"]["index"]]

            if (
                self.gui_dict["folder_structure"]["index"]
                == len(self.gui_dict["folder_structure"]["all"]) - 2
            ):
                self.nextButton.setText("Next")

            self.plot_setter()
            if (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["enh"]
            ):
                if len(self.Ep) < 2:
                    self.Ep = []
                    self.gui_dict["enhancement_plot"]["xdata"] = []
                    self.gui_dict["enhancement_plot"]["ydata"] = []
                else:
                    self.Ep = self.Ep[0 : len(self.Ep) - 1]
                    self.gui_dict["enhancement_plot"]["xdata"] = self.gui_dict[
                        "dnpLab_data"
                    ]["Epowers"][0 : len(self.Ep)]
                    self.gui_dict["enhancement_plot"]["ydata"] = self.Ep
                    self.gui_dict["enhancement_plot"]["ytick"] = [0, min(self.Ep)]
                    if min(self.Ep) <= -10:
                        self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                            "0",
                            str(int(min(self.Ep))),
                        ]
                    else:
                        self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                            "0",
                            str(round(min(self.Ep), 1)),
                        ]
                self.plot_enh()

            elif (
                self.gui_dict["rawdata_function"]["folder"]
                in self.gui_dict["folder_structure"]["T1"]
            ):
                if len(self.T1p) < 2:
                    self.T1p = []
                    self.T1p_stdd = []
                    self.gui_dict["t1_plot"]["xdata"] = []
                    self.gui_dict["t1_plot"]["ydata"] = []
                else:
                    self.T1p = self.T1p[0 : len(self.T1p) - 1]
                    self.T1p_stdd = self.T1p_stdd[0 : len(self.T1p_stdd) - 1]
                    self.gui_dict["t1_plot"]["xdata"] = self.gui_dict["dnpLab_data"][
                        "T1powers"
                    ][0 : len(self.T1p)]
                    self.gui_dict["t1_plot"]["ydata"] = self.T1p
                    self.gui_dict["t1_plot"]["ytick"] = [max(self.T1p)]
                    self.gui_dict["t1_plot"]["ytickLabel"] = [
                        str(round(max(self.T1p), 1))
                    ]
                    self.gui_dict["t1_plot"]["ymin"] = min(self.T1p) * 0.85
                    self.gui_dict["t1_plot"]["ymax"] = max(self.T1p) * 1.15
                self.plot_t1()

            self.dnpLab_workspace = self.import_create_workspace(
                self.gui_dict["rawdata_function"]["directory"],
                self.gui_dict["rawdata_function"]["folder"],
            )

            self.processData()

        else:
            pass

    def Auto_Process_Button(self):
        """Allow the correct phase and integration window to be automatically chosen and process the full ODNP dataset automatically."""
        if self.gui_dict["gui_function"]["buttons"]:
            try:
                self.optphsCheckbox.setChecked(True)
                self.optcentCheckbox.setChecked(True)
                print("Auto processing, please wait...")
                self.gui_dict["gui_function"]["autoProcess"] = True
                # t = time.time()
                for k in range(
                    self.gui_dict["folder_structure"]["index"] + 1,
                    len(self.gui_dict["folder_structure"]["all"]) + 1,
                ):
                    self.Next_Button()
                # elapsed = time.time() - t
                # print('AutoProcess Time = ' + str(elapsed))
            except:
                self.gui_dict["folder_structure"]["index"] = len(
                    self.gui_dict["folder_structure"]["all"]
                )
                self.Back_Button()
                print(
                    "Error in auto processing folder # "
                    + str(self.gui_dict["folder_structure"]["all"][k - 2])
                    + ", resetting to folder # "
                    + str(self.gui_dict["folder_structure"]["p0"])
                )
        else:
            pass

    def processData(self):

        self.processing_workspace = copy.deepcopy(self.dnpLab_workspace)
        dnplab.dnpNMR.remove_offset(self.processing_workspace)
        dnplab.dnpNMR.window(self.processing_workspace, linewidth=5)
        dnplab.dnpNMR.fourier_transform(self.processing_workspace, zero_fill_factor=2)

        if self.processing_workspace["proc"].ndim == 2:
            dnplab.dnpNMR.align(self.processing_workspace)

        if (
            self.optphsCheckbox.isChecked()
            or self.gui_dict["gui_function"]["autoProcess"]
        ):

            curve = self.processing_workspace["proc"].values

            self.gui_dict["processing_spec"]["original_phase"] = np.arctan(
                np.sum(np.imag(curve)) / np.sum(np.real(curve))
            )
            """
            phases = np.linspace(-np.pi / 2, np.pi / 2, 100).reshape(1, -1)
            rotated_data = (curve.reshape(-1, 1)) * np.exp(-1j * phases)
            success = (np.real(rotated_data) ** 2).sum(axis=0) / (
                (np.imag(rotated_data) ** 2).sum(axis=0)
            )
            bestindex = np.argmax(success)

            self.gui_dict["processing_spec"]["original_phase"] = phases[0, bestindex]
            """

        if (
            self.optcentCheckbox.isChecked()
            or self.gui_dict["gui_function"]["autoProcess"]
        ):
            self.optCenter(10)

        if self.optwidthCheckbox.isChecked():

            optwidth_workspace = copy.deepcopy(self.processing_workspace)
            optwidth_workspace = self.phs_workspace(
                optwidth_workspace, self.gui_dict["processing_spec"]["original_phase"]
            )

            ydata = abs(np.real(optwidth_workspace["proc"].values))
            qual_factor = 1 / 3
            if optwidth_workspace["proc"].ndim == 1:
                xdata = np.ravel(optwidth_workspace["proc"].coords)
                one_third = np.where(ydata > max(ydata) * qual_factor)
                one_third = np.ravel(one_third)

                best_width = xdata[one_third[-1]] - xdata[one_third[0]]

            else:
                xdata = np.ravel(optwidth_workspace["proc"].coords[0])
                min_x = []
                max_x = []
                for k in range(0, len(ydata[0, :])):
                    one_third = np.where(
                        ydata[
                            round(len(ydata[:, 0]) / 2)
                            - 75 : round(len(ydata[:, 0]) / 2)
                            + 75,
                            k,
                        ]
                        > max(
                            ydata[
                                round(len(ydata[:, 0]) / 2)
                                - 75 : round(len(ydata[:, 0]) / 2)
                                + 75,
                                k,
                            ]
                        )
                        * qual_factor
                    )
                    one_third = np.ravel(one_third)
                    min_x.append(xdata[one_third[0]])
                    max_x.append(xdata[one_third[len(one_third) - 1]])

                best_width = max(max_x) - min(min_x)

            self.gui_dict["processing_spec"]["integration_width"] = round(best_width)

            self.optCenter(self.gui_dict["processing_spec"]["integration_width"])
            self.optcentCheckbox.setChecked(True)

        if (
            self.gui_dict["gui_function"]["autoProcess"]
            or self.gui_dict["gui_function"]["T100_process"]
        ):
            self.gui_dict["processing_spec"]["phase"] = self.gui_dict[
                "processing_spec"
            ]["original_phase"]
        else:
            self.gui_dict["gui_function"]["sliders"] = False

            fac = np.pi / self.gui_dict["processing_spec"]["original_phase"]
            self.phaseSlider.setMinimum(round(-1000 * abs(fac)))
            self.phaseSlider.setMaximum(round(1000 * abs(fac)))

            if self.optphsCheckbox.isChecked():
                self.phaseSlider.setValue(
                    self.gui_dict["processing_spec"]["original_phase"]
                )

            self.intcenterSlider.setMinimum(
                self.gui_dict["processing_spec"]["integration_center"] - 50
            )
            self.intcenterSlider.setMaximum(
                self.gui_dict["processing_spec"]["integration_center"] + 50
            )

            if self.optcentCheckbox.isChecked():
                self.intcenterSlider.setValue(
                    self.gui_dict["processing_spec"]["integration_center"]
                )

            self.intwindowSlider.setValue(
                self.gui_dict["processing_spec"]["integration_width"]
            )
            self.intwindowEdit.setText(
                str(self.gui_dict["processing_spec"]["integration_width"])
            )

            self.gui_dict["gui_function"]["sliders"] = True

        if self.gui_dict["gui_function"]["T100_process"]:
            pass
        else:
            self.adjustSliders()

    def adjustSliders(self):

        adjslider_workspace = copy.deepcopy(self.processing_workspace)

        if self.gui_dict["gui_function"]["autoProcess"]:
            pass
        else:

            self.gui_dict["processing_spec"]["phase"] = self.gui_dict[
                "processing_spec"
            ]["original_phase"] + (
                self.gui_dict["processing_spec"]["phase_factor"]
                * self.gui_dict["processing_spec"]["original_phase"]
            )

            self.gui_dict["data_plot"]["xdata"] = np.reshape(
                adjslider_workspace["proc"].coords["t2"], -1
            )

            ydata = adjslider_workspace["proc"].values * np.exp(
                -1j * self.gui_dict["processing_spec"]["phase"]
            )
            self.gui_dict["data_plot"]["ydata"] = np.real(ydata)

        adjslider_workspace = self.phs_workspace(
            adjslider_workspace, self.gui_dict["processing_spec"]["phase"]
        )

        int_params = {
            "integrate_center": self.gui_dict["processing_spec"]["integration_center"],
            "integrate_width": self.gui_dict["processing_spec"]["integration_width"],
        }
        adjslider_workspace = self.int_workspace(adjslider_workspace, int_params)

        if len(adjslider_workspace["proc"].values) == 1:
            pass
        else:

            self.gui_dict["t1_fit"]["tau"] = np.reshape(
                adjslider_workspace["proc"].coords, -1
            )
            self.gui_dict["t1_fit"]["t1Amps"] = adjslider_workspace["proc"].values

            try:
                dnplab.dnpFit.t1Fit(adjslider_workspace)
            except:
                self.gui_dict["data_plot"]["xmin"] = int(
                    round(
                        self.gui_dict["processing_spec"]["integration_center"]
                        - np.abs(self.gui_dict["processing_spec"]["integration_width"])
                        / 2
                    )
                )
                self.gui_dict["data_plot"]["xmax"] = int(
                    round(
                        self.gui_dict["processing_spec"]["integration_center"]
                        + np.abs(self.gui_dict["processing_spec"]["integration_width"])
                        / 2
                    )
                )
                self.plot_data()
                self.gui_dict["t1_fit"]["xaxis"] = []
                self.gui_dict["t1_fit"]["t1Fit"] = []
                self.gui_dict["t1_fit"]["t1Val"] = 000
                self.gui_dict["enhancement_plot"]["title"] = "T1 Fit Error"
                self.plot_enh()
                return

            if self.gui_dict["gui_function"]["autoProcess"]:
                pass
            else:

                self.gui_dict["t1_fit"]["tau"] = np.reshape(
                    adjslider_workspace["proc"].coords["t1"], -1
                )
                self.gui_dict["t1_fit"]["t1Amps"] = adjslider_workspace["proc"].values

                self.gui_dict["t1_fit"]["xaxis"] = adjslider_workspace["fit"].coords[0]
                self.gui_dict["t1_fit"]["t1Fit"] = adjslider_workspace["fit"].values
                self.gui_dict["t1_fit"]["t1Val"] = adjslider_workspace["fit"].attrs[
                    "t1"
                ]

                self.plot_enh()

        if self.gui_dict["gui_function"]["autoProcess"]:
            pass
        else:

            if self.gui_dict["rawdata_function"]["folder"] == -1:
                print("---Standard Deviation in T1---")
                print(
                    "T1: "
                    + str(round(adjslider_workspace["fit"].attrs["t1"], 4))
                    + " +/- "
                    + str(round(adjslider_workspace["fit"].attrs["t1_stdd"], 4))
                )

            self.gui_dict["data_plot"]["xmin"] = int(
                round(
                    self.gui_dict["processing_spec"]["integration_center"]
                    - np.abs(self.gui_dict["processing_spec"]["integration_width"]) / 2
                )
            )
            self.gui_dict["data_plot"]["xmax"] = int(
                round(
                    self.gui_dict["processing_spec"]["integration_center"]
                    + np.abs(self.gui_dict["processing_spec"]["integration_width"]) / 2
                )
            )
            self.plot_data()

    def finishProcessing(self):

        self.gui_dict["gui_function"]["calculating"] = True
        self.show_hide_components()

        if (
            self.gui_dict["gui_function"]["isWorkup"]
            or self.gui_dict["gui_function"]["addWorkup"]
        ):

            wenh = np.array(
                [
                    self.gui_dict["workup_data"]["Epowers"],
                    self.gui_dict["workup_data"]["Ep"],
                ]
            )
            wenh = np.transpose(wenh)
            wenh = wenh[wenh[:, 0].argsort()]
            self.gui_dict["workup_data"]["Epowers"] = wenh[:, 0]
            self.gui_dict["workup_data"]["Ep"] = wenh[:, 1]

            wt1 = np.array(
                [
                    self.gui_dict["workup_data"]["T1powers"],
                    self.gui_dict["workup_data"]["T1p"],
                    self.gui_dict["workup_data"]["T1p_stdd"],
                ]
            )
            wt1 = np.transpose(wt1)
            wt1 = wt1[wt1[:, 0].argsort()]
            self.gui_dict["workup_data"]["T1powers"] = wt1[:, 0]
            self.gui_dict["workup_data"]["T1p"] = wt1[:, 1]
            self.gui_dict["workup_data"]["T1p_stdd"] = wt1[:, 2]

            if self.gui_dict["gui_function"]["addWorkup"]:
                self.show_wrkupCheckbox.setVisible(True)
                self.fit_wrkupCheckbox.setVisible(True)

        if self.gui_dict["gui_function"]["isWorkup"]:
            self.gui_dict["enhancement_plot"]["xdata"] = self.gui_dict["workup_data"][
                "Epowers"
            ]
            self.gui_dict["enhancement_plot"]["ydata"] = self.gui_dict["workup_data"][
                "Ep"
            ]

            self.gui_dict["enhancement_plot"]["ytick"] = [
                0,
                min(self.gui_dict["workup_data"]["Ep"]),
            ]

            if min(self.gui_dict["workup_data"]["Ep"]) <= -10:
                self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                    "0",
                    str(int(min(self.gui_dict["workup_data"]["Ep"]))),
                ]
            else:
                self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                    "0",
                    str(round(min(self.gui_dict["workup_data"]["Ep"]), 1)),
                ]

            self.gui_dict["t1_plot"]["xdata"] = self.gui_dict["workup_data"]["T1powers"]
            self.gui_dict["t1_plot"]["ydata"] = self.gui_dict["workup_data"]["T1p"]

        else:
            if self.gui_dict["rawdata_function"]["nopowers"]:
                self.gui_dict["dnpLab_data"]["Ep"] = self.Ep
                self.gui_dict["dnpLab_data"]["T1p"] = self.T1p
                self.gui_dict["dnpLab_data"]["T1p_stdd"] = self.T1p_stdd
            else:
                enh = np.array([self.gui_dict["dnpLab_data"]["Epowers"], self.Ep])
                enh = np.transpose(enh)
                enh = enh[enh[:, 0].argsort()]
                self.gui_dict["dnpLab_data"]["Epowers"] = enh[:, 0]
                self.gui_dict["dnpLab_data"]["Ep"] = enh[:, 1]

                t1 = np.array(
                    [self.gui_dict["dnpLab_data"]["T1powers"], self.T1p, self.T1p_stdd]
                )
                t1 = np.transpose(t1)
                t1 = t1[t1[:, 0].argsort()]
                self.gui_dict["dnpLab_data"]["T1powers"] = t1[:, 0]
                self.gui_dict["dnpLab_data"]["T1p"] = t1[:, 1]
                self.gui_dict["dnpLab_data"]["T1p_stdd"] = t1[:, 2]

                self.gui_dict["enhancement_plot"]["xdata"] = self.gui_dict[
                    "dnpLab_data"
                ]["Epowers"]
                self.gui_dict["enhancement_plot"]["ydata"] = self.gui_dict[
                    "dnpLab_data"
                ]["Ep"]

                self.gui_dict["enhancement_plot"]["ytick"] = [
                    0,
                    min(self.gui_dict["dnpLab_data"]["Ep"]),
                ]

                if min(self.gui_dict["dnpLab_data"]["Ep"]) <= -10:
                    self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                        "0",
                        str(int(min(self.gui_dict["dnpLab_data"]["Ep"]))),
                    ]
                else:
                    self.gui_dict["enhancement_plot"]["ytickLabel"] = [
                        "0",
                        str(round(min(self.gui_dict["dnpLab_data"]["Ep"]), 1)),
                    ]

                self.gui_dict["t1_plot"]["xdata"] = self.gui_dict["dnpLab_data"][
                    "T1powers"
                ]
                self.gui_dict["t1_plot"]["ydata"] = self.gui_dict["dnpLab_data"]["T1p"]

        self.gui_dict["enhancement_plot"]["title"] = "E[p]"
        self.gui_dict["enhancement_plot"]["xLabel"] = "microwave power"
        self.gui_dict["enhancement_plot"]["yLabel"] = "enhancement"
        self.gui_dict["t1_plot"]["title"] = r"$T_1[p]$"
        self.gui_dict["t1_plot"]["xLabel"] = "microwave power"
        self.gui_dict["t1_plot"]["yLabel"] = r"$T_1 (s)$"

        self.gui_dict["gui_function"]["sliders"] = False

        self.gui_dict["data_plot"]["plotksig"] = True
        self.gui_dict["t1_plot"]["plotT1interp"] = True
        self.gui_dict["enhancement_plot"]["plotT1fit"] = False
        self.gui_dict["enhancement_plot"]["plotEpfit"] = True

        if self.gui_dict["gui_function"]["isWorkup"]:
            print("---Workup Standard Deviations in T1s---")
            print(
                "T1(0): "
                + str(round(self.gui_dict["workup_data"]["T10"], 2))
                + " +/- "
                + str(round(self.gui_dict["workup_data"]["T10_stdd"], 4))
            )
            for k in range(0, len(self.gui_dict["workup_data"]["T1p"])):
                print(
                    str(round(self.gui_dict["workup_data"]["T1p"][k], 2))
                    + " +/- "
                    + str(round(self.gui_dict["workup_data"]["T1p_stdd"][k], 4))
                )

        else:
            print("---Standard Deviations in T1s---")
            print(
                "T1(0): "
                + str(round(self.gui_dict["dnpLab_data"]["T10"], 2))
                + " +/- "
                + str(round(self.gui_dict["dnpLab_data"]["T10_stdd"], 4))
            )
            for k in range(0, len(self.T1p)):
                print(
                    str(round(self.T1p[k], 2))
                    + " +/- "
                    + str(round(self.T1p_stdd[k], 4))
                )

            if self.gui_dict["gui_function"]["addWorkup"]:
                print("---Workup Standard Deviations in T1s---")
                print(
                    "T1(0): "
                    + str(round(self.gui_dict["workup_data"]["T10"], 2))
                    + " +/- "
                    + str(round(self.gui_dict["workup_data"]["T10_stdd"], 4))
                )
                for k in range(0, len(self.gui_dict["workup_data"]["T1p"])):
                    print(
                        str(round(self.gui_dict["workup_data"]["T1p"][k], 2))
                        + " +/- "
                        + str(round(self.gui_dict["workup_data"]["T1p_stdd"][k], 4))
                    )

        self.exclude1T1Checkbox.setChecked(False)

        self.Hydration_Calculator()

    def Hydration_Calculator(self):
        """Pass the processed data to the dnpHydration module.

        The GUI builds the input structure:

           dict =   {
                     'E' (numpy.array)        : signal enhancements,

                     'E_power' (numpy.array)  : microwave powers corresponding to                          the array 'E',

                     'T1' (numpy.array)       : T1 times,

                     'T1_power' (numpy.array) : microwave powers corresponding to                          the array 'T1',

                     'T10' (float)            : T1 time collected without microwave                        power,

                     'T100' (float)           : T1 time for a separate sample made                         without spin probe and collected                           without microwave power,

                     'spin_C' (float)         : concentration of                                           spin probe,

                     'field' (float)          : magnetic field setting for the                             experiment in units of mT,

                     'smax_model' (str)       : choice of model for setting s_max.                         Allowed values are 'tethered' where                        s_max=1 OR 'free' where s_max is                           calculated using spin_C,

                     't1_interp_method' (str) : choice of linear or second order                           interpolation of T1 onto E_power.                          Allowed values are 'linear' OR                             'second_order'.
                     }

        """

        self.dnpLab_errorLabel.setVisible(False)
        self.workup_errorLabel.setVisible(False)
        self.gui_dict["gui_function"]["hydrationEdits"] = True

        try:
            spin_C = float(self.slcEdit.text())
            field = float(self.fieldEdit.text())
            T100 = float(self.t100Edit.text())
            T10 = float(self.t10Edit.text())
        except:
            self.dnpLab_errorLabel.setVisible(True)
            print("Supply all parameters in numerical format")
            return

        if self.tetheredCheckbox.isChecked():
            smax_model = "tethered"
            self.wrkup_smax = 1
        else:
            smax_model = "free"
            self.wrkup_smax = 1 - (2 / (3 + (3 * (spin_C * 1e-6 * 198.7))))

        if self.linearfitCheckbox.isChecked():
            t1_interp_method = "linear"
        else:
            t1_interp_method = "second_order"

        if self.gui_dict["gui_function"]["isWorkup"]:
            self.workupt10Label.setVisible(True)
            self.workupt10Edit.setVisible(True)
            self.t10Label.setVisible(False)
            self.t10Edit.setVisible(False)
        else:

            if self.gui_dict["gui_function"]["addWorkup"]:
                self.workupt10Label.setVisible(True)
                self.workupt10Edit.setVisible(True)
            else:
                self.workupt10Label.setVisible(False)
                self.workupt10Edit.setVisible(False)

            if self.exclude1T1Checkbox.isChecked():

                T1p = self.gui_dict["dnpLab_data"]["T1p"][
                    1 : len(self.gui_dict["dnpLab_data"]["T1p"])
                ]
                T1powers = self.gui_dict["dnpLab_data"]["T1powers"][
                    1 : len(self.gui_dict["dnpLab_data"]["T1powers"])
                ]
            else:

                T1p = self.gui_dict["dnpLab_data"]["T1p"]
                T1powers = self.gui_dict["dnpLab_data"]["T1powers"]

            self.t10Label.setVisible(True)
            self.t10Edit.setVisible(True)
            self.gui_dict["dnpLab_data"]["T100"] = T100

            hydration = {
                "E": np.array(self.gui_dict["dnpLab_data"]["Ep"]),
                "E_power": np.array(self.gui_dict["dnpLab_data"]["Epowers"]),
                "T1": np.array(T1p),
                "T1_power": np.array(T1powers),
            }
            hydration.update(
                {
                    "T10": T10,
                    "T100": self.gui_dict["dnpLab_data"]["T100"],
                    "spin_C": spin_C,
                    "field": field,
                    "smax_model": smax_model,
                    "t1_interp_method": t1_interp_method,
                }
            )
            hyd = dnplab.create_workspace()
            hyd.add("hydration_inputs", hydration)

            try:
                self.gui_dict["hydration_results"] = dnplab.dnpHydration.hydration(hyd)
                self.addHyd_workspace = copy.deepcopy(hyd)
                self.addHyd_workspace.add(
                    "hydration_results", self.gui_dict["hydration_results"]
                )
            except:
                if T100 <= T10:
                    self.dnpLab_errorLabel.setText(
                        "DNPLab fit Error: T10(0) cannot be less than or equal to T1(0)"
                    )
                if spin_C <= 0:
                    self.dnpLab_errorLabel.setText(
                        "DNPLab fit Error: Spin concentration cannot be zero or negative"
                    )
                self.dataplt.axes.cla()
                self.dataplt.draw()
                self.dnpLab_errorLabel.setVisible(True)
                return

        if (
            self.gui_dict["gui_function"]["isWorkup"]
            or self.gui_dict["gui_function"]["addWorkup"]
        ):

            if (
                self.gui_dict["workup_function"]["fit"]
                or self.gui_dict["workup_function"]["show"]
            ):

                if self.exclude1T1Checkbox.isChecked():

                    wT1p = self.gui_dict["workup_data"]["T1p"][
                        1 : len(self.gui_dict["workup_data"]["T1p"])
                    ]
                    wT1powers = self.gui_dict["workup_data"]["T1powers"][
                        1 : len(self.gui_dict["workup_data"]["T1powers"])
                    ]

                else:

                    wT1p = self.gui_dict["workup_data"]["T1p"]
                    wT1powers = self.gui_dict["workup_data"]["T1powers"]

                try:
                    wT10 = float(self.workupt10Edit.text())
                except:
                    self.workup_errorLabel.setVisible(True)
                    print("Supply all parameters in numerical format")
                    return

                self.gui_dict["workup_data"]["T100"] = T100

                whydration = {
                    "E": np.array(self.gui_dict["workup_data"]["Ep"]),
                    "E_power": np.array(self.gui_dict["workup_data"]["Epowers"]),
                    "T1": np.array(wT1p),
                    "T1_power": np.array(wT1powers),
                }
                whydration.update(
                    {
                        "T10": wT10,
                        "T100": self.gui_dict["workup_data"]["T100"],
                        "spin_C": spin_C,
                        "field": field,
                        "smax_model": smax_model,
                        "t1_interp_method": t1_interp_method,
                    }
                )
                whyd = dnplab.create_workspace()
                whyd.add("hydration_inputs", whydration)

                try:
                    self.gui_dict[
                        "workup_hydration_results"
                    ] = dnplab.dnpHydration.hydration(whyd)
                    if (
                        self.gui_dict["workup_function"]["fit"]
                        or self.gui_dict["gui_function"]["isWorkup"]
                    ):
                        self.addHyd_workspace = copy.deepcopy(whyd)
                        self.addHyd_workspace.add(
                            "hydration_results",
                            self.gui_dict["workup_hydration_results"],
                        )
                except:
                    if T100 <= wT10:
                        self.workup_errorLabel.setText(
                            "Workup fit Error: T10(0) cannot be less than or equal to T1(0)"
                        )
                    if spin_C <= 0:
                        self.workup_errorLabel.setText(
                            "Workup fit Error: Spin concentration cannot be zero or negative"
                        )
                    self.dataplt.axes.cla()
                    self.dataplt.draw()
                    self.workup_errorLabel.setVisible(True)
                    return

        if self.gui_dict["gui_function"]["isWorkup"]:

            if min(wT1p) < 0.1:
                self.gui_dict["t1_plot"]["ymin"] = 0
            else:
                self.gui_dict["t1_plot"]["ymin"] = (
                    min(self.gui_dict["workup_data"]["T1p"]) * 0.85
                )

            if max(wT1p) > 5:
                self.gui_dict["t1_plot"]["ymax"] = 1
            else:
                self.gui_dict["t1_plot"]["ymax"] = (
                    max(self.gui_dict["workup_data"]["T1p"]) * 1.15
                )

        else:
            if min(T1p) < 0.1:
                self.gui_dict["t1_plot"]["ymin"] = 0
            else:
                self.gui_dict["t1_plot"]["ymin"] = (
                    min(self.gui_dict["dnpLab_data"]["T1p"]) * 0.85
                )

            if max(T1p) > 5:
                self.gui_dict["t1_plot"]["ymax"] = 1
            else:
                self.gui_dict["t1_plot"]["ymax"] = (
                    max(self.gui_dict["dnpLab_data"]["T1p"]) * 1.15
                )

        self.gui_dict["t1_plot"]["ytick"] = [
            self.gui_dict["t1_plot"]["ymin"],
            self.gui_dict["t1_plot"]["ymax"],
        ]
        self.gui_dict["t1_plot"]["ytickLabel"] = [
            str(round(self.gui_dict["t1_plot"]["ymin"], 1)),
            str(round(self.gui_dict["t1_plot"]["ymax"], 1)),
        ]

        self.gui_dict["data_plot"]["title"] = r"$k_\sigma[p]$"

        self.plot_data()
        self.plot_enh()
        self.plot_t1()

        print("-----Standard Deviation in ksigma-----")
        if self.gui_dict["gui_function"]["isWorkup"]:
            print(
                "Workup (dnpHydration): "
                + str(round(self.gui_dict["workup_hydration_results"]["ksigma"], 2))
                + " +/- "
                + str(
                    round(self.gui_dict["workup_hydration_results"]["ksigma_stdd"], 4)
                )
            )
            print(
                "Workup = "
                + str(
                    round(
                        self.gui_dict["workup_data"]["kSigma"]
                        / spin_C
                        / self.wrkup_smax,
                        2,
                    )
                )
                + " +/- "
                + str(
                    round(
                        self.gui_dict["workup_data"]["kSigma_stdd"]
                        / spin_C
                        / self.wrkup_smax,
                        4,
                    )
                )
            )
        else:
            print(
                "DNPLab = "
                + str(round(self.gui_dict["hydration_results"]["ksigma"], 2))
                + " +/- "
                + str(round(self.gui_dict["hydration_results"]["ksigma_stdd"], 4))
            )
            if self.gui_dict["workup_function"]["fit"]:
                print(
                    "Workup (dnpHydration) = "
                    + str(round(self.gui_dict["workup_hydration_results"]["ksigma"], 2))
                    + " +/- "
                    + str(
                        round(
                            self.gui_dict["workup_hydration_results"]["ksigma_stdd"], 4
                        )
                    )
                )
            if self.gui_dict["gui_function"]["addWorkup"]:
                print(
                    "Workup = "
                    + str(
                        round(
                            self.gui_dict["workup_data"]["kSigma"]
                            / spin_C
                            / self.wrkup_smax,
                            2,
                        )
                    )
                    + " +/- "
                    + str(
                        round(
                            self.gui_dict["workup_data"]["kSigma_stdd"]
                            / spin_C
                            / self.wrkup_smax,
                            4,
                        )
                    )
                )

    def Save_Results_Button(self):
        """Save the results of processing to a format that can be read by the hydrationGUI using the 'GUI Result' button or by the MATLAB App called xODNP."""

        pthnm1 = QFileDialog.getSaveFileName(self)
        if pthnm1[0]:
            pthnm = pthnm1[0]
        else:
            return

        if (
            self.gui_dict["workup_function"]["fit"]
            or self.gui_dict["gui_function"]["isWorkup"]
        ):
            self.addHyd_workspace["hydration_results"].update(
                {
                    "T1_stdd": self.gui_dict["workup_data"]["T1p_stdd"],
                    "T10_stdd": self.gui_dict["workup_data"]["T10_stdd"],
                }
            )
        else:
            self.addHyd_workspace["hydration_results"].update(
                {
                    "T1_stdd": self.gui_dict["dnpLab_data"]["T1p_stdd"],
                    "T10_stdd": self.gui_dict["dnpLab_data"]["T10_stdd"],
                }
            )

        odnpData = {
            "Epowers": self.addHyd_workspace["hydration_inputs"]["E_power"],
            "Ep": self.addHyd_workspace["hydration_inputs"]["E"],
            "T1powers": self.addHyd_workspace["hydration_inputs"]["T1_power"],
            "T1p": self.addHyd_workspace["hydration_inputs"]["T1"],
            "T1p_stdd": self.addHyd_workspace["hydration_results"]["T1_stdd"],
            "T10": self.addHyd_workspace["hydration_inputs"]["T10"],
            "T10_stdd": self.addHyd_workspace["hydration_results"]["T10_stdd"],
            "T100": self.addHyd_workspace["hydration_inputs"]["T100"],
        }

        odnpResults = {
            "kSigmas": self.addHyd_workspace["hydration_results"]["ksigma_array"],
            "kSigmas_fit": self.addHyd_workspace["hydration_results"]["ksigma_fit"],
        }

        spltpthnm = pthnm.split(os.sep)
        flnm = spltpthnm[-1]
        svpthnm = pthnm + " hydrationGUI Results"
        if os.path.isdir(svpthnm):
            svpthnm = pthnm + "_COPY" + " hydrationGUI Results"
        os.mkdir(svpthnm)

        print("Save name: " + flnm)
        print("Save path: " + svpthnm)

        dnplab.dnpImport.h5.saveh5(
            self.addHyd_workspace,
            os.path.join(svpthnm, flnm + " hydration_parameters.h5"),
        )

        savemat(
            os.path.join(svpthnm, flnm + " xODNP.mat"),
            {"odnp": odnpData, "ksig": odnpResults},
            oned_as="column",
        )

        dfE = np.vstack(
            (
                self.addHyd_workspace["hydration_inputs"]["E_power"],
                self.addHyd_workspace["hydration_inputs"]["E"],
                self.addHyd_workspace["hydration_results"]["ksigma_array"],
                self.addHyd_workspace["hydration_results"]["ksigma_fit"],
            )
        ).T
        np.savetxt(
            os.path.join(svpthnm, flnm + " E_ksig.csv"),
            dfE,
            fmt="%10.10f",
            delimiter=",",
            header="E powers,E(p),ksigma(p),ksigma(p) fit",
            comments="",
        )

        dfT1 = np.vstack(
            (
                self.addHyd_workspace["hydration_inputs"]["T1_power"],
                self.addHyd_workspace["hydration_inputs"]["T1"],
                self.addHyd_workspace["hydration_results"]["T1_stdd"][
                    0 : len(self.addHyd_workspace["hydration_inputs"]["T1_power"])
                ],
            )
        ).T
        np.savetxt(
            os.path.join(svpthnm, flnm + " T1s.csv"),
            dfT1,
            fmt="%10.10f",
            delimiter=",",
            header="T1 powers,T1(p),T1(p) Std dev",
            comments="",
        )

        dfPars = np.vstack(
            (
                self.addHyd_workspace["hydration_inputs"]["spin_C"],
                self.addHyd_workspace["hydration_inputs"]["field"],
                self.addHyd_workspace["hydration_inputs"]["T100"],
                self.addHyd_workspace["hydration_inputs"]["T10"],
                self.addHyd_workspace["hydration_results"]["T10_stdd"],
                self.addHyd_workspace["hydration_results"]["krho"],
                self.addHyd_workspace["hydration_results"]["ksigma"],
                self.addHyd_workspace["hydration_results"]["ksigma_stdd"],
                self.addHyd_workspace["hydration_results"]["klow"],
                self.addHyd_workspace["hydration_results"]["coupling_factor"],
                self.addHyd_workspace["hydration_results"]["tcorr"],
                self.addHyd_workspace["hydration_results"]["Dlocal"],
            )
        ).T
        np.savetxt(
            os.path.join(svpthnm, flnm + " params.csv"),
            dfPars,
            fmt="%10.2f,%10.2f,%10.4f,%10.4f,%10.4f,%10.2f,%10.2f,%10.4f,%10.2f,%1.5f,%10.2f,%1.3e",
            delimiter=",",
            header="spin concentration (uM),Field (mT),T10(0) (s),T1(0) (s),T1(0) stdd, krho (s^-1M^-1),ksigma (s^-1M^-1),ksigma stdd,klow (s^-1M^-1),coupling factor,tcorr (ps),Dlocal (m^2/s)",
            comments="",
        )

    def Spectrum_Phase_Slider(self, pvalue):
        """Slider to change the phase correction applied to the spectrum."""
        if self.gui_dict["gui_function"]["sliders"]:
            self.gui_dict["processing_spec"]["phase_factor"] = pvalue / 1000
            self.optphsCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass

    def Integration_Center_Slider(self, cvalue):
        """Slider to change the center of the spectrum integration window."""
        if self.gui_dict["gui_function"]["sliders"]:
            self.gui_dict["processing_spec"]["integration_center"] = cvalue
            self.optcentCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass

    def Integration_Window_Slider(self, wvalue):
        """Slider to change the width of the spectrum integration window."""
        if self.gui_dict["gui_function"]["sliders"]:
            self.gui_dict["processing_spec"]["integration_width"] = wvalue
            self.intwindowEdit.setText(str(wvalue))
            self.optwidthCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass

    def Integration_Window_Edit(self):
        """This function passes the text from the various edit boxes to dnpHydration as floats and re-calculates
        hydration parameters."""
        if self.gui_dict["gui_function"]["sliders"]:
            int_wind = float(self.intwindowEdit.text()) + 0.1
            self.gui_dict["processing_spec"]["integration_width"] = round(int_wind)
            self.intwindowSlider.setValue(
                self.gui_dict["processing_spec"]["integration_width"]
            )
            self.intwindowEdit.setText(
                str(self.gui_dict["processing_spec"]["integration_width"])
            )
        else:
            pass

    def Optimize_Phase_Checkbox(self):
        """Check this to have the GUI automatically choose the best phase."""
        if self.gui_dict["gui_function"]["sliders"]:
            if self.optphsCheckbox.isChecked():
                self.gui_dict["processing_spec"]["phase_factor"] = 0
                self.processData()
            else:
                pass
        else:
            pass

    def Optimize_Center_Checkbox(self):
        """Check this to have the GUI automatically choose the best integration center."""
        if self.gui_dict["gui_function"]["sliders"]:
            if self.optcentCheckbox.isChecked():
                self.processData()
            else:
                pass
        else:
            pass

    def Optimize_Width_Checkbox(self):
        """Check this to have the GUI automatically choose the best integration width."""
        if (
            self.gui_dict["gui_function"]["sliders"]
            and self.gui_dict["rawdata_function"]["folder"] != -2
        ):
            if self.optwidthCheckbox.isChecked():
                pass
            else:
                self.gui_dict["processing_spec"]["integration_width"] = 10

            self.processData()

        else:
            pass

    def Linear_Interpolation_Checkbox(self):
        """Choose a linear T1 interpolation."""
        if self.linearfitCheckbox.isChecked():
            self.order2fitCheckbox.setChecked(False)
        else:
            self.order2fitCheckbox.setChecked(True)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def SecondOrder_Interpolation_Checkbox(self):
        """Choose a second order T1 interpolation."""
        if self.order2fitCheckbox.isChecked():
            self.linearfitCheckbox.setChecked(False)
        else:
            self.linearfitCheckbox.setChecked(True)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Exclude_FirstT1_Checkbox(self):
        """Exclude the first T1 point from the interpolation if it deviates significantly from the trend of the other
        points."""
        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Smax_Tethered_Checkbox(self):
        """Choose s_max = 1"""
        if self.tetheredCheckbox.isChecked():
            self.freeCheckbox.setChecked(False)
        else:
            self.freeCheckbox.setChecked(True)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Smax_Free_Checkbox(self):
        """Choose to have s_max calculated based on spin probe concentration."""
        if self.freeCheckbox.isChecked():
            self.tetheredCheckbox.setChecked(False)
        else:
            self.tetheredCheckbox.setChecked(True)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Edit_Hydration_Inputs(self):
        """This function passes the text from the various edit boxes to dnpHydration as floats and re-calculates
        hydration parameters."""
        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Show_Workup_Checkbox(self):
        """Show or hide the Workup results if they were found in the data folder."""
        if self.show_wrkupCheckbox.isChecked():
            self.gui_dict["workup_function"]["show"] = True
        else:
            self.gui_dict["workup_function"]["show"] = False
            self.gui_dict["workup_function"]["fit"] = False
            self.fit_wrkupCheckbox.setChecked(False)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Fit_Workup_Checkbox(self):
        """Use dnpHydration to analyze the results of the workup code processing."""
        if self.fit_wrkupCheckbox.isChecked() and self.show_wrkupCheckbox.isChecked():
            self.gui_dict["workup_function"]["fit"] = True
        else:
            self.gui_dict["workup_function"]["fit"] = False
            self.fit_wrkupCheckbox.setChecked(False)

        if self.gui_dict["gui_function"]["hydrationEdits"]:
            self.Hydration_Calculator()
        else:
            pass

    def Only_T1p_Checkbox(self):
        """Rather than return to the beginning of the E(p) series, the Restart button will return to the first T1(p)
        point."""
        if self.onlyT1pCheckbox.isChecked():
            self.onlyT10Checkbox.setChecked(False)

    def Only_T10_Checkbox(self):
        """Rather than return to the beginning of the E(p) series, the Restart button will return to the T1(0) point."""
        if self.onlyT10Checkbox.isChecked():
            self.onlyT1pCheckbox.setChecked(False)

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

        if self.gui_dict["data_plot"]["plotksig"]:

            indx_h = (
                max(self.addHyd_workspace["hydration_results"]["ksigma_array"]) * 0.8
            )
            self.dataplt.axes.set_ylim(
                0, max(self.addHyd_workspace["hydration_results"]["ksigma_array"]) * 1.1
            )

            indexes = [0.11, 0.21, 0.31, 0.41, 0.51, 0.61]
            if (
                self.gui_dict["gui_function"]["addWorkup"]
                and self.gui_dict["workup_function"]["show"]
            ):
                self.dataplt.axes.plot(
                    self.gui_dict["workup_data"]["Epowers"],
                    self.gui_dict["workup_hydration_results"]["ksigma_array"],
                    color="#003660",
                    marker="o",
                    linestyle="none",
                    label=r"Workup $k_\sigma$[p]",
                )

                self.dataplt.axes.text(
                    max(self.gui_dict["dnpLab_data"]["Epowers"]) * 0.645,
                    indx_h - (0.21 * indx_h),
                    r"Workup $k_\sigma = $"
                    + str(
                        round(
                            self.gui_dict["workup_data"]["kSigma"]
                            / float(self.slcEdit.text())
                            / self.wrkup_smax,
                            2,
                        )
                    ),
                    fontsize=12,
                )
                indexes = [0, 0.11, 0.31, 0.41, 0.51, 0.61]

            if self.gui_dict["gui_function"]["isWorkup"]:
                self.dataplt.axes.plot(
                    self.gui_dict["workup_data"]["Epowers"],
                    self.gui_dict["workup_hydration_results"]["ksigma_array"],
                    color="#003660",
                    marker="o",
                    linestyle="none",
                    label=r"Workup $k_\sigma$[p]",
                )
                self.dataplt.axes.plot(
                    self.gui_dict["workup_data"]["Epowers"],
                    self.gui_dict["workup_hydration_results"]["ksigma_fit"],
                    color="#F37021",
                    label=r"dnpHydration Fit",
                )
            else:
                self.dataplt.axes.plot(
                    self.gui_dict["dnpLab_data"]["Epowers"],
                    self.gui_dict["hydration_results"]["ksigma_array"],
                    color="#46812B",
                    marker="o",
                    linestyle="none",
                    label=r"DNPLab $k_\sigma$[p]",
                )
                if self.gui_dict["workup_function"]["fit"]:
                    self.dataplt.axes.plot(
                        self.gui_dict["workup_data"]["Epowers"],
                        self.gui_dict["workup_hydration_results"]["ksigma_fit"],
                        color="#F37021",
                        label="Workup Fit",
                    )
                else:
                    self.dataplt.axes.plot(
                        self.gui_dict["dnpLab_data"]["Epowers"],
                        self.gui_dict["hydration_results"]["ksigma_fit"],
                        color="#F37021",
                        label="DNPLab Fit",
                    )

            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[0] * indx_h),
                r"$k_\rho = $"
                + str(round(self.addHyd_workspace["hydration_results"]["krho"], 2)),
                fontsize=12,
            )
            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[1] * indx_h),
                r"$k_\sigma = $"
                + str(round(self.addHyd_workspace["hydration_results"]["ksigma"], 2)),
                fontsize=12,
            )
            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[2] * indx_h),
                r"$k_{low} = $"
                + str(round(self.addHyd_workspace["hydration_results"]["klow"], 2)),
                fontsize=12,
            )
            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[3] * indx_h),
                r"$\xi = $"
                + str(
                    round(
                        self.addHyd_workspace["hydration_results"]["coupling_factor"], 4
                    )
                ),
                fontsize=12,
            )
            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[4] * indx_h),
                r"$t_{corr} = $"
                + str(round(self.addHyd_workspace["hydration_results"]["tcorr"], 2)),
                fontsize=12,
            )
            d_local = round(
                self.addHyd_workspace["hydration_results"]["Dlocal"] * 1e10, 2
            )
            self.dataplt.axes.text(
                max(self.addHyd_workspace["hydration_inputs"]["E_power"]) * 0.75,
                indx_h - (indexes[5] * indx_h),
                r"$D_{local} = $" + str(d_local) + r"$e^{-10}$",
                fontsize=12,
            )
            self.dataplt.axes.set_yticks(
                [0, max(self.addHyd_workspace["hydration_results"]["ksigma_array"])]
            )
            self.dataplt.axes.set_yticklabels(
                [
                    "0",
                    str(
                        round(
                            max(
                                self.addHyd_workspace["hydration_results"][
                                    "ksigma_array"
                                ]
                            ),
                            1,
                        )
                    ),
                ]
            )
            self.dataplt.axes.set_xticks([])
            self.dataplt.axes.set_xlabel("microwave power")
            self.dataplt.axes.set_ylabel(r"$k_\sigma[p]$")
            self.dataplt.axes.legend()
        else:
            self.dataplt.axes.plot(
                self.gui_dict["data_plot"]["xdata"], self.gui_dict["data_plot"]["ydata"]
            )
            self.dataplt.axes.set_xlim(
                self.gui_dict["data_plot"]["xmin"], self.gui_dict["data_plot"]["xmax"]
            )
            self.dataplt.axes.set_xticks(
                [self.gui_dict["data_plot"]["xmin"], self.gui_dict["data_plot"]["xmax"]]
            )
            self.dataplt.axes.set_xticklabels(
                [
                    str(self.gui_dict["data_plot"]["xmin"]),
                    str(self.gui_dict["data_plot"]["xmax"]),
                ]
            )
            self.dataplt.axes.set_xlabel("ppm")
            self.dataplt.axes.set_yticks([0])
            self.dataplt.axes.set_yticklabels("0")

        self.dataplt.axes.set_title(self.gui_dict["data_plot"]["title"])
        self.dataplt.draw()

    def plot_enh(self):

        self.enhplt.axes.cla()
        if self.gui_dict["enhancement_plot"]["plotT1fit"]:
            self.enhplt.axes.plot(
                self.gui_dict["t1_fit"]["tau"],
                self.gui_dict["t1_fit"]["t1Amps"],
                color="#46812B",
                marker="o",
                linestyle="none",
            )
            self.enhplt.axes.plot(
                self.gui_dict["t1_fit"]["xaxis"],
                self.gui_dict["t1_fit"]["t1Fit"],
                "#F37021",
            )
            self.enhplt.axes.text(
                max(self.gui_dict["t1_fit"]["tau"]) * 0.55,
                max(self.gui_dict["t1_fit"]["t1Amps"]) * 0.3,
                r"$T_1$ =" + str(round(self.gui_dict["t1_fit"]["t1Val"], 4)) + " s",
                fontsize=10,
            )
        else:
            if self.gui_dict["enhancement_plot"]["plotEpfit"]:
                if (
                    self.gui_dict["gui_function"]["addWorkup"]
                    and self.gui_dict["workup_function"]["show"]
                ):
                    self.enhplt.axes.plot(
                        self.gui_dict["workup_data"]["Epowers"],
                        self.gui_dict["workup_data"]["Ep"],
                        color="#003660",
                        marker="o",
                        linestyle="none",
                        label="Workup",
                    )

                if self.gui_dict["gui_function"]["isWorkup"]:
                    self.enhplt.axes.plot(
                        self.gui_dict["workup_data"]["Epowers"],
                        self.gui_dict["workup_data"]["Ep"],
                        color="#003660",
                        marker="o",
                        linestyle="none",
                        label="Workup",
                    )
                    self.enhplt.axes.plot(
                        self.gui_dict["workup_data"]["Epowers"],
                        self.gui_dict["workup_hydration_results"]["uncorrected_Ep"],
                        color="#F37021",
                        label="dnpHydration Fit",
                    )
                else:
                    self.enhplt.axes.plot(
                        self.gui_dict["dnpLab_data"]["Epowers"],
                        self.gui_dict["dnpLab_data"]["Ep"],
                        color="#46812B",
                        marker="o",
                        linestyle="none",
                        label="DNPLab",
                    )
                    if self.gui_dict["workup_function"]["fit"]:
                        self.enhplt.axes.plot(
                            self.gui_dict["workup_data"]["Epowers"],
                            self.gui_dict["workup_hydration_results"]["uncorrected_Ep"],
                            color="#F37021",
                            label="Workup Fit",
                        )
                    else:
                        self.enhplt.axes.plot(
                            self.gui_dict["dnpLab_data"]["Epowers"],
                            self.gui_dict["hydration_results"]["uncorrected_Ep"],
                            color="#F37021",
                            label="DNPLab Fit",
                        )

                self.enhplt.axes.legend()
            else:
                self.enhplt.axes.plot(
                    self.gui_dict["enhancement_plot"]["xdata"],
                    self.gui_dict["enhancement_plot"]["ydata"],
                    color="#46812B",
                    marker="o",
                    linestyle="none",
                )

        self.enhplt.axes.set_title(self.gui_dict["enhancement_plot"]["title"])
        self.enhplt.axes.set_xlabel(self.gui_dict["enhancement_plot"]["xLabel"])
        self.enhplt.axes.set_xticks([])
        self.enhplt.axes.set_ylabel(self.gui_dict["enhancement_plot"]["yLabel"])
        self.enhplt.axes.set_yticks(self.gui_dict["enhancement_plot"]["ytick"])
        self.enhplt.axes.set_yticklabels(
            self.gui_dict["enhancement_plot"]["ytickLabel"]
        )
        self.enhplt.draw()

    def plot_t1(self):

        self.t1plt.axes.cla()
        if self.gui_dict["t1_plot"]["plotT1interp"]:
            if (
                self.gui_dict["gui_function"]["addWorkup"]
                and self.gui_dict["workup_function"]["show"]
            ):
                self.t1plt.axes.plot(
                    self.gui_dict["workup_data"]["T1powers"],
                    self.gui_dict["workup_data"]["T1p"],
                    color="#003660",
                    marker="o",
                    linestyle="none",
                    label="Workup",
                )

            if self.gui_dict["gui_function"]["isWorkup"]:
                self.t1plt.axes.plot(
                    self.gui_dict["workup_data"]["T1powers"],
                    self.gui_dict["workup_data"]["T1p"],
                    color="#003660",
                    marker="o",
                    linestyle="none",
                    label="Workup",
                )
                self.t1plt.axes.plot(
                    self.gui_dict["workup_data"]["Epowers"],
                    self.gui_dict["workup_hydration_results"]["interpolated_T1"],
                    "#F37021",
                    label="Interpolation",
                )
            else:
                self.t1plt.axes.plot(
                    self.gui_dict["dnpLab_data"]["T1powers"],
                    self.gui_dict["dnpLab_data"]["T1p"],
                    color="#46812B",
                    marker="o",
                    linestyle="none",
                    label="DNPLab",
                )
                if self.gui_dict["workup_function"]["fit"]:
                    self.t1plt.axes.plot(
                        self.gui_dict["workup_data"]["Epowers"],
                        self.gui_dict["workup_hydration_results"]["interpolated_T1"],
                        "#F37021",
                        label="Interpolation",
                    )
                else:
                    self.t1plt.axes.plot(
                        self.gui_dict["dnpLab_data"]["Epowers"],
                        self.gui_dict["hydration_results"]["interpolated_T1"],
                        "#F37021",
                        label="Interpolation",
                    )
            self.t1plt.axes.legend()
        else:
            self.t1plt.axes.plot(
                self.gui_dict["t1_plot"]["xdata"],
                self.gui_dict["t1_plot"]["ydata"],
                color="#46812B",
                marker="o",
                linestyle="none",
            )

        self.t1plt.axes.set_title(self.gui_dict["t1_plot"]["title"])
        self.t1plt.axes.set_xlabel(self.gui_dict["t1_plot"]["xLabel"])
        self.t1plt.axes.set_xticks([])
        self.t1plt.axes.set_ylabel(self.gui_dict["t1_plot"]["yLabel"])
        self.t1plt.axes.set_ylim(
            self.gui_dict["t1_plot"]["ymin"], self.gui_dict["t1_plot"]["ymax"]
        )
        self.t1plt.axes.set_yticks(self.gui_dict["t1_plot"]["ytick"])
        self.t1plt.axes.set_yticklabels(self.gui_dict["t1_plot"]["ytickLabel"])
        self.t1plt.draw()


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
    ex = hydrationGUI()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main_func()
