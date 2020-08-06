import sys
import os

from PyQt5 import QtWidgets,QtGui,QtCore

from PyQt5.QtWidgets import QApplication, QMainWindow, QSizePolicy, QWidget, QPushButton, QLineEdit, QSlider, QLabel, QCheckBox, QFileDialog, QLineEdit
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np
from scipy.io import loadmat, savemat
from scipy import optimize
import pandas as pd
import copy
import re
import datetime
import time

sys.path.append('...dnpLab')
import dnpLab as odnp

class hydrationGUI(QMainWindow):
    """Interactive tool for processing Han lab ODNP NMR data and conversion to enhancement or T1 as functions of microwave power. The processed data is used to calculate hydration parameters.

    Uses the odnpImport, odnpData, odnpNMR, and odnpFit modules to import, process, integrate, and fit (T1 data) then passes the results to the hydration module for analysis. Also supplies an output of the results in .mat format for analysis in MATALB using xODNP.

    Args:
        Selected data folder: Select the base folder of Bruker format data collected using the Han lab dnp routine. Specify the folder structure using the dropdown.

    Returns:
        A GUI for interactive processing and optional *.mat output compatible with the MATLAB based xODNP App
    """
    
    def __init__(self):
    
        super().__init__()
        self.left = 10
        self.top = 10
        self.title = 'Han Lab ODNP Processing'
        self.width = 1050
        self.height = 625
        
        #self.setStyleSheet('background-color : rgb(255,255,255)')
        
        self.dataplt = PlotCanvas(self, width=7.2, height=4.8)
        self.enhplt = PlotCanvas(self, width=3.15, height=2)
        self.t1plt = PlotCanvas(self, width=3.15, height=2)
        
        self.initUI()

    def initUI(self):
        
        self.testmode = False # set to True for testing, False for normal use
        self.labpath = '...dnplab' # same as sys path to dnpLab
         
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setContentsMargins(0,0,0,0)
        
        # main plot
        self.dataplt.move(5,40)
        # Enh plot
        self.enhplt.move(730,40)
        # t1 plot
        self.t1plt.move(730,260)
        
        # Create a load hydrationGUI button
        self.hanlabButton = QPushButton('GUI Result', self)
        self.hanlabButton.setStyleSheet('font : bold ; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.hanlabButton.move(5,5)
        self.hanlabButton.resize(80,30)
        self.hanlabButton.clicked.connect(self.loadLab)
        
        # Create a load workup button
        self.workupButton = QPushButton('Workup', self)
        self.workupButton.setStyleSheet('font : bold ; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.workupButton.move(90,5)
        self.workupButton.resize(80,30)
        self.workupButton.clicked.connect(self.loadWorkup)
        
        # Create a load single button
        self.singleButton = QPushButton('Bruker', self)
        self.singleButton.setStyleSheet('font : bold ; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.singleButton.move(175,5)
        self.singleButton.resize(80,30)
        self.singleButton.clicked.connect(self.loadSingle)
        
        # Create a load button
        self.rawdataButton = QPushButton('Han Lab', self)
        self.rawdataButton.setStyleSheet('font : bold ; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.rawdataButton.move(260,5)
        self.rawdataButton.resize(80,30)
        self.rawdataButton.clicked.connect(self.loadRaw)

        # Dataname label
        self.pathLabel = QLabel(self)
        self.pathLabel.setStyleSheet('font : bold 14px; color : rgb(0, 77, 159)')
        self.pathLabel.move(345, 13)
        self.pathLabel.resize(700,20)
        self.pathLabel.setText('Data folder path')
        
        # int center slider label
        self.intcenterLabel = QLabel(self)
        self.intcenterLabel.setStyleSheet('font : bold 14px') #; color : rgb(255,184,20)') ; background-color : rgb(0, 77, 159)')
        self.intcenterLabel.move(123, 520)
        self.intcenterLabel.resize(490,30)
        self.intcenterLabel.setText(' Window center:')
        # int center slider
        self.intcenterSlider = QSlider(Qt.Horizontal, self)
        #self.intcenterSlider.setStyleSheet('background-color : rgb(0, 77, 159)')
        self.intcenterSlider.setGeometry(250, 525, 360, 20)
        
        # int window slider label
        self.intwindowLabel = QLabel(self)
        self.intwindowLabel.setStyleSheet('font : bold 14px') #; color : rgb(255,184,20)') ; background-color : rgb(0, 77, 159)')
        self.intwindowLabel.move(123, 555)
        self.intwindowLabel.resize(490,30)
        self.intwindowLabel.setText('  Window width:')
        # int window slider
        self.intwindowSlider = QSlider(Qt.Horizontal, self)
        #self.intwindowSlider.setStyleSheet('background-color : rgb(0, 77, 159)')
        self.intwindowSlider.setGeometry(250, 560, 360, 20)
        
        # Phase slider label
        self.phaseLabel = QLabel(self)
        self.phaseLabel.setStyleSheet('font : bold 14px') #; color : rgb(255,184,20)') ; background-color : rgb(0, 77, 159)')
        self.phaseLabel.move(123, 590)
        self.phaseLabel.resize(490,30)
        self.phaseLabel.setText('   Phase Adjust:')
        # Phase slider
        self.phaseSlider = QSlider(Qt.Horizontal, self)
        #self.phaseSlider.setStyleSheet('background-color : rgb(0, 77, 159)')
        self.phaseSlider.setGeometry(250, 595, 360, 20)
        # autophase checkbox
        self.autophaseCheckbox = QCheckBox(self)
        self.autophaseCheckbox.setStyleSheet('font : bold 14px')
        self.autophaseCheckbox.move(612,595)
        self.autophaseCheckbox.resize(100,20)
        self.autophaseCheckbox.setText('Optimize')
        self.autophaseCheckbox.setChecked(False)
        
        # show workup checkbox
        self.show_wrkupCheckbox = QCheckBox(self)
        self.show_wrkupCheckbox.setStyleSheet('font : bold 14px')
        self.show_wrkupCheckbox.move(627,565)
        self.show_wrkupCheckbox.resize(130,20)
        self.show_wrkupCheckbox.setText('Show workup')
        # fit workup checkbox
        self.fit_wrkupCheckbox = QCheckBox(self)
        self.fit_wrkupCheckbox.setStyleSheet('font : bold 14px')
        self.fit_wrkupCheckbox.move(627,585)
        self.fit_wrkupCheckbox.resize(130,20)
        self.fit_wrkupCheckbox.setText('Fit workup')
        
        # autophase checkbox
        self.onlyT1pCheckbox = QCheckBox(self)
        self.onlyT1pCheckbox.setStyleSheet('font : bold 14px')
        self.onlyT1pCheckbox.move(10,565)
        self.onlyT1pCheckbox.resize(100,20)
        self.onlyT1pCheckbox.setText('Only T1(p)')
        self.onlyT1pCheckbox.setChecked(False)
        # autophase checkbox
        self.onlyT10Checkbox = QCheckBox(self)
        self.onlyT10Checkbox.setStyleSheet('font : bold 14px')
        self.onlyT10Checkbox.move(10,585)
        self.onlyT10Checkbox.resize(100,20)
        self.onlyT10Checkbox.setText('Only T1(0)')
        self.onlyT10Checkbox.setChecked(False)
        
        # Create a next button
        self.nextButton = QPushButton('Next Plot', self)
        self.nextButton.setStyleSheet('font : bold 14px; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.nextButton.move(625,525)
        self.nextButton.resize(100,40)
        
        # Create a next button
        self.autoButton = QPushButton('Auto Process', self)
        self.autoButton.setStyleSheet('font : bold 14px; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.autoButton.move(740,525)
        self.autoButton.resize(100,40)
        
        # Create a back button
        self.backButton = QPushButton('Back', self)
        self.backButton.setStyleSheet('font : bold 14px; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.backButton.move(5,525)
        self.backButton.resize(100,40)
        
        # T1(0) label
        self.errorLabel = QLabel(self)
        self.errorLabel.setStyleSheet('font : bold 14px')
        self.errorLabel.move(730,535)
        self.errorLabel.resize(200,20)
        #self.errorLabel.setFont(QtGui.QFont('Helvetica', 16))
        self.errorLabel.setText('')
        
        # T1(0) label
        self.t10Label = QLabel(self)
        self.t10Label.setStyleSheet('font : bold 14px')
        self.t10Label.move(345,525)
        self.t10Label.resize(80,20)
        self.t10Label.setText('T1(0) (s):')
        # Create a T10(0) text edit
        self.t10Edit = QLineEdit(self)
        self.t10Edit.move(410,525)
        self.t10Edit.resize(65,25)
        self.t10Edit.setText('2.5')
        
        # T1(0) label
        self.workupt10Label = QLabel(self)
        self.workupt10Label.setStyleSheet('font : bold 14px')
        self.workupt10Label.move(345,560)
        self.workupt10Label.resize(150,20)
        self.workupt10Label.setText('workup T1(0) (s):')
        # Create a T10(0) text edit
        self.workupt10Edit = QLineEdit(self)
        self.workupt10Edit.move(465,560)
        self.workupt10Edit.resize(65,25)
        self.workupt10Edit.setText('2.5')

        # T1 interpolation label
        self.t1fitLabel = QLabel(self)
        self.t1fitLabel.setStyleSheet('font : bold 14px')
        self.t1fitLabel.move(750,470)
        self.t1fitLabel.resize(230,20)
        self.t1fitLabel.setText('T1 interpolation:')
        # linear interpolation checkbox
        self.linearfitCheckbox = QCheckBox(self)
        self.linearfitCheckbox.setStyleSheet('font : bold 14px')
        self.linearfitCheckbox.move(865,470)
        self.linearfitCheckbox.resize(100,20)
        self.linearfitCheckbox.setText('Linear')
        # 2nd order interpolation checkbox
        self.order2fitCheckbox = QCheckBox(self)
        self.order2fitCheckbox.setStyleSheet('font : bold 14px')
        self.order2fitCheckbox.move(930,470)
        self.order2fitCheckbox.resize(100,20)
        self.order2fitCheckbox.setText('2nd Order')
        # Exclude first T1(p) checkbox
        self.exclude1T1Checkbox = QCheckBox(self)
        self.exclude1T1Checkbox.setStyleSheet('font : bold 14px')
        self.exclude1T1Checkbox.move(865,500)
        self.exclude1T1Checkbox.resize(150,20)
        self.exclude1T1Checkbox.setText('Exclude first T1(p)')
        self.exclude1T1Checkbox.setChecked(False)
        
        # T10(0) label
        self.t100Label = QLabel(self)
        self.t100Label.setStyleSheet('font : bold 14px')
        self.t100Label.move(183,525)
        self.t100Label.resize(80,20)
        self.t100Label.setText('T1<sub>0</sub>(0) (s):')
        # Create a T10(0) text edit
        self.t100Edit = QLineEdit(self)
        self.t100Edit.move(255,525)
        self.t100Edit.resize(65,25)
        self.t100Edit.setText('2.5')
        
        # concentration label
        self.slcLabel = QLabel(self)
        self.slcLabel.setStyleSheet('font : bold 14px')
        self.slcLabel.move(123,560)
        self.slcLabel.resize(150,20)
        self.slcLabel.setText('Concentration (uM):')
        # Create a concentration text edit
        self.slcEdit = QLineEdit(self)
        self.slcEdit.move(265,560)
        self.slcEdit.resize(65,25)
        self.slcEdit.setText('100')
        
        # magnetic field label
        self.fieldLabel = QLabel(self)
        self.fieldLabel.setStyleSheet('font : bold 14px')
        self.fieldLabel.move(123,595)
        self.fieldLabel.resize(150,20)
        self.fieldLabel.setText('Magnetic Field (mT):')
        # Create a magnetic field text edit
        self.fieldEdit = QLineEdit(self)
        self.fieldEdit.move(265,595)
        self.fieldEdit.resize(65,25)
        self.fieldEdit.setText('348.5')

        # smax label
        self.smaxLabel = QLabel(self)
        self.smaxLabel.setStyleSheet('font : bold 14px')
        self.smaxLabel.move(345,595)
        self.smaxLabel.resize(100,20)
        self.smaxLabel.setText('s<sub>max</sub> model:')
        # tethered checkbox
        self.tetheredCheckbox = QCheckBox(self)
        self.tetheredCheckbox.setStyleSheet('font : bold 14px')
        self.tetheredCheckbox.move(425,595)
        self.tetheredCheckbox.resize(100,20)
        self.tetheredCheckbox.setText('Tethered')
        # bulk checkbox
        self.bulkCheckbox = QCheckBox(self)
        self.bulkCheckbox.setStyleSheet('font : bold 14px')
        self.bulkCheckbox.move(510,595)
        self.bulkCheckbox.resize(100,20)
        self.bulkCheckbox.setText('Bulk-like')
        
        # MATLAB output button
        self.matoutButton = QPushButton('MATLAB output', self)
        self.matoutButton.setStyleSheet('font : bold 14px; color : rgb(255,184,20) ; background-color : rgb(0, 77, 159)')
        self.matoutButton.move(925,590)
        self.matoutButton.resize(120,30)
        
        self.guiDict = {'guiFunction' : {}, 'folderStruct' : {}, 'rawdataFunction' : {}, 'processSpectrum' : {}, 'workupFunction' : {},  'dnpLabFunction' : {}, 'workupData' : {}, 'dnpLabData' : {}, 'hydrationResults' : {}, 'dataPlot' : {}, 'enhancementPlot' : {}, 't1Plot' : {}}
        
        self.wrkup_smax = 1
        
        self.intwindowSlider.setMinimum(5)
        self.intwindowSlider.setMaximum(50)
        self.guiDict['processSpectrum']['integrationWidth'] = 20
        self.intwindowSlider.setValue(self.guiDict['processSpectrum']['integrationWidth'])
        
        self.guiDict['processSpectrum']['integrationCenter'] = 0
        self.intcenterSlider.setMinimum(self.guiDict['processSpectrum']['integrationCenter']-50)
        self.intcenterSlider.setMaximum(self.guiDict['processSpectrum']['integrationCenter']+50)
        self.intcenterSlider.setValue(self.guiDict['processSpectrum']['integrationCenter'])
        
        #set blank plots
        self.reset_plots()
        
        self.guiDict['guiFunction']['buttons'] = False
        self.guiDict['guiFunction']['sliders'] = False
        self.connect_widgets()
        
        self.show()

    def show_hide_components(self):
        
        if self.guiDict['guiFunction']['calculating']:
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
            self.bulkCheckbox.setVisible(True)
            self.matoutButton.setVisible(True)
            self.nextButton.setText('Hydration')
            self.autoButton.setVisible(False)
            self.backButton.setText('Restart')
            self.onlyT1pCheckbox.setVisible(True)
            self.onlyT10Checkbox.setVisible(True)
            self.intcenterLabel.setVisible(False)
            self.intcenterSlider.setVisible(False)
            self.intwindowLabel.setVisible(False)
            self.intwindowSlider.setVisible(False)
            self.phaseLabel.setVisible(False)
            self.phaseSlider.setVisible(False)
            self.autophaseCheckbox.setVisible(False)
            self.backButton.setVisible(True)

            if self.guiDict['workupFunction']['isWorkup'] or self.guiDict['dnpLabFunction']['isLab']:
                self.backButton.setVisible(False)
                self.onlyT1pCheckbox.setVisible(False)
                self.onlyT10Checkbox.setVisible(False)
                
            if self.guiDict['workupFunction']['isWorkup']:
                self.workupt10Label.setVisible(True)
                self.workupt10Edit.setVisible(True)
            else:
                self.t10Label.setVisible(True)
                self.t10Edit.setVisible(True)
                
            if self.guiDict['workupFunction']['addWorkup']:
                self.show_wrkupCheckbox.setVisible(True)
                self.fit_wrkupCheckbox.setVisible(True)
                self.guiDict['workupFunction']['show'] = True
                self.guiDict['workupFunction']['fit'] = False
                self.show_wrkupCheckbox.setChecked(True)
                self.fit_wrkupCheckbox.setChecked(False)
                
        else:
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
            self.bulkCheckbox.setVisible(False)
            self.matoutButton.setVisible(False)
            self.intcenterLabel.setVisible(True)
            self.intcenterSlider.setVisible(True)
            self.intwindowLabel.setVisible(True)
            self.intwindowSlider.setVisible(True)
            self.phaseLabel.setVisible(True)
            self.phaseSlider.setVisible(True)
            self.autophaseCheckbox.setVisible(True)
            self.show_wrkupCheckbox.setVisible(False)
            self.fit_wrkupCheckbox.setVisible(False)
            self.backButton.setText('Back')
            self.onlyT1pCheckbox.setVisible(False)
            self.onlyT10Checkbox.setVisible(False)
            self.nextButton.setText('Next Plot')
            self.autoButton.setVisible(True)
            self.backButton.setVisible(True)
            self.nextButton.setVisible(True)
        
        self.t1plt.setVisible(True)
        self.enhplt.setVisible(True)
        self.errorLabel.setVisible(False)

    def connect_widgets(self):
        
        self.intcenterSlider.valueChanged[int].connect(self.changeCenterValue)
        self.intwindowSlider.valueChanged[int].connect(self.changeWindowValue)
        self.phaseSlider.valueChanged[int].connect(self.changePhaseValue)
        self.autophaseCheckbox.clicked.connect(self.autoPhaseCheck)
        self.autophaseCheckbox.setChecked(False)
        self.nextButton.clicked.connect(self.on_click_next)
        self.autoButton.clicked.connect(self.autoProcess)
        self.backButton.clicked.connect(self.on_click_back)
        self.show_wrkupCheckbox.setChecked(True)
        self.show_wrkupCheckbox.clicked.connect(self.show_wrkupCheck)
        self.fit_wrkupCheckbox.setChecked(False)
        self.fit_wrkupCheckbox.clicked.connect(self.fit_wrkupCheck)
        self.t100Edit.editingFinished.connect(self.edithydText)
        self.t10Edit.editingFinished.connect(self.edithydText)
        self.workupt10Edit.editingFinished.connect(self.edithydText)
        self.linearfitCheckbox.clicked.connect(self.linearCheck)
        self.linearfitCheckbox.setChecked(False)
        self.order2fitCheckbox.clicked.connect(self.order2Check)
        self.order2fitCheckbox.setChecked(True)
        self.exclude1T1Checkbox.clicked.connect(self.exfT1pCheck)
        self.slcEdit.editingFinished.connect(self.edithydText)
        self.fieldEdit.editingFinished.connect(self.edithydText)
        self.tetheredCheckbox.clicked.connect(self.tetheredCheck)
        self.tetheredCheckbox.setChecked(True)
        self.bulkCheckbox.clicked.connect(self.bulkCheck)
        self.bulkCheckbox.setChecked(False)
        self.matoutButton.clicked.connect(self.matoutput)
    
    def reset_plots(self):
        
        self.guiDict['guiFunction']['autoProcess'] = False
        
        self.guiDict['dataPlot']['xdata'] = []
        self.guiDict['dataPlot']['ydata'] = []
        self.guiDict['dataPlot']['xmin'] = -1
        self.guiDict['dataPlot']['xmax'] = 1
        self.guiDict['dataPlot']['plotksig'] = False
        self.guiDict['dataPlot']['title'] = 'Spectrum'
        self.plot_data()
        
        self.guiDict['enhancementPlot']['xdata'] = []
        self.guiDict['enhancementPlot']['ydata'] = []
        self.guiDict['enhancementPlot']['xmin'] = 0
        self.guiDict['enhancementPlot']['xmax'] = 1
        self.guiDict['enhancementPlot']['title'] = 'E[p]'
        self.guiDict['enhancementPlot']['ytick'] = [0]
        self.guiDict['enhancementPlot']['ytickLabel'] = ['0']
        self.guiDict['enhancementPlot']['tau'] = []
        self.guiDict['enhancementPlot']['tauAmps'] = []
        self.guiDict['enhancementPlot']['t1Fit'] = []
        self.guiDict['enhancementPlot']['plotT1fit'] = False
        self.guiDict['enhancementPlot']['plotEpfit'] = False
        self.plot_enh()
        
        self.guiDict['t1Plot']['xdata'] = []
        self.guiDict['t1Plot']['ydata'] = []
        self.guiDict['t1Plot']['xmin'] = 0
        self.guiDict['t1Plot']['xmax'] = 1
        self.guiDict['t1Plot']['ymin'] = 0
        self.guiDict['t1Plot']['ymax'] = 4
        self.guiDict['t1Plot']['title'] = 'T1[p]'
        self.guiDict['t1Plot']['ytick'] = [1,3]
        self.guiDict['t1Plot']['ytickLabel'] = ['1','3']
        self.guiDict['t1Plot']['plotT1interp'] = False
        self.plot_t1()
        
        self.guiDict['guiFunction']['sliders'] = True
        self.guiDict['guiFunction']['hydrationEdits'] = False
        self.guiDict['guiFunction']['calculating'] = False
        
        self.guiDict['folderStruct']['index'] = 0
        self.guiDict['processSpectrum']['phaseFactor'] = 0

        self.show_hide_components()
    
    def plot_setter(self):
       
        if self.guiDict['rawdataFunction']['folder'] == -1:
            self.guiDict['dataPlot']['title'] = 'T1 Measurement, Folder # ' + self.singlefolder
            self.guiDict['enhancementPlot']['title'] = 'T1 Fit'
            
        elif self.guiDict['rawdataFunction']['folder'] == -2:
            self.guiDict['dataPlot']['title'] = '1D Data, Folder # ' + self.singlefolder
        else:
        
            if self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['p0']:
                self.guiDict['dataPlot']['title'] = 'Signal with power=0, Folder # ' + str(self.guiDict['rawdataFunction']['folder'])
                self.backButton.setText('Back')
                self.onlyT1pCheckbox.setVisible(False)
                self.onlyT10Checkbox.setVisible(False)
                
            elif self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['T10']:
                self.guiDict['dataPlot']['title'] = 'T1 with power=0, Folder # ' + str(self.guiDict['rawdataFunction']['folder'])
                
            elif self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['enh']:
                self.guiDict['dataPlot']['title'] = 'Enhanced Signal, Folder # ' + str(self.guiDict['rawdataFunction']['folder'])
                
            elif self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['T1']:
                self.guiDict['dataPlot']['title'] = 'T1 measurement, Folder # ' + str(self.guiDict['rawdataFunction']['folder'])
                
            if self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['T1'] or self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['T10']:
                self.guiDict['enhancementPlot']['title'] = 'T1 Fit'
                self.guiDict['enhancementPlot']['ytick'] = [0]
                self.guiDict['enhancementPlot']['ytickLabel'] = ['0']
                self.guiDict['enhancementPlot']['plotT1fit'] = True
                self.guiDict['enhancementPlot']['plotEpfit'] = False
                self.guiDict['enhancementPlot']['plotT1interp'] = False
            else:
                self.guiDict['enhancementPlot']['title'] = 'E[p]'
                self.guiDict['enhancementPlot']['plotT1fit'] = False
                self.guiDict['enhancementPlot']['plotEpfit'] = False
        
            self.guiDict['processSpectrum']['phaseFactor'] = 0
            self.guiDict['guiFunction']['sliders'] = False
            self.phaseSlider.setValue(0)
            self.guiDict['guiFunction']['sliders'] = True
    
    def loadLab(self):
        
        try:
            if self.testmode:
                flname = self.labpath + os.sep + 'data' + os.sep + 'topspin' + os.sep + 'GUI_output.mat'
            else:
                dir = QFileDialog.getOpenFileName(self)
                flname = dir[0]
            
            print('GUI Results: ' + flname)
            x = flname.split(os.sep)
            self.pathLabel.setText('GUI RESULTS DIRECTORY: ' + x[len(x)-2] + ' ' + os.sep + ' ' + x[len(x)-1])
            
            matin = loadmat(flname)

            self.ksiglabel='dnpLab'
            self.guiDict['rawdataFunction']['folder'] = -3
            
            self.guiDict['dnpLabFunction']['isLab'] = True
            self.guiDict['workupFunction']['isWorkup'] = False
            self.guiDict['workupFunction']['addWorkup'] = False
            self.guiDict['workupFunction']['show'] = False
            self.guiDict['workupFunction']['fit'] = False
            
            self.reset_plots()
            self.t10Edit.setText(str(round(float(matin['odnp']['T10']),4)))
            
            self.guiDict['dnpLabData']['T10'] = float(matin['odnp']['T10'])
            self.guiDict['dnpLabData']['T10_error'] = float(matin['odnp']['T10_error'])
            epows=matin['odnp']['Epowers'][0]
            self.guiDict['dnpLabData']['Epowers'] = np.ravel(epows[0])
            ep = matin['odnp']['Ep'][0]
            self.Ep = np.ravel(ep[0])
            t1pows = matin['odnp']['T1powers'][0]
            self.guiDict['dnpLabData']['T1powers'] = np.ravel(t1pows[0])
            t1p = matin['odnp']['T1p'][0]
            self.T1p = np.ravel(t1p[0])
            t1perr = matin['odnp']['T1p_error'][0]
            self.T1p_error = np.ravel(t1perr[0])
            
            self.guiDict['rawdataFunction']['nopowers'] = False
            
            self.finishProcessing()
            
            self.guiDict['guiFunction']['buttons'] = True
        
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText('GUI DATA ERROR')
            self.guiDict['guiFunction']['buttons'] = False
        
    def loadSingle(self):
        
        try:
            if self.testmode:
                dirname = dirname = self.labpath + os.sep + 'data' + os.sep + 'topspin' + os.sep + '304'
            else:
                dirname = QFileDialog.getExistingDirectory(self)
                
            dirname = os.path.join(dirname + os.sep)
            
            x = dirname.split(os.sep)
            self.pathLabel.setText('DATA DIRECTORY: ' + x[len(x)-3] + ' ' + os.sep + ' ' + x[len(x)-2])

            self.singlefolder = x[len(x)-2]
            path = dirname.replace(str(self.singlefolder) + os.sep, '')
            
            data = odnp.dnpImport.topspin.import_topspin(path,self.singlefolder)
            self.dnpLabWS = odnp.create_workspace('raw',data)
            self.dnpLabWS.copy('raw','proc')
            
            if self.dnpLabWS['proc'].ndim == 2:
                print('T1 Measurement: ' + dirname)
                self.guiDict['rawdataFunction']['folder'] = -1
            elif self.dnpLabWS['proc'].ndim == 1:
                print('1D Data: ' + dirname)
                self.guiDict['rawdataFunction']['folder'] = -2
                
            self.reset_plots()
            self.plot_setter()
            
            self.guiDict['guiFunction']['buttons'] = False
            self.guiDict['workupFunction']['isWorkup'] = False
            self.guiDict['workupFunction']['addWorkup'] = False
            self.guiDict['dnpLabFunction']['isLab'] = False
            self.guiDict['workupFunction']['fit'] = False
            self.guiDict['workupFunction']['show'] = False
            self.guiDict['enhancementPlot']['plotT1fit'] = True
            self.backButton.setVisible(False)
            self.onlyT1pCheckbox.setVisible(False)
            self.onlyT10Checkbox.setVisible(False)
            self.nextButton.setVisible(False)
            self.autoButton.setVisible(False)
            self.t1plt.setVisible(False)
            
            self.processData()
            
            if self.guiDict['rawdataFunction']['folder'] == -2:
                self.enhplt.setVisible(False)
            elif self.guiDict['rawdataFunction']['folder'] == -1:
                self.enhplt.setVisible(True)
   
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText('T1 DATA ERROR')
            self.guiDict['guiFunction']['sliders'] = False

    def loadWorkup(self):
    
        try:
            if self.testmode:
                wrkname = self.labpath + os.sep + 'data' + os.sep + 'topspin' + os.sep + 'Workup'
            else:
                wrkname = QFileDialog.getExistingDirectory(self)
                
            wrkname = os.path.join(wrkname + os.sep)
            self.guiDict['workupFunction']['directory'] = wrkname
            print('Workup: ' + wrkname)
            x = wrkname.split(os.sep)
            self.pathLabel.setText('WORKUP DIRECTORY: ' + x[len(x)-3] + ' ' + os.sep + ' ' + x[len(x)-2])
            
            self.ksiglabel='workup'
            self.guiDict['rawdataFunction']['folder'] = -3

            self.guiDict['workupFunction']['isWorkup'] = True
            self.guiDict['dnpLabFunction']['isLab'] = False
            self.guiDict['workupFunction']['addWorkup'] = False
            self.guiDict['workupFunction']['show'] = False
            self.guiDict['workupFunction']['fit'] = False

            self.reset_plots()
            self.processWorkup()
            self.guiDict['dnpLabData'] = copy.deepcopy(self.guiDict['workupData'])
            self.Ep = self.guiDict['workupData']['Ep']
            self.T1p = self.guiDict['workupData']['T1p']
            self.T1p_error = self.guiDict['workupData']['T1p_error']
            self.guiDict['dnpLabData']['T10'] = self.guiDict['workupData']['T10']
            self.guiDict['dnpLabData']['T10_error'] = self.guiDict['workupData']['T10_error']
            
            self.guiDict['rawdataFunction']['nopowers'] = False
            
            self.finishProcessing()

            self.guiDict['guiFunction']['buttons'] = True
            
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText('WORKUP ERROR')
            self.guiDict['guiFunction']['buttons'] = False
    
    def processWorkup(self):

        # load enhancementPowers.csv
        Etest = np.loadtxt(self.guiDict['workupFunction']['directory'] + 'enhancementPowers.csv', delimiter = ',', usecols = range(0,3), max_rows=22, skiprows = 1)
        if Etest[0,2] == 0:
            Eraw = Etest
        else:
            try:
                Eraw = np.loadtxt(self.guiDict['workupFunction']['directory'] + 'enhancementPowers.csv', delimiter = ',', usecols = range(0,2), max_rows=30, skiprows = 1)
            except IndexError:
                Eraw = np.loadtxt(self.guiDict['workupFunction']['directory'] + 'enhancementPowers.csv', delimiter = ',', usecols = range(0,2), max_rows=22, skiprows = 1)
        
        # load t1Powers.csv
        T1test = np.loadtxt(self.guiDict['workupFunction']['directory'] + 't1Powers.csv', delimiter = ',', usecols = range(0,4), max_rows=6, skiprows = 1)
        if T1test[5,3] == 304:
            T1raw = T1test
        else:
            T1test = np.loadtxt(self.guiDict['workupFunction']['directory'] + 't1Powers.csv', delimiter = ',', usecols = range(0,4), max_rows=9, skiprows = 1)
            if T1test[8,3] == 36:
                T1raw = T1test
            else:
                T1test = np.loadtxt(self.guiDict['workupFunction']['directory'] + 't1Powers.csv', delimiter = ',', usecols = range(0,4), max_rows=10, skiprows = 1)
                if T1test[9,3] == 37:
                    T1raw = T1test
                else:
                    T1test = np.loadtxt(self.guiDict['workupFunction']['directory'] + 't1Powers.csv', delimiter = ',', usecols = range(0,4), max_rows=11, skiprows = 1)
                    if T1test[10,3] == 304:
                        T1raw = T1test
    
        ePows = Eraw[:,0].reshape(-1)
        eP = Eraw[:,1].reshape(-1)
        self.guiDict['workupData']['Epowers'] = ePows[1:len(ePows)]
        # take real enhancement points
        self.guiDict['workupData']['Ep'] = eP[1:len(ePows)]

        t1Pows = T1raw[:,0].reshape(-1)
        t1P = T1raw[:,1].reshape(-1)
        t1P_error = T1raw[:,2].reshape(-1)
        self.guiDict['workupData']['T1powers'] = t1Pows[0:len(t1Pows)-1]
        self.guiDict['workupData']['T1p'] = t1P[0:len(t1P)-1]
        self.guiDict['workupData']['T1p_error'] = t1P_error[0:len(t1P_error)-1]
        self.guiDict['workupData']['T10'] = t1P[len(t1P)-1]
        self.guiDict['workupData']['T10_error'] = t1P_error[len(t1P_error)-1]
        self.workupt10Edit.setText(str(round(self.guiDict['workupData']['T10'],4)))
        
        wrkupksig = np.loadtxt(self.guiDict['workupFunction']['directory'] + 'kSigma.csv', delimiter = ',', usecols = range(0,2), max_rows=1, skiprows = 1)
        
        self.guiDict['workupData']['kSigma'] = wrkupksig[0] * 1e6
        self.guiDict['workupData']['kSigma_error'] = wrkupksig[1] * 1e6
        
        """
        wrkupksig_array = np.loadtxt(self.guiDict['workupFunction']['directory'] + 'kSigma.csv', delimiter = ',', usecols = range(0,1), max_rows=21, skiprows = 6)
        ksig = np.array([self.guiDict['workupData']['Epowers'],wrkupksig_array.reshape(-1)])
        ksig = np.transpose(ksig)
        ksig = ksig[ksig[:,0].argsort()]
        self.workup_ksig_array = ksig[:,1]
        """

    def workupExpTimes(self):#{{{
        """ This reads the bruker meta data file audita.txt and returns the time elapsed for the given experiment.

        Args:
        fullPath - (string) path to experiment file
        exps - (array) experiment numbers

        Kwargs:
        dnpExp - True = enhancements or False = T1s

        Returns:
        array of experiment times, and a list of the absolute experiment start and stop times

        """
        fullPath = self.guiDict['rawdataFunction']['directory']
        exps = self.ExpNums
        dnpExp = self.dnpExp
        
        expTime = []
        self.absTime = []
        for exp in exps:
            opened = open(fullPath + str(exp) + os.sep + 'audita.txt')
            lines = opened.readlines()
            absStart = lines[8].split(' ')[2] + ' ' + lines[8].split(' ')[3]
            splitup = re.findall(r"[\w']+",absStart)
            absStart = datetime.datetime(int(splitup[0]),int(splitup[1]),int(splitup[2]),int(splitup[3]),int(splitup[4]),int(splitup[5]),int(splitup[6]))
            absStart = time.mktime(absStart.utctimetuple()) # this returns seconds since the epoch
            start = lines[8].split(' ')[3]
            start = start.split(':') # hours,min,second
            hour = int(start[0],10)*3600
            minute = int(start[1],10)*60
            second = int(start[2].split('.')[0],10)
            start = second+minute+hour # in seconds
            absStop = lines[6].split('<')[1].split('>')[0].split(' ')
            absStop = absStop[0] + ' ' + absStop[1]
            splitup = re.findall(r"[\w']+",absStop)
            absStop = datetime.datetime(int(splitup[0]),int(splitup[1]),int(splitup[2]),int(splitup[3]),int(splitup[4]),int(splitup[5]),int(splitup[6]))
            absStop = time.mktime(absStop.utctimetuple()) # this returns seconds since the epoch
            stop = lines[6].split(' ')[4]
            stop = stop.split(':')
            hour = int(stop[0],10)*3600
            minute = int(stop[1],10)*60
            second = int(stop[2].split('.')[0],10)
            stop = second+minute+hour # in seconds
            expTime.append(stop-start)
            self.absTime.append((absStart,absStop))
        
        expTime = list(expTime)
        for count,timeVal in enumerate(expTime):
            if timeVal < 0:
                expTime.pop(count)
    
    def workupSplitPowers(self):
        """ Reads power file from odnp experiment and returns a list of the determined power steps.

        Args:
        fullPath - (string) 'path/to/expDirectory'
        powerfile - (string) 'fileName', do not include extension .mat or .csv
        absTime - (list of tuples) the absolute time values for the (start, stop) of each NMR experiment. This is returned from returnExpTimes in nmr.py
        threshold - (double) threshold value for the power steps. Units = (d dBm / d t) - differential

        Returns:
        powerList - (array) This is an array of the average power value measured between each start and stop time of the absTimes list.

        This code works by aligning the time values given in absTime to the last spike of the derivative powers spectrum. The last spike corresponds to when the amplifier is turned off and typically has the largest spike. Once the times are aligned I average the power between each start and stop time. I return the average of the power values.
        """
        fullPath = self.guiDict['rawdataFunction']['directory']
        powerfile = self.pwrfile
        bufferVal = self.buff
        threshold = 20
        
        if os.path.isfile(fullPath + powerfile + '.mat'): # This is a matlab file from cnsi
            print('Extracted powers from ' + self.pwrfile + '.mat file')
            openfile = loadmat(fullPath + os.sep + powerfile + '.mat')
            power = openfile.pop('powerlist')
            power = np.array([x for i in power for x in i])
            exptime = openfile.pop('timelist')
            exptime = np.array([x for i in exptime for x in i])
        elif os.path.isfile(fullPath + powerfile + '.csv'): # This is a csv file
            print('Extracted powers from ' + self.pwrfile + '.csv file')
            openfile = open(fullPath + os.sep + powerfile + '.csv','r')
            lines = openfile.readlines()
            if len(lines) == 1:
                lines = lines[0].split('\r') # this might not be what I want to do...
            lines.pop(0)
            timeList = []
            powerList = []
            for line in lines:
                exptime,power = line.split('\r')[0].split(',')
                timeList.append(float(exptime))
                powerList.append(float(power))
            exptime = np.array(timeList)
            power = np.array(powerList)

        #### Take the derivative of the power list
        step = exptime[1]-exptime[0]
        dp = []
        for i in range(len(power) - 1):
            dp.append((power[i+1] - power[i])/step)
        dp = abs(np.array(dp))
        ### Go through and threshold the powers
        timeBreak = []
        
        for i in range(len(dp)):
            if dp[i] >= threshold:
                timeBreak.append(exptime[i])
        timeBreak.sort()

        self.absTime.sort(key=lambda tup: tup[0])

        # align to the last spike
        offSet = self.absTime[-1][1] - timeBreak[-1] + bufferVal

        self.powerList = []
        for timeVals in self.absTime:
            start = int(timeVals[0] - offSet + bufferVal)
            stop = int(timeVals[1] - offSet - bufferVal)
            cutPower = []
            for k in range(0,len(exptime)-1):
                if start <= exptime[k] <= stop:
                    cutPower.append(power[k])
            powers = round(np.average(cutPower),3)
            self.powerList.append(float(powers))
    
    def workupPhaseOpt(self):
    
           optDict = copy.deepcopy(self.procDict)
           
           curve = optDict['proc'].values
           #{{{ find bestindex once
           phases = np.linspace(-np.pi/2,np.pi/2,100).reshape(1,-1) # this should work w/out altering the sign
           rotated_data = (curve.reshape(-1,1))*np.exp(-1j*phases)
           success = (np.real(rotated_data)**2).sum(axis=0)/((np.imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
           bestindex = np.argmax(success)
           #}}}
           #{{{ find the bestindex based on that
           if (bestindex>phases.shape[1]-2):
               bestindex = phases.shape[1]-2
               phases = np.linspace(
                   phases[0,bestindex-1],
                   phases[0,bestindex+1],
                   100).reshape(1,-1)
               rotated_data = (curve.reshape(-1,1))*np.exp(-1j*phases)
               success = (np.real(rotated_data)**2).sum(axis=0)/((np.imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
               bestindex = np.argmax(success)
               
           self.guiDict['processSpectrum']['originalPhase'] = phases[0,bestindex]
             
    def loadRaw(self):
    
        try:
            if self.testmode:
                dirname = self.labpath + os.sep + 'data' + os.sep + 'topspin'
            else:
                dirname = QFileDialog.getExistingDirectory(self)
                
            dirname = os.path.join(dirname + os.sep)
            self.guiDict['rawdataFunction']['directory'] = dirname
            print('Data: ' + dirname)
            x = dirname.split(os.sep)
            self.pathLabel.setText('DATA DIRECTORY: ' + x[len(x)-3] + ' ' + os.sep + ' ' + x[len(x)-2])

            self.guiDict['folderStruct'] = {}
            self.guiDict['workupFunction']['isWorkup'] = False
            self.guiDict['dnpLabFunction']['isLab'] = False

            if os.path.exists(self.guiDict['rawdataFunction']['directory'] + '40'):

                self.guiDict['folderStruct']['p0'] = 5
                self.guiDict['folderStruct']['enh'] = list(range(6,30))
                self.guiDict['folderStruct']['T1'] = range(31,41)
                self.guiDict['folderStruct']['T10'] = 304

            else:
                
                self.guiDict['folderStruct']['p0'] = 5
                self.guiDict['folderStruct']['enh'] = range(6,27)
                self.guiDict['folderStruct']['T1'] = range(28,33)
                self.guiDict['folderStruct']['T10'] = 304
                                
            self.guiDict['workupFunction']['addWorkup'] = False
            self.guiDict['workupFunction']['show'] = False
            self.guiDict['workupFunction']['fit'] = False
            
            self.guiDict['rawdataFunction']['nopowers'] = True
            
            if os.path.exists(self.guiDict['rawdataFunction']['directory'] + 'Workup') and os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 'Workup' + os.sep + 'enhancementPowers.csv') and os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 'Workup' + os.sep + 'kSigma.csv') and os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 'Workup' + os.sep + 't1Powers.csv'):
                
                self.guiDict['workupFunction']['addWorkup'] = True
                self.guiDict['workupFunction']['show'] = True
                
                self.guiDict['workupFunction']['directory'] = os.path.join(dirname + 'Workup' + os.sep)
                
                self.processWorkup()
                
                if len(self.guiDict['workupData']['Epowers']) == len(self.guiDict['folderStruct']['enh']) and len(self.guiDict['workupData']['T1powers']) == len(self.guiDict['folderStruct']['T1']):
                    
                    Epowers = self.guiDict['workupData']['Epowers']
                    T1powers = self.guiDict['workupData']['T1powers']
                
                    print('Found workup output, using power values from workup.')
                    self.guiDict['rawdataFunction']['nopowers'] = False
            
            if self.guiDict['rawdataFunction']['nopowers']:
                if os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 'power.mat') or os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 'power.csv'):
                
                    if os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 't1_powers.mat') or os.path.isfile(self.guiDict['rawdataFunction']['directory'] + 't1_powers.csv'):
                        
                        print('No workup output found, using power readings files.')
                        
                        self.ExpNums = self.guiDict['folderStruct']['enh']
                        self.dnpExp = True
                        self.workupExpTimes()
                        self.pwrfile = 'power'
                        self.buff = 2.5
                        self.workupSplitPowers()
                        
                        Epowers = np.add(self.powerList, 21.9992)
                        Epowers = np.divide(Epowers, 10)
                        Epowers = np.power(10, Epowers)
                        Epowers = np.multiply((1e-3),Epowers)
                        
                        self.ExpNums = self.guiDict['folderStruct']['T1']
                        self.dnpExp = False
                        self.workupExpTimes()
                        self.pwrfile = 't1_powers'
                        self.buff = 20*2.5
                        self.workupSplitPowers()
                        
                        T1powers = np.add(self.powerList, 21.9992)
                        T1powers = np.divide(T1powers, 10)
                        T1powers = np.power(10, T1powers)
                        T1powers = np.multiply((1e-3),T1powers)
                        
                        self.guiDict['rawdataFunction']['nopowers'] = False
                          
            if self.guiDict['rawdataFunction']['nopowers']:
                print('No power readings found.')
                print('Trying to find power settings in experiment titles...')
                try:
                    Eplist = []
                    for k in self.guiDict['folderStruct']['enh']:
                        title = odnp.dnpImport.bruker.loadTitle(self.guiDict['rawdataFunction']['directory'],expNum = k)
                        splitTitle = title.split(' ')
                        Eplist.append(float(splitTitle[-1]))
                        
                    T1plist = []
                    for k in self.guiDict['folderStruct']['T1']:
                        title = odnp.dnpImport.bruker.loadTitle(self.guiDict['rawdataFunction']['directory'],expNum = k)
                        splitTitle = title.split(' ')
                        T1plist.append(float(splitTitle[-1]))

                    Epowers = np.multiply(-1,Eplist)
                    Epowers = np.add(Epowers, 29.01525)
                    Epowers = np.divide(Epowers, 10)
                    Epowers = np.power(10, Epowers)
                    Epowers = np.multiply((1e-3),Epowers)

                    T1powers = np.multiply(-1,T1plist)
                    T1powers = np.add(T1powers, 29.01525)
                    T1powers = np.divide(T1powers, 10)
                    T1powers = np.power(10, T1powers)
                    T1powers = np.multiply((1e-3),T1powers)
                    
                    print('Powers taken from experiment titles. *WARNING: this is not accurate!')
                    self.guiDict['rawdataFunction']['nopowers'] = False
                    
                except:
                    print('No power readings available. E[p] and T1[p] are indexed by folder #.')
                    Epowers = self.guiDict['folderStruct']['enh']
                    T1powers = self.guiDict['folderStruct']['T1']
            
            self.guiDict['folderStruct']['all'] = []
            self.guiDict['folderStruct']['all'].append(self.guiDict['folderStruct']['p0'])
            for k in self.guiDict['folderStruct']['enh']:
                self.guiDict['folderStruct']['all'].append(k)
            for k in self.guiDict['folderStruct']['T1']:
                self.guiDict['folderStruct']['all'].append(k)
            self.guiDict['folderStruct']['all'].append(self.guiDict['folderStruct']['T10'])
            
            self.Ep = []
            self.T1p = []
            self.T1p_error = []
            self.guiDict['dnpLabData']['Epowers'] = Epowers
            self.guiDict['dnpLabData']['T1powers'] = T1powers
            self.originalEPowers = self.guiDict['dnpLabData']['Epowers']
            self.originalT1Powers = self.guiDict['dnpLabData']['T1powers']
            self.guiDict['guiFunction']['buttons'] = True
            
            self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['p0']
            self.ksiglabel='dnpLab'
            
            self.reset_plots()
            self.plot_setter()
            
            data = odnp.dnpImport.topspin.import_topspin(self.guiDict['rawdataFunction']['directory'],self.guiDict['rawdataFunction']['folder'])
            self.dnpLabWS = odnp.create_workspace('raw',data)
            self.dnpLabWS.copy('raw','proc')

            self.processData()
            
        except:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            self.pathLabel.setText('RAW DATA ERROR')
            self.guiDict['guiFunction']['buttons'] = False
            self.guiDict['guiFunction']['sliders'] = False

    def on_click_next(self):
        
        if self.guiDict['guiFunction']['buttons']:
        
            if self.guiDict['workupFunction']['isWorkup']:
                self.run_hydration()
            else:
            
                if self.guiDict['rawdataFunction']['folder'] == -3 or self.guiDict['folderStruct']['index'] >= len(self.guiDict['folderStruct']['all']):
                    self.run_hydration()
                else:
                
                    self.nextProcessing()
        else:
            pass
            
    def t1Func(self, tau, x1, x2, x3):
        return x2 - x3 * np.exp(-1.*tau/x1)
            
    def nextProcessing(self):
        
        nextprocDict = copy.deepcopy(self.procDict)
        
        nextprocDict['proc'] *= np.exp(-1j*self.guiDict['processSpectrum']['phase'])
        nextprocDict = odnp.dnpNMR.integrate(nextprocDict,{'integrate_center' : self.guiDict['processSpectrum']['integrationCenter'] , 'integrate_width' : self.guiDict['processSpectrum']['integrationWidth']})
            
        if self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['p0']:
            self.guiDict['dnpLabData']['p0'] = np.real(nextprocDict['proc'].values[0])
        elif self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['enh']:
            Ep = np.real(nextprocDict['proc'].values[0]) / self.guiDict['dnpLabData']['p0']
            self.Ep.append(np.real(Ep))
            if self.guiDict['guiFunction']['autoProcess']:
                pass
            else:
                self.guiDict['enhancementPlot']['xdata'] = self.guiDict['dnpLabData']['Epowers'][0:len(self.Ep)]
                self.guiDict['enhancementPlot']['ydata'] = self.Ep
                self.guiDict['enhancementPlot']['ytick'] = [0,min(self.Ep)]
                if min(self.Ep) <= -10:
                    self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(int(min(self.Ep)))]
                else:
                    self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(round(min(self.Ep),1))]
                self.plot_enh()
        elif self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['T1'] or self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['T10']:
        
            tau = np.reshape(nextprocDict['proc'].coords,-1)
            Mdata = np.real(nextprocDict['proc'].values)
            popt, pcov = optimize.curve_fit(self.t1Func, tau, Mdata, p0=[1., Mdata[-1], Mdata[-1]], method='lm')
            stdd = np.sqrt(np.diag(pcov))

            if self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['T1']:
                self.T1p.append(popt[0])
                self.T1p_error.append(stdd[0])
            elif self.guiDict['rawdataFunction']['folder'] == self.guiDict['folderStruct']['T10']:
                self.guiDict['dnpLabData']['T10'] = popt[0]
                self.guiDict['dnpLabData']['T10_error'] = stdd[0]
                self.t10Edit.setText(str(round(self.guiDict['dnpLabData']['T10'],4)))

            if self.guiDict['guiFunction']['autoProcess']:
                pass
            else:
                new_tau = np.r_[np.min(self.guiDict['t1Plot']['tau']):np.max(self.guiDict['t1Plot']['tau']):100j]
                self.guiDict['t1Plot']['t1Fit'] = self.t1Func(new_tau, popt[0],popt[1],popt[2])
                self.guiDict['t1Plot']['t1Val'] = popt[0]
            
                self.guiDict['t1Plot']['xdata'] = self.guiDict['dnpLabData']['T1powers'][0:len(self.T1p)]
                self.guiDict['t1Plot']['ydata'] = self.T1p
                self.guiDict['t1Plot']['ymin'] = min(self.guiDict['t1Plot']['ydata'])*.9
                self.guiDict['t1Plot']['ymax'] = max(self.guiDict['t1Plot']['ydata'])*1.1
                self.guiDict['t1Plot']['ytick'] = [max(self.T1p)]
                self.guiDict['t1Plot']['ytickLabel'] = [str(round(max(self.T1p),1))]
                
                self.guiDict['t1Plot']['tau'] = tau
                self.guiDict['t1Plot']['t1Amps'] = Mdata
                self.plot_t1()
                self.plot_enh()
        
        self.guiDict['folderStruct']['index'] += 1
        if self.guiDict['guiFunction']['autoProcess']:
            print('Finished with Folder #' + str(self.guiDict['folderStruct']['index']) + ' of ' + str(len(self.guiDict['folderStruct']['all'])))
        
        if self.guiDict['folderStruct']['index'] >= len(self.guiDict['folderStruct']['all']):
            self.finishProcessing()
        else:
            self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['all'][self.guiDict['folderStruct']['index']]
            
            if self.guiDict['guiFunction']['autoProcess']:
                pass
            else:
                if self.guiDict['folderStruct']['index'] == len(self.guiDict['folderStruct']['all'])-1:
                    self.nextButton.setText('Finish')
                self.plot_setter()
                
            data = odnp.dnpImport.topspin.import_topspin(self.guiDict['rawdataFunction']['directory'],self.guiDict['rawdataFunction']['folder'])
            self.dnpLabWS = odnp.create_workspace('raw',data)
            self.dnpLabWS.copy('raw','proc')
            
            self.processData()
            
    def on_click_back(self):
           
        if self.guiDict['guiFunction']['buttons']:
           
           self.guiDict['folderStruct']['index'] -= 1

           if self.guiDict['folderStruct']['index'] <= 0:
               self.guiDict['folderStruct']['index'] = 0
               self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['p0']

           if self.guiDict['folderStruct']['index'] >= len(self.guiDict['folderStruct']['all'])-1:
              self.reset_plots()
              self.guiDict['dnpLabData']['Epowers'] = self.originalEPowers
              self.guiDict['dnpLabData']['T1powers'] = self.originalT1Powers
              if self.onlyT10Checkbox.isChecked():
                 self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['T10']
                 self.guiDict['folderStruct']['index'] = len(self.guiDict['folderStruct']['all'])-1
                 self.nextButton.setText('Finish')
              else:
                  if self.onlyT1pCheckbox.isChecked():
                     self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['T1'][0]
                     self.guiDict['folderStruct']['index'] = len(self.guiDict['folderStruct']['all'])-1-len(self.guiDict['folderStruct']['T1'])
                     self.T1p = []
                     self.T1p_error = []
                  else:
                     self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['p0']
                     self.Ep = []
                     self.T1p = []
                     self.T1p_error = []
           else:
                self.guiDict['rawdataFunction']['folder'] = self.guiDict['folderStruct']['all'][self.guiDict['folderStruct']['index']]
                
           if self.guiDict['folderStruct']['index'] == len(self.guiDict['folderStruct']['all'])-2:
                self.nextButton.setText('Next Plot')
                
           self.plot_setter()
           if self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['enh']:
               if len(self.Ep) < 2:
                   self.Ep =[]
                   self.guiDict['enhancementPlot']['xdata'] = []
                   self.guiDict['enhancementPlot']['ydata'] = []
               else:
                    self.Ep =self.Ep[0:len(self.Ep)-1]
                    self.guiDict['enhancementPlot']['xdata'] = self.guiDict['dnpLabData']['Epowers'][0:len(self.Ep)]
                    self.guiDict['enhancementPlot']['ydata'] = self.Ep
                    self.guiDict['enhancementPlot']['ytick'] = [0,min(self.Ep)]
                    if min(self.Ep) <= -10:
                       self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(int(min(self.Ep)))]
                    else:
                       self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(round(min(self.Ep),1))]
               self.plot_enh()
                   
           elif self.guiDict['rawdataFunction']['folder'] in self.guiDict['folderStruct']['T1']:
               if len(self.T1p) < 2:
                   self.T1p = []
                   self.T1p_error = []
                   self.guiDict['t1Plot']['xdata'] = []
                   self.guiDict['t1Plot']['ydata'] = []
               else:
                    self.T1p = self.T1p[0:len(self.T1p)-1]
                    self.T1p_error = self.T1p_error[0:len(self.T1p_error)-1]
                    self.guiDict['t1Plot']['xdata']  = self.guiDict['dnpLabData']['T1powers'][0:len(self.T1p)]
                    self.guiDict['t1Plot']['ydata'] = self.T1p
                    self.guiDict['t1Plot']['ytick'] = [max(self.T1p)]
                    self.guiDict['t1Plot']['ytickLabel'] = [str(round(max(self.T1p),1))]
                    self.guiDict['t1Plot']['ymin']  = min(self.T1p)*.85
                    self.guiDict['t1Plot']['ymax'] = max(self.T1p)*1.15
               self.plot_t1()
               
           data = odnp.dnpImport.topspin.import_topspin(self.guiDict['rawdataFunction']['directory'],self.guiDict['rawdataFunction']['folder'])
           self.dnpLabWS = odnp.create_workspace('raw',data)
           self.dnpLabWS.copy('raw','proc')
           
           self.processData()
               
        else:
            pass
    
    def autoProcess(self):
    
        if self.guiDict['guiFunction']['buttons']:
            try:
                print('Auto processing, please wait...')
                self.guiDict['guiFunction']['autoProcess'] = True
                #t = time.time()
                for k in range(self.guiDict['folderStruct']['index']+1,len(self.guiDict['folderStruct']['all'])+1):
                    self.nextProcessing()
                #elapsed = time.time() - t
                #print('AutoProcess Time = ' + str(elapsed))
            except:
                print('Error in auto processing, resetting to folder # ' + str(self.guiDict['folderStruct']['p0']))
                self.guiDict['folderStruct']['index'] = len(self.guiDict['folderStruct']['all'])
                self.on_click_back()
        else:
            pass
        
    def processData(self):
        
        self.procDict = copy.deepcopy(self.dnpLabWS)
        
        self.procDict = odnp.dnpNMR.remove_offset(self.procDict,{})
        self.procDict = odnp.dnpNMR.window(self.procDict,{})
        self.procDict = odnp.dnpNMR.fourier_transform(self.procDict,{})
        phsDict = copy.deepcopy(self.procDict)
        self.guiDict['processSpectrum']['originalPhase'] = phsDict['proc'].phase()
        if self.autophaseCheckbox.isChecked() or self.guiDict['guiFunction']['autoProcess']:
            self.guiDict['processSpectrum']['integrationWidth'] = 15
            self.workupPhaseOpt()
        self.optCenter()
        
        if self.guiDict['guiFunction']['autoProcess']:
            self.guiDict['processSpectrum']['phase'] = self.guiDict['processSpectrum']['originalPhase']
        else:
            self.guiDict['guiFunction']['sliders'] = False
            
            fac = (np.pi/self.guiDict['processSpectrum']['originalPhase'])
            self.phaseSlider.setMinimum(round(-1000*abs(fac)))
            self.phaseSlider.setMaximum(round(1000*abs(fac)))
            self.phaseSlider.setValue(self.guiDict['processSpectrum']['originalPhase'])
            
            self.intcenterSlider.setValue(self.guiDict['processSpectrum']['integrationCenter'])
            self.intcenterSlider.setMinimum(self.guiDict['processSpectrum']['integrationCenter']-50)
            self.intcenterSlider.setMaximum(self.guiDict['processSpectrum']['integrationCenter']+50)
            
            self.intwindowSlider.setValue(self.guiDict['processSpectrum']['integrationWidth'])
            
            self.guiDict['guiFunction']['sliders'] = True
        
        xdata = self.procDict['proc'].coords
        self.guiDict['dataPlot']['xdata'] = np.reshape(xdata['t2'], -1)
        
        self.adjustSliders()
            
    def optCenter(self):
    
        optCentDict = copy.deepcopy(self.procDict)
        
        intgrlarray = []
        indx = range(-50,51)
        optCentDict['proc'] *= np.exp(-1j*self.guiDict['processSpectrum']['originalPhase'])
        for k in indx:
            iteroptCentDict = copy.deepcopy(optCentDict)
            iteroptCentDict = odnp.dnpNMR.integrate(iteroptCentDict, {'integrate_center' : k , 'integrate_width' : 10})
            iteroptCentDict['proc'].values = np.real(iteroptCentDict['proc'].values)
            intgrlarray.append(abs(iteroptCentDict['proc'].values[0]))
        cent = np.argmax(intgrlarray)
        self.guiDict['processSpectrum']['integrationCenter'] = indx[cent]
        
    def adjustSliders(self):
        
        adjsliderDict = copy.deepcopy(self.procDict)
        
        self.guiDict['processSpectrum']['phase'] =  self.guiDict['processSpectrum']['originalPhase'] + (self.guiDict['processSpectrum']['phaseFactor'] * self.guiDict['processSpectrum']['originalPhase'])

        ydata = adjsliderDict['proc'].values * np.exp(-1j*self.guiDict['processSpectrum']['phase'])
        self.guiDict['dataPlot']['ydata'] = np.real(ydata)
        
        xdata = adjsliderDict['proc'].coords
        self.guiDict['dataPlot']['xdata'] = np.reshape(xdata['t2'], -1)
        adjsliderDict['proc'] *= np.exp(-1j*self.guiDict['processSpectrum']['phase'])
        adjsliderDict = odnp.dnpNMR.integrate(adjsliderDict,{'integrate_center' : self.guiDict['processSpectrum']['integrationCenter'] , 'integrate_width' : self.guiDict['processSpectrum']['integrationWidth']})
        adjsliderDict['proc'].values = np.real(adjsliderDict['proc'].values)
        
        if len(adjsliderDict['proc'].values) == 1:
            pass
        else:
            
            self.guiDict['t1Plot']['t1Amps'] = adjsliderDict['proc'].values
            self.guiDict['t1Plot']['tau'] = np.reshape(adjsliderDict['proc'].coords,-1)
            
            popt, pcov = optimize.curve_fit(self.t1Func, self.guiDict['t1Plot']['tau'], self.guiDict['t1Plot']['t1Amps'], p0=[1., self.guiDict['t1Plot']['t1Amps'][-1], self.guiDict['t1Plot']['t1Amps'][-1]], method='lm')
            stdd = np.sqrt(np.diag(pcov))
            
            tau = np.r_[np.min(self.guiDict['t1Plot']['tau']):np.max(self.guiDict['t1Plot']['tau']):100j]
            self.guiDict['t1Plot']['t1Fit'] = self.t1Func(tau, popt[0],popt[1],popt[2])
            self.guiDict['t1Plot']['t1Val'] = popt[0]
            self.plot_enh()
            
            if self.guiDict['rawdataFunction']['folder'] == -1:

                print('---Error in T1---')
                print('T1: ' + str(round(popt[0],4)) + ' +/- ' + str(round(stdd[0],4)))
        
        self.guiDict['dataPlot']['xmin'] = int(round(self.guiDict['processSpectrum']['integrationCenter'] - np.abs(self.guiDict['processSpectrum']['integrationWidth'])/2))
        self.guiDict['dataPlot']['xmax'] = int(round(self.guiDict['processSpectrum']['integrationCenter'] + np.abs(self.guiDict['processSpectrum']['integrationWidth'])/2))
        self.plot_data()

    def finishProcessing(self):
        
        self.guiDict['guiFunction']['calculating'] = True
        self.show_hide_components()
        
        if self.guiDict['rawdataFunction']['nopowers']:
            self.guiDict['dnpLabData']['Ep'] = self.Ep
            self.guiDict['dnpLabData']['T1p'] = self.T1p
            self.guiDict['dnpLabData']['T1p_error'] = self.T1p_error
        else:
            enh = np.array([self.guiDict['dnpLabData']['Epowers'], self.Ep])
            enh = np.transpose(enh)
            enh = enh[enh[:,0].argsort()]
            self.guiDict['dnpLabData']['Epowers'] = enh[:,0]
            self.guiDict['dnpLabData']['Ep'] = enh[:,1]
            
            t1 = np.array([self.guiDict['dnpLabData']['T1powers'], self.T1p, self.T1p_error])
            t1 = np.transpose(t1)
            t1 = t1[t1[:,0].argsort()]
            self.guiDict['dnpLabData']['T1powers'] = t1[:,0]
            self.guiDict['dnpLabData']['T1p'] = t1[:,1]
            self.guiDict['dnpLabData']['T1p_error'] = t1[:,2]
        
        if self.guiDict['workupFunction']['isWorkup'] or self.guiDict['workupFunction']['addWorkup']:

            wenh = np.array([self.guiDict['workupData']['Epowers'], self.guiDict['workupData']['Ep']])
            wenh = np.transpose(wenh)
            wenh = wenh[wenh[:,0].argsort()]
            self.guiDict['workupData']['Epowers'] = wenh[:,0]
            self.guiDict['workupData']['Ep'] = wenh[:,1]
            
            wt1 = np.array([self.guiDict['workupData']['T1powers'], self.guiDict['workupData']['T1p'], self.guiDict['workupData']['T1p_error']])
            wt1 = np.transpose(wt1)
            wt1 = wt1[wt1[:,0].argsort()]
            self.guiDict['workupData']['T1powers'] = wt1[:,0]
            self.guiDict['workupData']['T1p'] = wt1[:,1]
            self.guiDict['workupData']['T1p_error'] = wt1[:,2]
            
            if self.guiDict['workupFunction']['addWorkup']:
                self.show_wrkupCheckbox.setVisible(True)
                self.fit_wrkupCheckbox.setVisible(True)

        self.guiDict['enhancementPlot']['xdata'] = self.guiDict['dnpLabData']['Epowers']
        self.guiDict['enhancementPlot']['ydata'] = self.guiDict['dnpLabData']['Ep']
        self.guiDict['enhancementPlot']['title'] = 'E[p]'
        self.guiDict['enhancementPlot']['ytick'] = [0,min(self.guiDict['dnpLabData']['Ep'])]

        if min(self.guiDict['dnpLabData']['Ep']) <= -10:
            self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(int(min(self.guiDict['dnpLabData']['Ep'])))]
        else:
            self.guiDict['enhancementPlot']['ytickLabel'] = ['0',str(round(min(self.guiDict['dnpLabData']['Ep']),1))]
        
        self.guiDict['t1Plot']['xdata'] = self.guiDict['dnpLabData']['T1powers']
        self.guiDict['t1Plot']['ydata'] = self.guiDict['dnpLabData']['T1p']
        self.guiDict['t1Plot']['title'] = 'T1[p]'
        self.guiDict['guiFunction']['sliders'] = False
        
        self.guiDict['dataPlot']['plotksig'] = True
        self.guiDict['t1Plot']['plotT1interp'] = True
        self.guiDict['enhancementPlot']['plotT1fit'] = False
        self.guiDict['enhancementPlot']['plotEpfit'] = True
        
        if self.guiDict['workupFunction']['isWorkup']:
            print('---workup Errors in T1---')
            print('T10: ' + str(round(self.guiDict['workupData']['T10'],2)) + ' +/- ' + str(round(self.guiDict['workupData']['T10_error'],4)))
            for k in range(0,len(self.guiDict['workupData']['T1p'])):
                print(str(round(self.guiDict['workupData']['T1p'][k],2)) + ' +/- ' + str(round(self.guiDict['workupData']['T1p_error'][k],4)))
        else:
            print('---Errors in T1---')
            print('T10: ' + str(round(self.guiDict['dnpLabData']['T10'],2)) + ' +/- ' + str(round(self.guiDict['dnpLabData']['T10_error'],4)))
            for k in range(0,len(self.T1p)):
                print(str(round(self.T1p[k],2)) + ' +/- ' + str(round(self.T1p_error[k],4)))
            
            if self.guiDict['workupFunction']['addWorkup']:
                print('---workup Errors in T1---')
                print('T10: ' + str(round(self.guiDict['workupData']['T10'],2)) + ' +/- ' + str(round(self.guiDict['workupData']['T10_error'],4)))
                for k in range(0,len(self.guiDict['workupData']['T1p'])):
                    print(str(round(self.guiDict['workupData']['T1p'][k],2)) + ' +/- ' + str(round(self.guiDict['workupData']['T1p_error'][k],4)))
                
        self.run_hydration()

    def run_hydration(self):
        self.err = False
        self.errorLabel.setVisible(False)
        
        spin_C = float(self.slcEdit.text())
        field = float(self.fieldEdit.text())
        T100 = float(self.t100Edit.text())
        self.guiDict['dnpLabData']['T100'] = float(self.t100Edit.text())
        if self.guiDict['workupFunction']['isWorkup']:
            T10 = float(self.workupt10Edit.text())
        else:
            T10 = float(self.t10Edit.text())
            
        if self.tetheredCheckbox.isChecked():
            smax_model = 'tethered'
            self.wrkup_smax = 1
        else:
            smax_model = 'free'
            self.wrkup_smax = 1-(2/(3+(3*(spin_C*1e-6*198.7))))
            
        if self.linearfitCheckbox.isChecked():
            t1_interp_method = 'linear'
        else:
            t1_interp_method = 'second_order'

        if self.exclude1T1Checkbox.isChecked():
        
            T1p = self.guiDict['dnpLabData']['T1p'][1:len(self.guiDict['dnpLabData']['T1p'])]
            T1powers = self.guiDict['dnpLabData']['T1powers'][1:len(self.guiDict['dnpLabData']['T1powers'])]
        else:
            
            T1p = self.guiDict['dnpLabData']['T1p']
            T1powers = self.guiDict['dnpLabData']['T1powers']
        
        
        hydration = {'E' : np.array(self.guiDict['dnpLabData']['Ep']), 'E_power' : np.array(self.guiDict['dnpLabData']['Epowers']), 'T1' : np.array(T1p), 'T1_power' : np.array(T1powers)}
        hydration.update({
            'T10': T10,
            'T100': self.guiDict['dnpLabData']['T100'],
            'spin_C': spin_C,
            'field': field,
            'smax_model': smax_model,
            't1_interp_method': t1_interp_method
        })
        hyd = odnp.create_workspace()
        hyd.add('hydration', hydration)

        try:
            hydresults = odnp.dnpHydration.hydration(hyd)
            self.guiDict['hydrationResults'] = copy.deepcopy(hydresults)
            self.plothydresults = copy.deepcopy(hydresults)
        except:
            self.err = True
            self.errorLabel.setVisible(True)
            if self.guiDict['workupFunction']['isWorkup']:
                self.errorLabel.setText('workup data Error')
            else:
                self.errorLabel.setText('dnpLab data Error')

        if self.guiDict['workupFunction']['isWorkup']:
            self.workupt10Label.setVisible(True)
            self.workupt10Edit.setVisible(True)
            self.t10Label.setVisible(False)
            self.t10Edit.setVisible(False)
        else:
            self.workupt10Label.setVisible(False)
            self.workupt10Edit.setVisible(False)
            self.t10Label.setVisible(True)
            self.t10Edit.setVisible(True)
            
        if self.guiDict['workupFunction']['isWorkup'] or self.guiDict['workupFunction']['addWorkup']:
        
                if self.guiDict['workupFunction']['fit'] or self.guiDict['workupFunction']['show']:
                    
                    if self.exclude1T1Checkbox.isChecked():
                    
                        wT1p = self.guiDict['workupData']['T1p'][1:len(self.guiDict['workupData']['T1p'])]
                        wT1powers = self.guiDict['workupData']['T1powers'][1:len(self.guiDict['workupData']['T1powers'])]
                
                    else:
                    
                        wT1p = self.guiDict['workupData']['T1p']
                        wT1powers = self.guiDict['workupData']['T1powers']
                    
                    wT10 = float(self.workupt10Edit.text())
                    
                    whydration = {'E' : np.array(self.guiDict['workupData']['Ep']), 'E_power' : np.array(self.guiDict['workupData']['Epowers']), 'T1' : np.array(wT1p), 'T1_power' : np.array(wT1powers)}
                    whydration.update({
                        'T10': wT10,
                        'T100': self.guiDict['dnpLabData']['T100'],
                        'spin_C': spin_C,
                        'field': field,
                        'smax_model': smax_model,
                        't1_interp_method': t1_interp_method
                    })
                    whyd = odnp.create_workspace()
                    whyd.add('hydration', whydration)
                    
                    try:
                        whydresults = odnp.dnpHydration.hydration(whyd)
                        self.guiDict['whydrationResults'] = copy.deepcopy(whydresults)
                        if self.guiDict['workupFunction']['fit']:
                            self.plothydresults = copy.deepcopy(whydresults)
                    except:
                        self.err = True
                        self.errorLabel.setVisible(True)
                        self.errorLabel.setText('workup data Error')
                    
                    self.workupt10Label.setVisible(True)
                    self.workupt10Edit.setVisible(True)

        self.guiDict['guiFunction']['hydrationEdits'] = True
        
        if self.err:
            self.dataplt.axes.cla()
            self.dataplt.draw()
            pass
        else:

            self.guiDict['dataPlot']['title'] = r'$k_\sigma[p]$'
            self.plot_data()
            
            self.plot_enh()
            
            if min(T1p)<0.1:
                self.guiDict['t1Plot']['ymin'] = 0
            else:
                self.guiDict['t1Plot']['ymin'] =  min(self.guiDict['dnpLabData']['T1p'])*.85
                
            if max(T1p)>5:
                self.guiDict['t1Plot']['ymax'] = 1
            else:
                self.guiDict['t1Plot']['ymax'] =  max(self.guiDict['dnpLabData']['T1p'])*1.15
            
            self.guiDict['t1Plot']['ytick'] = [self.guiDict['t1Plot']['ymin'],self.guiDict['t1Plot']['ymax']]
            self.guiDict['t1Plot']['ytickLabel'] = [str(round(self.guiDict['t1Plot']['ymin'],1)),str(round(self.guiDict['t1Plot']['ymax'],1))]
            
            self.plot_t1()
            
            print('-----Errors in ksigma-----')
            if self.guiDict['workupFunction']['isWorkup']:
                print('workup (hydration): ' + str(round(hydresults['ksigma'],2)) + ' +/- ' +  str(round(hydresults['ksigma_error'],4)))
                print('workup = ' + str(round(self.guiDict['workupData']['kSigma']/spin_C/self.wrkup_smax,2)) + ' +/- ' +  str(round(self.guiDict['workupData']['kSigma_error']/spin_C/self.wrkup_smax,4)))
            else:
                print('dnpLab (hydration): = ' + str(round(hydresults['ksigma'],2)) + ' +/- ' +  str(round(hydresults['ksigma_error'],4)))
                if self.guiDict['workupFunction']['fit']:
                    print('workup (hydration): ' + str(round(whydresults['ksigma'],2)) + ' +/- ' +  str(round(whydresults['ksigma_error'],4)))
                if self.guiDict['workupFunction']['addWorkup']:
                    print('workup = ' + str(round(self.guiDict['workupData']['kSigma']/spin_C/self.wrkup_smax,2)) + ' +/- ' +  str(round(self.guiDict['workupData']['kSigma_error']/spin_C/self.wrkup_smax,4)))
            
    def matoutput(self):

        flnm1 = QFileDialog.getSaveFileName(self)
        if flnm1[0]:
            flnm = flnm1[0]

        if flnm:
            
            odnpData = {'Epowers' : self.guiDict['dnpLabData']['Epowers'], 'Ep' : self.guiDict['dnpLabData']['Ep'], 'T1powers' : self.guiDict['dnpLabData']['T1powers'], 'T1p' : self.guiDict['dnpLabData']['T1p'], 'T1p_error' : self.guiDict['dnpLabData']['T1p_error'], 'T10' : self.guiDict['dnpLabData']['T10'], 'T10_error' : self.guiDict['dnpLabData']['T10_error'], 'T100' : self.guiDict['dnpLabData']['T100']}
            
            odnpResults = {'kSigmas' : self.guiDict['hydrationResults']['ksigma_array'], 'kSigmas_fit' : self.guiDict['hydrationResults']['ksigma_fit']}

            savemat(flnm + '.mat', {'odnp' : odnpData, 'ksig' : odnpResults}, oned_as='column')
            
            e = []
            e.append(1)
            epow = []
            epow.append(0)
            ksig = []
            ksig.append(0)
            ksig_fit = []
            ksig_fit.append(0)
            for k in range(0,len(self.guiDict['dnpLabData']['Epowers'])):
                e.append(self.guiDict['dnpLabData']['Ep'][k])
                epow.append(self.guiDict['dnpLabData']['Epowers'][k])
                ksig.append(self.guiDict['hydrationResults']['ksigma_array'][k])
                ksig_fit.append(self.guiDict['hydrationResults']['ksigma_fit'][k])
            
            dfE = pd.DataFrame(list(zip(*[list(epow), list(e), list(ksig), list(ksig_fit)])))
            dfE.to_csv(flnm + '_E_ksig.csv', index=False, header=['Epowers','Ep','ksigma','ksigma_fit'])
            
            t = []
            t.append(self.guiDict['dnpLabData']['T10'])
            terr = []
            terr.append(self.guiDict['dnpLabData']['T10_error'])
            tpow = []
            tpow.append(0)
            for k in range(0,len(self.guiDict['dnpLabData']['T1powers'])):
                t.append(self.guiDict['dnpLabData']['T1p'][k])
                terr.append(self.guiDict['dnpLabData']['T1p_error'][k])
                tpow.append(self.guiDict['dnpLabData']['T1powers'][k])
                
            dfT1 = pd.DataFrame(list(zip(*[list(tpow), list(t), list(terr)])))
            dfT1.to_csv(flnm + '_T1.csv', index=False, header=['T1powers','T1p','T1p error'])
        

    def changeCenterValue(self, cvalue):
        
        if self.guiDict['guiFunction']['sliders']:
            self.guiDict['processSpectrum']['integrationCenter'] = cvalue
            self.autophaseCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass
        
    def changeWindowValue(self, wvalue):
    
        if self.guiDict['guiFunction']['sliders']:
            self.guiDict['processSpectrum']['integrationWidth'] = wvalue
            self.autophaseCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass
            
    def changePhaseValue(self, pvalue):
        
        if self.guiDict['guiFunction']['sliders']:
            self.guiDict['processSpectrum']['phaseFactor'] = pvalue/1000
            self.autophaseCheckbox.setChecked(False)
            self.adjustSliders()
        else:
            pass
    
    def autoPhaseCheck(self):
        
        if self.guiDict['guiFunction']['sliders']:
            if self.autophaseCheckbox.isChecked():
                self.guiDict['processSpectrum']['phaseFactor'] = 0
                self.processData()
            else:
                pass
        else:
            pass
        
    def linearCheck(self):
        
        if self.linearfitCheckbox.isChecked():
            self.order2fitCheckbox.setChecked(False)
        else:
            self.order2fitCheckbox.setChecked(True)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
        
    def order2Check(self):
    
        if self.order2fitCheckbox.isChecked():
            self.linearfitCheckbox.setChecked(False)
        else:
            self.linearfitCheckbox.setChecked(True)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
            
    def exfT1pCheck(self):
    
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
            
    def tetheredCheck(self):
    
        if self.tetheredCheckbox.isChecked():
            self.bulkCheckbox.setChecked(False)
        else:
            self.bulkCheckbox.setChecked(True)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
    
    def bulkCheck(self):
        
        if self.bulkCheckbox.isChecked():
            self.tetheredCheckbox.setChecked(False)
        else:
            self.tetheredCheckbox.setChecked(True)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
    
    def edithydText(self):
    
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
            
    def show_wrkupCheck(self):

        if self.show_wrkupCheckbox.isChecked():
            self.guiDict['workupFunction']['show'] = True
        else:
            self.guiDict['workupFunction']['show'] = False
            self.guiDict['workupFunction']['fit'] = False
            self.fit_wrkupCheckbox.setChecked(False)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
    
    def fit_wrkupCheck(self):
    
        if self.fit_wrkupCheckbox.isChecked() and self.show_wrkupCheckbox.isChecked():
            self.guiDict['workupFunction']['fit'] = True
        else:
            self.guiDict['workupFunction']['fit'] = False
            self.fit_wrkupCheckbox.setChecked(False)
            
        if self.guiDict['guiFunction']['hydrationEdits']:
            self.run_hydration()
        else:
            pass
    
    """
    --Plot Colors Key--
    dark_green = '#46812B'
    light_green = '#67AE3E'
    dark_grey = '#4D4D4F'
    light_grey = '#A7A9AC'
    orange = '#F37021'
    ucsb blue = '#004D9F'
    color='#004D9F',marker='o',linestyle='none'
    """
    
    def plot_data(self):

        self.dataplt.axes.cla()

        if self.guiDict['dataPlot']['plotksig']:
            
            self.dataplt.axes.plot(self.guiDict['dnpLabData']['Epowers'],self.guiDict['hydrationResults']['ksigma_array'],color='#46812B',marker='o',linestyle='none',label= self.ksiglabel + r' $k_\sigma$[p]')
            self.dataplt.axes.plot(self.guiDict['dnpLabData']['Epowers'],self.plothydresults['ksigma_fit'],'#F37021',label=r'hydration fit')
            
            indx_h = max(self.plothydresults['ksigma_array'])*.8
            self.dataplt.axes.set_ylim(0, max(self.plothydresults['ksigma_array'])*1.1)
            
            if self.guiDict['workupFunction']['show']:
                self.dataplt.axes.plot(self.guiDict['workupData']['Epowers'],self.guiDict['whydrationResults']['ksigma_array'],color='#004D9F',marker='o',linestyle='none',label= r'workup $k_\sigma$[p]')
                self.dataplt.axes.plot(self.guiDict['dnpLabData']['Epowers'],self.guiDict['whydrationResults']['ksigma_fit'],color='#004D9F',linestyle=':',label = 'workup fit')

                indexes = [0, .11, .31, .41, .51, .61]
                self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.645, indx_h-(.21*indx_h), r'workup $k_\sigma = $' + str(round(self.guiDict['workupData']['kSigma']/float(self.slcEdit.text())/self.wrkup_smax,2)), fontsize=12)
            else:
                indexes = [.11, .21, .31, .41, .51, .61]
                
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[0]*indx_h), r'$k_\rho = $' + str(round(self.plothydresults['krho'],2)), fontsize=12)
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[1]*indx_h), r'$k_\sigma = $' + str(round(self.plothydresults['ksigma'],2)), fontsize=12)
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[2]*indx_h), r'$k_{low} = $' + str(round(self.plothydresults['klow'],2)), fontsize=12)
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[3]*indx_h), r'$\xi = $' + str(round(self.plothydresults['coupling_factor'],4)), fontsize=12)
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[4]*indx_h), r'$t_{corr} = $' + str(round(self.plothydresults['tcorr'],2)), fontsize=12)
            d_local = round(self.plothydresults['Dlocal']*1e10,2)
            self.dataplt.axes.text(max(self.guiDict['dnpLabData']['Epowers'])*.75,indx_h-(indexes[5]*indx_h), r'$D_{local} = $' + str(d_local) + r'$e^{-10}$', fontsize=12)
            self.dataplt.axes.set_yticks([0 , max(self.plothydresults['ksigma_array'])])
            self.dataplt.axes.set_yticklabels(['0' , str(round(max(self.plothydresults['ksigma_array']),1))])
            self.dataplt.axes.legend()
        else:
            self.dataplt.axes.plot(self.guiDict['dataPlot']['xdata'],self.guiDict['dataPlot']['ydata'])
            self.dataplt.axes.set_xlim(self.guiDict['dataPlot']['xmin'],self.guiDict['dataPlot']['xmax'])
            self.dataplt.axes.set_yticks([0])
            self.dataplt.axes.set_yticklabels('0')
            
        self.dataplt.axes.set_title(self.guiDict['dataPlot']['title'])
        self.dataplt.axes.set_xticks([])
        self.dataplt.draw()
        
    def plot_enh(self):

        self.enhplt.axes.cla()
        if self.guiDict['enhancementPlot']['plotT1fit']:
            self.enhplt.axes.plot(self.guiDict['t1Plot']['tau'],self.guiDict['t1Plot']['t1Amps'],color='#46812B',marker='o',linestyle='none')
            new_xax = np.r_[np.min(self.guiDict['t1Plot']['tau']):np.max(self.guiDict['t1Plot']['tau']):100j]
            self.enhplt.axes.plot(new_xax,self.guiDict['t1Plot']['t1Fit'],'#F37021')
            self.enhplt.axes.text(max(self.guiDict['t1Plot']['tau'])*.65, max(self.guiDict['t1Plot']['t1Amps'])*.3, 'T1 = ' + str(round(self.guiDict['t1Plot']['t1Val'], 4)), fontsize=10)
        else:
            if self.guiDict['enhancementPlot']['plotEpfit']:
                self.enhplt.axes.plot(self.guiDict['enhancementPlot']['xdata'],self.guiDict['enhancementPlot']['ydata'],color='#46812B',marker='o',linestyle='none',label=self.ksiglabel)
                self.enhplt.axes.plot(self.guiDict['enhancementPlot']['xdata'],self.plothydresults['uncorrected_Ep'],'#F37021',label='Hydration Fit')
                if self.guiDict['workupFunction']['show']:
                    self.enhplt.axes.plot(self.guiDict['workupData']['Epowers'],self.guiDict['workupData']['Ep'],color='#004D9F',marker='o',linestyle='none',label='workup')
                self.enhplt.axes.legend()
            else:
                self.enhplt.axes.plot(self.guiDict['enhancementPlot']['xdata'],self.guiDict['enhancementPlot']['ydata'],color='#46812B',marker='o',linestyle='none')
                
        self.enhplt.axes.set_title(self.guiDict['enhancementPlot']['title'])
        self.enhplt.axes.set_xticks([])
        self.enhplt.axes.set_yticks(self.guiDict['enhancementPlot']['ytick'])
        self.enhplt.axes.set_yticklabels(self.guiDict['enhancementPlot']['ytickLabel'])
        self.enhplt.draw()
        
    def plot_t1(self):

        self.t1plt.axes.cla()
        if self.guiDict['t1Plot']['plotT1interp']:
            self.t1plt.axes.plot(self.guiDict['t1Plot']['xdata'],self.guiDict['t1Plot']['ydata'],color='#46812B',marker='o',linestyle='none',label=self.ksiglabel)
            self.t1plt.axes.plot(self.guiDict['dnpLabData']['Epowers'],self.plothydresults['interpolated_T1'],'#F37021',label='Interpolation')
            if self.guiDict['workupFunction']['show']:
                self.t1plt.axes.plot(self.guiDict['workupData']['T1powers'],self.guiDict['workupData']['T1p'],color='#004D9F',marker='o',linestyle='none',label='workup')
            self.t1plt.axes.legend()
        else:
            self.t1plt.axes.plot(self.guiDict['t1Plot']['xdata'],self.guiDict['t1Plot']['ydata'],color='#46812B',marker='o',linestyle='none')
        self.t1plt.axes.set_ylim(self.guiDict['t1Plot']['ymin'], self.guiDict['t1Plot']['ymax'])
        self.t1plt.axes.set_title(self.guiDict['t1Plot']['title'])
        self.t1plt.axes.set_xticks([])
        self.t1plt.axes.set_yticks(self.guiDict['t1Plot']['ytick'])
        self.t1plt.axes.set_yticklabels(self.guiDict['t1Plot']['ytickLabel'])
        self.t1plt.draw()
            
class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=2, height=1, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

def main_func():
    app = QApplication(sys.argv)
    ex = hydrationGUI()
    sys.exit(app.exec_())
if __name__ == '__main__':
    main_func()
