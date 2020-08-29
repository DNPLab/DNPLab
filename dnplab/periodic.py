import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QWidget, QLabel, QFrame
import numpy as np

light_green = QtGui.QColor(103, 174, 62)
dark_green = QtGui.QColor(70, 129, 43)

light_grey = QtGui.QColor(167, 169, 172)
dark_grey = QtGui.QColor(88, 89, 91)

orange = QtGui.QColor(243, 112, 33)

element_width = 40
element_height = 50
#symbol_height = 80
atomic_number_height = 15
border_width = 5
border_height = 5
font = "Arial"
symbol_point_size = 14
atomic_number_point_size = 8
field_point_size = 8

# 4k display
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

with open('isotopes_data.csv', 'r') as f:

    elements = {}

    header_line = f.readline()
    header_list = header_line.strip().split(',')

    while True:
        isotope = {}
        line = f.readline()
        if line == '':
            break
        data = line.strip().split(',')
        for column, header in enumerate(header_list):
            isotope[header] = data[column]

        if isotope['Symbol'] not in elements:
            elements[isotope['Symbol']] = []

        elements[isotope['Symbol']].append(isotope)
 
with open('atomic_number_symbol_row_column.csv', 'r') as f:

    periodic_table_data = {}

    header_line = f.readline()
    header_list = header_line.strip().split(',')

    line = ''
    while True:

        line = f.readline()
        if line == '':
            break

        data = line.strip().split(',')

        info = {}
        for column, header in enumerate(header_list):
            info[header] = data[column]

        periodic_table_data[info['Symbol']] = info


class ElementWidget(QWidget):
    def __init__(self, parent = None, symbol = 'H', atomic_number = 1, row_ix = 0, column_ix = 0, color = light_green, isotopes = {}):
        QWidget.__init__(self, parent = parent)

        self.parent = parent
        self.resize(element_width, element_height)
        self.setGeometry(QtCore.QRect(int(column_ix*element_width + border_width * (column_ix + 1)),
                                     int(row_ix*element_height + border_height * (row_ix + 1)),
                                     int(element_width),
                                     int(element_height)))

        self.column_ix = column_ix
        self.row_ix = row_ix
        self.color = color

        self.set_symbol(symbol)
        self.set_atomic_number(atomic_number)
        p = self.palette()
        p.setColor(self.backgroundRole(), self.color)
        self.setPalette(p)
        self.setAutoFillBackground(True)

        self.isotopes = isotopes

    def set_symbol(self, symbol):
        self.symbol = QLabel(self.parent)
        self.symbol.setGeometry(QtCore.QRect(int(border_width * (self.column_ix+1) + element_width*self.column_ix), 
                                             int(border_height * (self.row_ix+1) + element_height*self.row_ix + atomic_number_height), 
                                             int(element_width),
                                             int(element_height - atomic_number_height)))

        symbol_font = QtGui.QFont()
        symbol_font.setFamily(font)
        symbol_font.setPointSize(symbol_point_size)
        self.symbol.setFont(symbol_font)
        self.symbol.setAlignment(QtCore.Qt.AlignCenter)
        self.symbol.setObjectName('symbol')
        self.symbol.setText(symbol)

        self.symbol.mousePressEvent = self.mousePressEvent


    def set_atomic_number(self, atomic_number):
        self.atomic_number = QLabel(self.parent)
        self.atomic_number.setGeometry(QtCore.QRect(int(border_width * (self.column_ix+1) + element_width*self.column_ix), 
                                                    int(border_height * (self.row_ix+1) + element_height*self.row_ix),
                                                    int(element_width),
                                                    int(atomic_number_height)))
        atomic_number_font = QtGui.QFont()
        atomic_number_font.setFamily(font)
        atomic_number_font.setPointSize(atomic_number_point_size)
        self.atomic_number.setFont(atomic_number_font)
        self.atomic_number.setAlignment(QtCore.Qt.AlignCenter)
        self.atomic_number.setObjectName('atomic_number')
        self.atomic_number.setText(str(atomic_number))

        self.atomic_number.mousePressEvent = self.mousePressEvent

    def mousePressEvent(self, QMouseEvent):
        self.parent.set_selected_element(self)
        self.parent.statusBar().showMessage('Clicked on %s'%str(self.symbol.text()))

    def set_color(self, color):
        self.color = color
        p = self.palette()
        p.setColor(self.backgroundRole(), self.color)
        self.setPalette(p)


class Ui_MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setWindowTitle('Periodic Table')
        self.resize(element_width*18 + border_width*20,
                    element_height * 13 + border_height*14)

        self.frequency_header = 'Frequency (MHz)'
        self.table_column_header = ['Mass Number', 'Symbol', 'Isotope', 'Spin', self.frequency_header, 'Magnetogyric ratio (10^7 rad / (T s))', 'Quadrupole Moment (Q/fm^2)','Frequency Ratio', 'Natural Abundance', 'Relative Receptivity 1H', 'Relative Receptivity 13C']
        self.table_column_header_alias = ['', 'Symbol','Isotope', 'Spin', 'Frequency (MHz)', 'Magnetogyric ratio (10^7 rad / (T s))', 'Quadrupole Moment (Q/fm^2)', 'Frequency Ratio', 'Natural Abundance', 'Relative Receptivity 1H', 'Relative Receptivity 13C']

        self.frequency_ix = self.table_column_header.index(self.frequency_header)
        self.gyro_ix = self.table_column_header.index('Magnetogyric ratio (10^7 rad / (T s))')

        self.statusBar().showMessage('Ready')

        self.list_data = QtWidgets.QTableWidget(self)
        self.list_data.setGeometry(QtCore.QRect(int(element_width * 0 + border_width * 1),
                                                int(element_height * 10.5),
                                                int(element_width * 18 + border_width * 17),
                                                int(element_height * 3 + border_height*4)))

        self.list_data.verticalHeader().setVisible(False)
        self.list_data.horizontalHeader().setVisible(False)
        self.list_data.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)


        self.field_label = QtWidgets.QLabel(self)
        self.field_label.setGeometry(QtCore.QRect(int(element_width * 0 + border_width * 1),
                                                  int(element_height * 7.5 + border_height*8),
                                                  int(element_width * 2 + border_width * 3),
                                                  int(element_height * 0.5 + border_height*1)))
        self.field_label.setText('Field (T)')
        field_font = QtGui.QFont()
        field_font.setFamily(font)
        field_font.setPointSize(field_point_size)
        self.field_label.setFont(field_font)
        self.field_label.setAlignment(QtCore.Qt.AlignCenter)

        self.field_spinbox = QtWidgets.QDoubleSpinBox(self)
        self.field_spinbox.setGeometry(QtCore.QRect(int(element_width * 0 + border_width * 1),
                                                    int(element_height * 8 + border_height*8.5),
                                                    int(element_width * 2 + border_width * 3),
                                                    int(element_height * 0.5 + border_height*1)))

        self.field_spinbox.setValue(0.35)
        self.field_spinbox.setSingleStep(0.0001)
        self.field_spinbox.setDecimals(6)
        self.field_spinbox.setRange(0,1000)

        self.field_spinbox.valueChanged.connect(self.update_table)

        selected_element_symbol = 'H'

        for element_symbol in periodic_table_data:

            data = periodic_table_data[element_symbol]

            atomic_number = data['Atomic Number']
            symbol = element_symbol

            row = data['Row']
            column = data['Column']

            if element_symbol in elements:
                isotopes = elements[element_symbol]
            else:
                isotopes = {}
            if symbol == selected_element_symbol:
                self.selected_element = ElementWidget(self, symbol = symbol, atomic_number = atomic_number, row_ix = float(row), column_ix = float(column), color = light_green, isotopes = isotopes)
            else:
                self.element = ElementWidget(self, symbol = symbol, atomic_number = atomic_number, row_ix = float(row), column_ix = float(column), color = light_green, isotopes = isotopes)

        self.set_selected_element(self.selected_element)

        self.setup_table()

    def set_selected_element(self, element):
        self.selected_element.set_color(light_green)
        self.selected_element = element
        self.selected_element.set_color(orange)
        self.setup_table()

    def update_table(self):
        for row in range(self.list_data.rowCount()):
            try:
                gyro = float(self.list_data.item(row, self.gyro_ix).text())
                frequency = gyro * float(self.field_spinbox.value()) / (2.*np.pi) * 10.
                self.list_data.setItem(row, self.frequency_ix, QtWidgets.QTableWidgetItem('%0.05f'%frequency))
            except:
                pass

        header = self.list_data.horizontalHeader()
        for column in range(self.list_data.columnCount()):
            header.setSectionResizeMode(column, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)

    def setup_table(self):
        isotopes = self.selected_element.isotopes
        self.list_data.clear()
        while (self.list_data.rowCount() > 0):
            self.list_data.removeRow(0)

        while (self.list_data.columnCount() > 0):
            self.list_data.removeColumn(0)

        for ix in range(len(self.table_column_header)):
            self.list_data.insertColumn(ix)
        
        font = QtGui.QFont()
        font.setBold(True)

        self.list_data.insertRow(self.list_data.rowCount())
        for header_ix, header in enumerate(self.table_column_header_alias):
            item = QtWidgets.QTableWidgetItem(str(header))
            item.setFont(font)
            self.list_item = self.list_data.setItem(self.list_data.rowCount() - 1, header_ix, item)

        for row_ix, isotope in enumerate(isotopes):
            self.list_data.insertRow(self.list_data.rowCount())

            for header_ix, header in enumerate(self.table_column_header):
                if header != self.frequency_header:
                    item = QtWidgets.QTableWidgetItem(str(isotopes[row_ix][header]))
                    if header_ix == 0:
                        item.setTextAlignment(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

                    self.list_data.setItem(self.list_data.rowCount() - 1, header_ix, item)

        self.update_table()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    main = Ui_MainWindow()
    main.show()
    sys.exit(app.exec_())
