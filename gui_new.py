#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 15:11:45 2018

@author: alexk
"""

import sys
import platform
import numpy as np # this has to be imported before the ones in line 11 and 12
from math import floor
import control
import control.matlab as matlab

from PyQt5.QtWidgets import (QMainWindow, QApplication, QDialog, QLineEdit, 
                             QVBoxLayout, QAction, QMessageBox, QFileDialog,
                             QSizePolicy, QComboBox, QPushButton, QHBoxLayout,
                             QGridLayout, QLabel)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR
from PyQt5.QtGui import QIcon

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure 

class MainWindow(QMainWindow) :
    
    def __init__(self, parent=None) :
        super(MainWindow, self).__init__(parent)
        
        # Add an Icon
        # self.setWindowIcon(QIcon('troll.png'))

        #######################################################################
        # ADD MENU ITEMS
        #######################################################################
        
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs = QAction("&Save As", self)
        self.actionSaveAs.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs, self.actionQuit])
        
        # Create the Help menu
        self.menuHelp = self.menuBar().addMenu("&Help")
        self.actionAbout = QAction("&About",self)
        self.actionAbout.triggered.connect(self.about)
        self.menuHelp.addActions([self.actionAbout])
        
        #######################################################################
        # CREATE CENTRAL WIDGET
        #######################################################################

        self.widget = QDialog()
        self.plot = MatplotlibCanvas()
        
        # Input and Output Boxes
        self.paramEdit = QLineEdit("np.linspace(0,10,11)")
        self.output = QLineEdit("Output Values")
        
        # Horizontal Time entries
        TimeStep = QLabel('Time Step =')
        SimLen = QLabel('Final Time =')
        unit0 = QLabel('sec. ')
        unit1= QLabel('sec. ')
        self.dtEdit = QLineEdit('0.05')
        self.dtEdit.setFixedWidth(50)
#        TimeStep.setBuddy(self.dtEdit)
        self.simEdit = QLineEdit('20')
        self.simEdit.setFixedWidth(50)
        
        # System Parameter
        mass = QLabel('Mass =')
        Ks = QLabel('Ks =')
        Kd = QLabel('Kd =')
        unitMass = QLabel('kg ')
        unitKs = QLabel('N/m ')
        unitKd = QLabel('N-sec/m ')
        self.massEdit = QLineEdit('1')
        self.KsEdit = QLineEdit('4')
        self.KdEdit = QLineEdit('1')      
        self.massEdit.setFixedWidth(50)
        self.KsEdit.setFixedWidth(50)
        self.KdEdit.setFixedWidth(50)

        # State Space
        FGrid = QGridLayout()
        F = QLabel('F =')
        self.FEdit1, self.FEdit2 = QLineEdit('0.0'), QLineEdit('1.0 ')
        self.FEdit3, self.FEdit4 = QLineEdit('-Ks/m'), QLineEdit('-Kd/m ')
        self.FEdit1.setFixedWidth(40)
        self.FEdit2.setFixedWidth(40)
        self.FEdit3.setFixedWidth(40)
        self.FEdit4.setFixedWidth(40)
        FGrid.addWidget(F, 3, 0)
        FGrid.addWidget(self.FEdit1, 3, 1, 1, 1)
        FGrid.addWidget(self.FEdit2, 3, 2, 1, 1)
        FGrid.addWidget(self.FEdit3, 4, 1, 1, 1)
        FGrid.addWidget(self.FEdit4, 4, 2, 1, 1)
        
        GGrid = QGridLayout()
        G = QLabel('G =')
        self.GEdit1, self.GEdit2 = QLineEdit('0.0'), QLineEdit('1/m')
        self.GEdit1.setFixedWidth(40)
        self.GEdit2.setFixedWidth(40)
        GGrid.addWidget(G, 3, 3)
        GGrid.addWidget(self.GEdit1, 3, 4)
        GGrid.addWidget(self.GEdit2, 4, 4)
        
        
        grid = QGridLayout()
#        TimeStep = QLabel('Time Step =')
#        SimLen = QLabel('Final Time =')
#        self.dtEdit = QLineEdit('0.05')
#        self.dtEdit.setFixedWidth(50)
#        TimeStep.setBuddy(self.dtEdit)
#        self.simEdit = QLineEdit('20')
#        self.simEdit.setFixedWidth(50)
#        mass = QLabel('Mass =')
#        Ks = QLabel('Ks =')
#        Kd = QLabel('Kd =')
#        self.massEdit = QLineEdit('1')
#        self.KsEdit = QLineEdit('4')
##        self.KdEdit = QLineEdit('1')
#        F = QLabel('F =')
#        G = QLabel('G =')
#        self.FEdit1, self.FEdit2 = QLineEdit('0.0'), QLineEdit('1.0')
#        self.FEdit3, self.FEdit4 = QLineEdit('-Ks/m'), QLineEdit('-Kd/m')
#        self.GEdit1, self.GEdit2 = QLineEdit('0.0'), QLineEdit('1/m')
        distNoise = QLabel('Disturbance Cov. =')
        self.distNoiseEdit = QLineEdit('0.0005')
        AccelNoise = QLabel('Accel. Meas. Noise =')
        self.AccelNoiseEdit = QLineEdit('0.02')
        measNoise = QLabel('Meas. Noise =')
        self.measNoiseEdit = QLineEdit('0.001')
        inputSignal = QLabel('Input =')
        self.inputSignalEdit = QLineEdit('10*np.sin(2*np.pi*0.05*t)')

#        grid.addWidget(TimeStep, 1, 0)
#        grid.addWidget(self.dtEdit, 1, 1)
#        grid.addWidget(SimLen, 1, 4)
#        grid.addWidget(self.simEdit, 1, 5)
#        grid.addWidget(mass, 2, 0)
#        grid.addWidget(self.massEdit, 2, 1)
#        grid.addWidget(Ks, 2, 2)
#        grid.addWidget(self.KsEdit, 2, 3)
#        grid.addWidget(Kd, 2, 4)
#        grid.addWidget(self.KdEdit, 2, 5)
#        grid.addWidget(F, 3, 0)
#        grid.addWidget(self.FEdit1, 3, 1, 1, 1)
#        grid.addWidget(self.FEdit2, 3, 2, 1, 1)
#        grid.addWidget(G, 3, 3)
#        grid.addWidget(self.GEdit1, 3, 4)
#        grid.addWidget(self.GEdit2, 3, 5)
#        grid.addWidget(self.FEdit3, 4, 1, 1, 1)
#        grid.addWidget(self.FEdit4, 4, 2, 1, 1)
        grid.addWidget(distNoise, 5, 0)
        grid.addWidget(self.distNoiseEdit, 5, 1)
        grid.addWidget(measNoise, 5, 2)
        grid.addWidget(self.measNoiseEdit, 5, 3)
        grid.addWidget(AccelNoise, 5, 4)
        grid.addWidget(self.AccelNoiseEdit, 5, 5)
        grid.addWidget(inputSignal, 6, 0)
        grid.addWidget(self.inputSignalEdit, 6, 1)
        
        # Horizontal Buttons Layout
        self.b1 = QPushButton('Run')
        self.b2 = QPushButton('Clear')
        
        Buttons = QHBoxLayout()
        Buttons.addWidget(self.b1)
        Buttons.addWidget(self.b2)
        Buttons.addStretch(1)
        
        TimeLayout = QHBoxLayout()
        TimeLayout.addWidget(TimeStep)
        TimeLayout.addWidget(self.dtEdit)
        TimeLayout.addWidget(unit0)
        TimeLayout.addWidget(SimLen)
        TimeLayout.addWidget(self.simEdit)
        TimeLayout.addWidget(unit1)
        TimeLayout.addStretch(1)
        
        SysParam = QHBoxLayout()
        SysParam.addWidget(mass)
        SysParam.addWidget(self.massEdit)
        SysParam.addWidget(unitMass)
        SysParam.addWidget(Ks)
        SysParam.addWidget(self.KsEdit)
        SysParam.addWidget(unitKs)
        SysParam.addWidget(Kd)
        SysParam.addWidget(self.KdEdit)
        SysParam.addWidget(unitKd)
        SysParam.addStretch(1)
        
        States = QHBoxLayout()
        States.addLayout(FGrid)
        States.addLayout(GGrid)        
        States.addStretch(1)
        
        # GUI Layout
        layout = QVBoxLayout()
        layout.addLayout(TimeLayout)
        layout.addLayout(SysParam)
        layout.addLayout(States)
        layout.addLayout(Buttons)
        layout.addStretch(1)
        layout2 = QHBoxLayout()
        layout2.addWidget(self.plot)
        layout2.addLayout(layout)
#        layout.addWidget(self.output)
        
        self.widget.setLayout(layout2)
        self.setCentralWidget(self.widget)      

        # Signals:
        # Pressing Enter returns Output
        self.paramEdit.returnPressed.connect(self.runButton)
        # Button clicks to run and clear
        self.b1.clicked.connect(self.runButton)
        self.b2.clicked.connect(self.clearPlot)      

    def kalmanfilter(self):
        """
        """
        m = self.massEdit.text()
        Ks = self.KsEdit.text()
        Kd = self.KdEdit.text()
        F1 = self.FEdit1.text()
        F2 = self.FEdit2.text()
        F3 = self.FEdit3.text()
        F4 = self.FEdit4.text()
        G1 = self.GEdit1.text()
        G2 = self.GEdit2.text()
        dt = self.dtEdit.text()
        time = self.simEdit.text()
        m = eval(m)
        Ks = eval(Ks)
        Kd = eval(Kd)
        F1 = eval(F1)
        F2 = eval(F2)
        F3 = eval(F3)
        F4 = eval(F4)
        G1 = eval(G1)
        G2 = eval(G2)
        dt = eval(dt)
        time = eval(time)
        
        F = np.array([[F1, F2], [F3, F4]])
        G = np.array([[G1], [G2]])
        H = np.array([1.0, 0.0])
        J = 0
        sysc = control.ss(F,G,H,J)
                
        # initialize simulations and Kalman Filter
        tf = max(time, floor(1000*dt))
        n = 0
        t = np.zeros((int(tf/dt),1))
        for i in np.arange(0,tf,dt):
            t[n] = i
            n+=1
        nSteps = len(t)
        Qc2 = self.AccelNoiseEdit.text()
        Qc2 = eval(Qc2)
        Qd2 = Qc2/dt
        w = np.sqrt(Qd2)*np.random.randn(nSteps,1)
        R = self.measNoiseEdit.text()
        R = eval(R)
        v = np.sqrt(R)*np.random.randn(nSteps,1)
        X0 = [0.25, 0]
        sysd = control.c2d(sysc,dt)
        [Phi,Gamma,H,J] = control.ssdata(sysd)
        K = np.zeros((2,nSteps))
        xp = np.zeros((2,nSteps))
        Pp = np.identity(2)
        xp[:,0]=[0, 0]
        
        # initialize INS/GPS
        Fi = np.array([[0.0, 1.0],[0.0, 0.0]])
        Gi = np.array([[0],[1]])
        Hi = H*1
        Ji = 0
        Qdi2 = self.AccelNoiseEdit.text()
        Qdi2 = eval(Qdi2)
        Qdi2 = Qdi2*Qdi2
        wi = np.sqrt(Qdi2)*np.random.randn(nSteps,1)
        Ri = R*1
        sysdi = control.c2d(control.ss(Fi,Gi,Hi,Ji),dt)
        Phii, Gammai, Hi, Ji = control.ssdata(sysdi)
        Ki = K*1
        xpi = xp*1
        Ppi = Pp*1
        
        u = self.inputSignalEdit.text()
        u = eval(u)
        u = u.T
        yt, t, xt = matlab.lsim(sysc, w+u.T, t, X0)
        xt = xt.T
        y = yt+v
        
        accel = np.array([0, 1])@(F@xt + G*(u+np.transpose(w)))
        accelmeas = accel + np.transpose(wi)
        for k in range(1, nSteps):
            xm = Phi@xp[:,[(k-1)]] + Gamma*u[:,k-1]
            Pm = Phi@Pp@np.transpose(Phi) + Gamma*Qd2*np.transpose(Gamma)
            
            K[:,[k]] = Pm@np.transpose(H)*np.linalg.inv(H@Pm@np.transpose(H) + R)
            xp[:,[k]] = xm + K[:,[k]]*(y[0,k]-H*xm)
            Pp = Pm - K[:,[k]]@H@Pm
            
            # INS/GPS Filtering
            xmi = Phii@xpi[:,[(k-1)]] + Gammai*accelmeas[:, k-1]
            Pmi = Phii@Ppi@np.transpose(Phii) + Gammai*Qdi2*np.transpose(Gammai)
            
            Ki[:,[k]] = Pmi@np.transpose(Hi)*np.linalg.inv(Hi@Pmi@np.transpose(Hi) + Ri)
            xpi[:,[k]] = xmi + Ki[:,[k]]*(y[0,k]-Hi*xmi)
            Ppi = Pmi - Ki[:,[k]]@Hi@Pmi
        
        self.plot.redraw(t, yt, t, xpi[0,:])

    def saveas(self):
        """Save input and output to a text file as seperate columns
        """
        
        name = QFileDialog.getSaveFileName(self, "saveas")[0]
        f = open(name, 'w')
        x = self.paramEdit.text()
        x = eval(x)
        if len(x) > 1 :
            x = np.array(x)
        output = eval(str(self.funcEdit.currentText()))
        f.write("  x    f(x)\n")
        for i in range(len(x)):
            f.write("%5.2f %5.2f\n" % (x[i], output[i]))
        f.close()
                
    def about(self):
        QMessageBox.about(self, 
            "About Function Evaluator",
            """<b>Function Evaluator</b>
               <p>Copyright &copy; 2016 Jeremy Roberts, All Rights Reserved.
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))
        
    def runButton(self):
        try:
            self.kalmanfilter()

        except:
            self.output.setText(
                    "Input error! Check your input parameters and function.")
        
    def clearPlot(self):
        self.plot.clear()
        self.output.setText("Output Values")            
        

class MatplotlibCanvas(FigureCanvas):
    """ This is borrowed heavily from the matplotlib documentation;
        specifically, see:
        http://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html
    """
    def __init__(self):
        
        # Initialize the figure and axes
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        
        # Give it some default empty plot
        x = 0
        y = 0
        self.axes.plot(x, y)
        self.axes.set_xlabel('Time (s)')
        self.axes.set_ylabel('y(x)')
        self.axes.set_title('Kalman Filter Simulation')
        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
    def redraw(self, x, y, u, v):
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.set_title('Kalman Filter Simulation')
        self.axes.plot(x, y, u, v)
        self.draw()
        
    def clear(self):
        self.axes.clear()
        self.axes.set_title('Kalman Filter Simulation') 
        self.draw()
        
app = QApplication(sys.argv)
form = MainWindow()
form.resize(950, 700)
form.show()
app.exec_()
