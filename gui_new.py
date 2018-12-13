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
#import matplotlib.animation as animation

from PyQt5.QtWidgets import (QMainWindow, QApplication, QDialog, QLineEdit, 
                             QVBoxLayout, QAction, QMessageBox, QFileDialog,
                             QSizePolicy, QPushButton, QHBoxLayout, QLabel,
                             QGridLayout, QShortcut)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR, Qt
from PyQt5.QtGui import QIcon, QKeySequence

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure 

class MainWindow(QMainWindow) :
    
    def __init__(self, parent=None) :
        super(MainWindow, self).__init__(parent)
        
        # Add an Icon
        self.setWindowIcon(QIcon('guiIcon.png'))
        
        # Add Shortcuts
        self.shortPlot1 = QShortcut(QKeySequence('Ctrl+p'), self)
        self.shortPlot1.activated.connect(self.runButton1)
        self.shortPlot2 = QShortcut(QKeySequence('Ctrl+v'), self)
        self.shortPlot2.activated.connect(self.runButton2)
        self.shortPlot3 = QShortcut(QKeySequence('Ctrl+c'), self)
        self.shortPlot3.activated.connect(self.clearPlot)
        self.shortPlot4 = QShortcut(QKeySequence('Ctrl+z'), self)
        self.shortPlot4.activated.connect(self.zoom)
        self.shortPlot5 = QShortcut(QKeySequence('Ctrl+a'), self)
        self.shortPlot5.activated.connect(self.pan)
        self.shortPlot6 = QShortcut(QKeySequence('Ctrl+h'), self)
        self.shortPlot6.activated.connect(self.home)
        self.shortSave = QShortcut(QKeySequence('Ctrl+s'), self)
        self.shortSave.activated.connect(self.saveas)
        self.shortClose = QShortcut(QKeySequence('Ctrl+w'), self)
        self.shortClose.activated.connect(self.close)

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
        self.toolbar = NavigationToolbar2QT(self.plot, self)
        self.toolbar.hide()
        
        # Input and Output Boxes
        self.paramEdit = QLineEdit("np.linspace(0,10,11)")
        self.output = QLineEdit("Output Values")
        
        # Horizontal Time entries
        TimeStep = QLabel('Time Step =')
        SimLen = QLabel('   Final Time =')
        unit0 = QLabel('sec. ')
        unit1= QLabel('sec. ')
        self.dtEdit = QLineEdit('0.05')
        self.dtEdit.setFixedWidth(50)
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
        G = QLabel('    G =')
        self.GEdit1, self.GEdit2 = QLineEdit('0.0'), QLineEdit('1/m')
        self.GEdit1.setFixedWidth(40)
        self.GEdit2.setFixedWidth(40)
        GGrid.addWidget(G, 3, 3)
        GGrid.addWidget(self.GEdit1, 3, 4)
        GGrid.addWidget(self.GEdit2, 4, 4)
        
        # Initial Condition
        InitCond = QGridLayout()
        IC = QLabel('      I.C. X0 =')
        self.ICEdit1, self.ICEdit2 = QLineEdit('0.25'), QLineEdit('0.0')
        self.ICEdit1.setFixedWidth(40)
        self.ICEdit2.setFixedWidth(40)
        InitCond.addWidget(IC, 3, 5)
        InitCond.addWidget(self.ICEdit1, 3, 6)
        InitCond.addWidget(self.ICEdit2, 3, 7)
        InitCond.setAlignment(Qt.AlignTop)
        
        # Noise Parameter
        distNoise = QLabel('  Disturbance Cov. =')
        self.distNoiseEdit = QLineEdit('0.0005')
        self.distNoiseEdit.setFixedWidth(80)
        distNoise.setBuddy(self.distNoiseEdit)
        AccelNoise = QLabel('Accel. Meas. Noise =')
        self.AccelNoiseEdit = QLineEdit('0.02')
        self.AccelNoiseEdit.setFixedWidth(80)
        measNoise = QLabel('          Meas. Noise =')
        self.measNoiseEdit = QLineEdit('0.001')
        self.measNoiseEdit.setFixedWidth(80)
        
        # Input Signal
        inputSignal = QLabel('Input Signal =')
        self.inputSignalEdit = QLineEdit('10*np.sin(2*np.pi*0.05*t)')

        # Horizontal Buttons Layout
        self.b1 = QPushButton('&Position Tracking')
        self.b2 = QPushButton('&Velocity Tracking')
        self.b3 = QPushButton('&Clear')
        
        # Horizontal Plot Toolbar
        self.b4 = QPushButton('&Zoom')
        self.b5 = QPushButton('P&an')
        self.b6 = QPushButton('&Home')
        
        # Plot window
        Graph = QVBoxLayout()
        Graph.addWidget(self.plot)
        
        # Plot Toolbar
        PlotButtons = QHBoxLayout()
        PlotButtons.addWidget(self.b6)
        PlotButtons.addWidget(self.b4)
        PlotButtons.addWidget(self.b5)
        
        # Run and clear buttons        
        Buttons = QHBoxLayout()
        Buttons.addWidget(self.b1)
        Buttons.addWidget(self.b2)
        Buttons.addWidget(self.b3)
        
        # Buttons and Toolbar
        Button = QVBoxLayout()
        Button.addLayout(Buttons)
        Button.addLayout(PlotButtons)
        
        # Time parameter        
        TimeLayout = QHBoxLayout()
        TimeLayout.addWidget(TimeStep)
        TimeLayout.addWidget(self.dtEdit)
        TimeLayout.addWidget(unit0)
        TimeLayout.addWidget(SimLen)
        TimeLayout.addWidget(self.simEdit)
        TimeLayout.addWidget(unit1)

        # System parameter        
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
        
        # Grid layouts   
        States = QHBoxLayout()
        States.addLayout(FGrid)
        States.addLayout(GGrid)
        States.addLayout(InitCond)
        States.addStretch(1)
        
        # Noise parameter        
        NoiseParam1 = QHBoxLayout()
        NoiseParam1.addWidget(distNoise)
        NoiseParam1.addWidget(self.distNoiseEdit)
        NoiseParam1.addStretch()
        NoiseParam2 = QHBoxLayout()
        NoiseParam2.addWidget(measNoise)
        NoiseParam2.addWidget(self.measNoiseEdit)
        NoiseParam2.addStretch()
        NoiseParam3 = QHBoxLayout()
        NoiseParam3.addWidget(AccelNoise)
        NoiseParam3.addWidget(self.AccelNoiseEdit)
        NoiseParam3.addStretch(1)

        InputParam = QHBoxLayout()
        InputParam.addWidget(inputSignal)
        InputParam.addWidget(self.inputSignalEdit)
        
        # Vertical layouts
        layoutV = QVBoxLayout()
        layoutV.addLayout(TimeLayout)
        layoutV.addLayout(SysParam)
        layoutV.addLayout(States)
        layoutV.addLayout(NoiseParam1)
        layoutV.addLayout(NoiseParam2)
        layoutV.addLayout(NoiseParam3)
        layoutV.addLayout(InputParam)
        layoutV.addLayout(Button)
        layoutV.addStretch(1)
        
        # Fianl GUI Layout        
        layout = QHBoxLayout()
        layout.addLayout(Graph, 5)
        layout.addLayout(layoutV, 1)
        self.widget.setLayout(layout)
        self.setCentralWidget(self.widget)

        # Signals:
        # Pressing Enter returns Output
        self.paramEdit.returnPressed.connect(self.runButton1)
        # Button clicks to run and clear
        self.b1.clicked.connect(self.runButton1)
        self.b2.clicked.connect(self.runButton2)
        self.b3.clicked.connect(self.clearPlot)
        self.b4.clicked.connect(self.zoom)
        self.b5.clicked.connect(self.pan)
        self.b6.clicked.connect(self.home)
    
    # Plot Toolbar functions
    def zoom(self):
        self.toolbar.zoom()
    
    def pan(self):
        self.toolbar.pan()
        
    def home(self):
        self.toolbar.home()

    def kalmanfilterInit(self, mode=1):
        """
        """
        m = eval(self.massEdit.text())
        Ks = eval(self.KsEdit.text())
        Kd = eval(self.KdEdit.text())
        F1 = eval(self.FEdit1.text())
        F2 = eval(self.FEdit2.text())
        F3 = eval(self.FEdit3.text())
        F4 = eval(self.FEdit4.text())
        G1 = eval(self.GEdit1.text())
        G2 = eval(self.GEdit2.text())
        X0_0 = eval(self.ICEdit1.text())
        X0_1 = eval(self.ICEdit2.text())
        dt = eval(self.dtEdit.text())
        time = eval(self.simEdit.text())
        
        self.m = m
        self.Ks = Ks
        self.Kd = Kd
        self.F1 = F1
        self.F2 = F2
        self.F3 = F3
        self.F4 = F4
        self.G1 = G1
        self.G2 = G2
        self.X0 = [X0_0, X0_1]
        self.dt = dt
        self.time = time        
        
        F = np.array([[self.F1, self.F2], [self.F3, self.F4]])
        G = np.array([[self.G1], [self.G2]])
        H = np.array([1.0, 0.0])
        J = 0
        sysc = control.ss(F,G,H,J)
                
        # initialize simulations and Kalman Filter
        tf = max(self.time, floor(1000*self.dt))
        n = 0
        t = np.zeros((int(tf/self.dt),1))
        for i in np.arange(0,tf,self.dt):
            t[n] = i
            n += 1
        nSteps = len(t)
        Qc2 = eval(self.AccelNoiseEdit.text())
        Qd2 = Qc2/self.dt
        w = np.sqrt(Qd2)*np.random.randn(nSteps,1)
        R = eval(self.measNoiseEdit.text())
        v = np.sqrt(R)*np.random.randn(nSteps,1)
        sysd = control.c2d(sysc,self.dt)
        [Phi, Gamma, H, J] = control.ssdata(sysd)
        K = np.zeros((2,nSteps))
        xp = np.zeros((2,nSteps))
        Pp = np.identity(2)
        xp[:,0]=[0, 0]
        
        # initialize INS/GPS
        Fi = np.array([[0.0, 1.0],[0.0, 0.0]])
        Gi = np.array([[0],[1]])
        Hi = H*1
        Ji = 0
        Qdi2 = eval(self.AccelNoiseEdit.text())
        Qdi2 = Qdi2*Qdi2
        wi = np.sqrt(Qdi2)*np.random.randn(nSteps,1)
        Ri = R*1
        sysdi = control.c2d(control.ss(Fi,Gi,Hi,Ji),self.dt)
        [Phii, Gammai, Hi, Ji] = control.ssdata(sysdi)
        Ki = K*1
        xpi = xp*1
        Ppi = Pp*1
        
        u = eval(self.inputSignalEdit.text())
        u = 1*u.T
        yt, t, xt = matlab.lsim(sysc, u.T+w, t, self.X0) # simulate with input noise
        xt = xt.T
        y = yt+v
        
        accel = np.array([0, 1])@(F@xt + G*(u+w.T))
        accelmeas = accel + wi.T
        
        for k in range(1, nSteps):
            # MSD Filtering
            xm = Phi@xp[:,[(k-1)]] + Gamma*u[:,k-1]
            Pm = Phi@Pp@Phi.T + Gamma*Qd2*Gamma.T
            
            K[:,[k]] = Pm@H.T*np.linalg.inv(H@Pm@H.T + R) # Kalman gain K
            xp[:,[k]] = xm + K[:,[k]]*(y[0,k]-H*xm)
            Pp = Pm - K[:,[k]]@H@Pm
            
            # INS/GPS Filtering
            xmi = Phii@xpi[:,[(k-1)]] + Gammai*accelmeas[:, k-1]
            Pmi = Phii@Ppi@Phii.T + Gammai*Qdi2*Gammai.T
            
            Ki[:,[k]] = Pmi@Hi.T*np.linalg.inv(Hi@Pmi@Hi.T + Ri)    
            xpi[:,[k]] = xmi + Ki[:,[k]]*(y[0,k]-Hi*xmi)
            Ppi = Pmi - Ki[:,[k]]@Hi@Pmi
        
        if mode == 1: # returns position, estimate, and INS/GPS estimate
            return(t, xt[0,:], t, xp[0,:], t, xpi[0,:])
        elif mode == 2: # returns velocity
            return(t, xt[1,:], t, xp[1,:], t, xpi[1,:])
            
    def saveas(self):
        """Save input and output to a text file as seperate columns
        """
        name = QFileDialog.getSaveFileName(self, "saveas")[0]
        f = open(name, 'w')
        out = self.kalmanfilterInit()
        out = np.array(out)
        t = out[0]   # time
        xt = out[1]  # X
        xp = out[3]  # MSD Estimate
        xpi = out[5] # INS/GPS Estimate
        f.write(" Time     X    MSD   INS/GPS\n")
        for i in range(len(t)):
            f.write("%5.4f %5.4f %5.4f %5.4f\n" % (t[i], xt[i], xp[i], xpi[i]))
        f.close()
        
        self.plot.savePlot()
                
    def about(self):
        QMessageBox.about(self, 
            "About Function Evaluator",
            """<b>Function Evaluator</b>
               <p>Copyright &copy; 2016 Jeremy Roberts, All Rights Reserved.
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))
        
    def runButton1(self):
        kf = self.kalmanfilterInit()
        try:
            self.plot.redraw(kf[0], kf[1], kf[2], kf[3], kf[4], kf[5], 1)

        except:
            print('Input Error! Check the input')
        
    def runButton2(self):
        kf = self.kalmanfilterInit(mode=2)
        try:
            self.plot.redraw(kf[0], kf[1], kf[2], kf[3], kf[4], kf[5], 2)

        except:
            print('Input Error! Check the input')
        
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
        self.axes.set_title('Kalman Filter Simulation')
        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
    def redraw(self, x, y, u, v, i, j, mode=1): # plots position vs time
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.set_xlabel('Time (s)')
        if mode==1:
            self.axes.set_title('Position: MSD VS. INS/GPS Estimates')
            self.axes.plot(x, y, '--', label='True Position')
            self.axes.plot(u, v, '.', label='MSD Kalman Estimate')
            self.axes.plot(i, j, label='INS/GPS Estimate')
            self.axes.set_ylabel('Position (m)')
        elif mode==2:
            self.axes.set_title('Velocity: MSD VS. INS/GPS Estimates')
            self.axes.plot(x, y, '--', label='True Position')
            self.axes.plot(u, v, '.', label='MSD Kalman Estimate')
            self.axes.plot(i, j, label='INS/GPS Estimate')
            self.axes.set_ylabel('Velocity (m/s)')            
        self.axes.legend()
        self.draw()
        
    def clear(self):
        self.axes.clear()
        self.axes.set_title('Kalman Filter Simulation')
        self.axes.set_xlabel('Time (s)')
        self.draw()
        
    def savePlot(self):
        self.fig.savefig('result.png', bbpx_inches='tight')
        
app = QApplication(sys.argv)
form = MainWindow()
form.resize(1000, 300)
form.show()
app.exec_()
