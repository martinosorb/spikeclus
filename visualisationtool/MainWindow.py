import PyQt4
from PyQt4.QtCore import QThread
from WindowDesign import *
from FlagConfiguration import *
from mplwidget import *
import numpy as np
from PyQt4 import QtCore
from Controller import ActionController
from PyQt4.QtGui import QColor, QDialog, QRadioButton, QTableWidgetItem, QCheckBox, QPixmap, QLabel, QAbstractItemView, \
    QApplication, QTableView, QBrush, QItemSelectionModel
from PyQt4.QtCore import Qt
from DialogConf import Dialog
#import scipy.special._ufuncs_cxx
import time
import matplotlib as mpl
import os

#_fromUtf8 = QtCore.QString.fromUtf8


class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        # set some defaults (should be done here?)
        #mpl.rcParams['lines.linewidth'] = 1
        mpl.rcParams['font.size'] = 6
        mpl.rcParams['axes.titlesize'] = 7   # fontsize of the axes title
        mpl.rcParams['axes.labelsize'] = 7 # fontsize of the x any y labels

        QtGui.QWidget.__init__(self, parent)
        #_fromUtf8 = QtCore.QString.fromUtf8
        self.ui = Ui_MainWindow()
        self.cw = Ui_Dialog()
        self.ui.setupUi(self)
        self.diag = Dialog(self)
        self.cw.setupUi(self.diag)

        self.listLocClu = []
        self.listDataPlot = []
        self.listCl =[]
        self.dataCursor = None
        self.dataConfig = False
        self.listSpikePlot = []
        self.previousbutton = None
        self.abc = None

        self.textc = ""
        self.texts = ""
        self.hx = []
        self.selectedh = None
        self.selectedspk = None
        self.timesselected=None
        self.previousstartx = -1
        self.previousline = None
        self.linestart = None
        self.previouslineend = None
        self.previousendx = -1
        self.lineend = None
        self.listRowSelected2 = np.array([],dtype=int)
        self.listWaveSelected = []
        self.listSpikesSelected = []
        self.visiblePCA = True

        self.spikefile = ""
        self.lineth = None
        self.colorspike = 'red'
        self.colorpoint = 'blue'
        self.errorn = ''
        self.maxv = 0

        self.initbuttons()

        self.ct = ActionController()
        self.ui.mpl.ntb.setController(self.ct)
        self.ui.mpl.ntb.setMainWindow(self)

        self.ui.widget_5.ntb.setController(self.ct)
        self.ui.widget_5.ntb.setMainWindow(self)


        self.ct.setView(self)

        self.diag.setupConfWindow(self.cw, self, self.ct)
        self.thread = MyThread(self, self.ct)
        self.enableClusterControls(False)
        self.enablePCAControls(False)


        #Add listeners to all buttons
        QtCore.QObject.connect(self.ui.open_raw, QtCore.SIGNAL('clicked()'), self.select_file_raw)
        QtCore.QObject.connect(self.ui.push_plot, QtCore.SIGNAL('clicked()'), self.openfileraw)
        QtCore.QObject.connect(self.ui.run_cluster, QtCore.SIGNAL('clicked()'), self.run_cluster)
        QtCore.QObject.connect(self.ui.checkbutt, QtCore.SIGNAL('clicked()'), self.checkClusters)
        QtCore.QObject.connect(self.ui.uncheckbutt, QtCore.SIGNAL('clicked()'), self.uncheckClusters)
        QtCore.QObject.connect(self.ui.pushButton_2, QtCore.SIGNAL('clicked()'), self.filterTime)
        QtCore.QObject.connect(self.ui.pushButton_6, QtCore.SIGNAL('clicked()'), self.ct.goafterwards)
        QtCore.QObject.connect(self.ui.pushButton_5, QtCore.SIGNAL('clicked()'), self.ct.gobackwards)

        QtCore.QObject.connect(self.ui.pushButton_7, QtCore.SIGNAL('clicked()'), self.printPCA)
        QtCore.QObject.connect(self.ui.pushButton_3, QtCore.SIGNAL('clicked()'), self.flagClusters)
        QtCore.QObject.connect(self.ui.pushButton_11, QtCore.SIGNAL('clicked()'), self.unflagClusters)
        QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL('clicked()'), self.changecolor)
        QtCore.QObject.connect(self.ui.pushButton_4, QtCore.SIGNAL('clicked()'), self.hidewindow)
        QtCore.QObject.connect(self.ui.pushButton_8, QtCore.SIGNAL('clicked()'), self.showconfwindow)
        QtCore.QObject.connect(self.ui.pushButton_10, QtCore.SIGNAL('clicked()'), self.alltime)
        QtCore.QObject.connect(self.ui.pushButton_9, QtCore.SIGNAL('clicked()'), self.saveFlagged)
        QtCore.QObject.connect(self.ui.pushButton_12, QtCore.SIGNAL('clicked()'), self.zoomcluster)

        #Add listeners to radio buttons
        self.ui.radioButton_2.toggled.connect(self.clusterview)
        self.ui.radioButton.toggled.connect(self.spikesview)
        self.ui.radioButton_5.toggled.connect(self.completepcaview)
        self.ui.radioButton_3.toggled.connect(self.linearview)


        #Called when the PCA computation has finished
        self.connect(self.thread, QtCore.SIGNAL("finished()"), self.updateUi)

        self.ui.widget_12.canvas.mpl_connect('button_press_event', self.onpick3)
        self.ui.widget_9.canvas.mpl_connect('pick_event', self.onpickwave)

        #Add listeners to the edit (start and end) lines.
        self.ui.line_start.editingFinished.connect(self.miEvent)
        self.ui.line_end.editingFinished.connect(self.miEvent2)
        self.ui.lineEdit.editingFinished.connect(self.changeBinsEvent)
        self.ui.line_start.mousePressEvent = self.linestartevent
        self.ui.line_end.mousePressEvent = self.lineendevent

        #Add listeners to the label(start and end)lines
        #The start line is filled out with the start time of the recording session
        #The end line is filled out with the end time of the recording session
        self.ui.label_2.mousePressEvent = self.fillstartline
        self.ui.label_3.mousePressEvent = self.fillendline

        self.configureAxesPCA()
        self.configureISIPlot()
        self.configureWaveformPlot()
        self.configureSpikeTrainPlot()

        self.ui.pushButton.setStyleSheet("QWidget { background-color : blue}")
        self.ui.pushButton_9.setToolTip("Saves the flagged clusters in the current directory")


    def configureSpatialmapPlot(self):
        #Initializes the axes information of the Spatial map
        self.ui.mpl.canvas.ax.set_xlabel('X (Electrode ID)', fontsize=10)
        self.ui.mpl.canvas.ax.set_ylabel('Y (Electrode ID)', fontsize=10)

    def configureISIPlot(self):
        #Initializes the axes information of the ISI window
        self.ui.widget_12.canvas.ax.set_title('ISI Histogram',fontsize=10)
        self.ui.widget_12.canvas.ax.set_xscale('log')
        self.ui.widget_12.canvas.ax.set_xlabel('s', fontsize=10)
        self.ui.widget_12.canvas.fig.tight_layout(rect=[0.15, 0.35, .9, .85])

    def configureSpikeTrainPlot(self):
        #Initializes the axes information of the Spike train window
        self.ui.widget_11.canvas.ax.set_title('Spike train',fontsize=10)
        self.ui.widget_11.canvas.ax.get_yaxis().set_visible(False)
        self.ui.widget_11.canvas.fig.tight_layout(rect=[0.01, 0.4, .99, .8])
        self.ui.widget_11.canvas.ax.set_xlabel('s', fontsize=10)

    def configureWaveformPlot(self):
        #Initializes the axes information of the Waveform window
        self.ui.widget_9.canvas.ax.set_title('Waveforms', fontsize=10 )
        self.ui.widget_9.canvas.ax.set_ylabel('mV', fontsize=10)
        self.ui.widget_9.canvas.ax.set_xlabel('ms', fontsize=10)
        self.ui.widget_9.canvas.fig.tight_layout(rect=[0.15, 0.17, .9, .85])

    def initbuttons(self):
        self.ui.checkbutt.setDisabled(True)
        self.ui.uncheckbutt.setDisabled(True)
        self.ui.pushButton_3.setDisabled(True)
        self.ui.pushButton_2.setDisabled(True)
        self.ui.pushButton_10.setDisabled(True)
        self.ui.pushButton_8.setDisabled(True)
        self.ui.pushButton_9.setDisabled(True)
        self.ui.pushButton_5.setDisabled(True)
        self.ui.pushButton_6.setDisabled(True)
        self.ui.pushButton_7.setDisabled(True)
        self.ui.pushButton_11.setDisabled(True)
        self.ui.pushButton_12.setDisabled(True)
        self.ui.line_start.setDisabled(True)
        self.ui.line_end.setDisabled(True)

    def zoomcluster(self):
        #Zooms into the current selected cluster in the Spatial map
        if self.ct.getCurrentCluster() is not None:
            minx,maxx,miny,maxy = self.ct.getClusterInfo()
            self.ui.mpl.ntb.zoom_point(minx, maxx, miny, maxy)
        else:
            self.printLog("You haven't selected any clusters yet")

    def saveFlagged(self):
        #Saves the flagged clusters in the file flagged.txt
        f=open('flagged.txt','w')
        for i in self.ct.getFlaggedList():
            f.write(str(i) + "\n")
        f.close()

    def alltime(self):
        #Selects all the complete hdf5 recording session
        self.printStartLine(self.ct.mint)
        # make sure data is in range (hack, not nice)
        self.ui.line_start.setText(str(self.ct.mint/self.ct.sampling+0.0001))
        self.printEndLine(self.ct.maxt)
        self.ui.line_end.setText(str(self.ct.maxt/self.ct.sampling-0.0001))
        self.ui.pushButton_2.setFocus()



    def flagClusters(self):
        #Flags all clusters currently visible in the Cluster table
        self.ct.flagClusters(self.ct.getFilteredRows())

    def unflagClusters(self):
        #Unflags all clusters currently visible in the Cluster table
        self.ct.unflagClusters(self.ct.getFilteredRows())


    def getConfWindow(self):
        return self.diag

    def showconfwindow(self):
        #Shows the Configuration Window
        self.diag.exec_()


    def hidewindow(self):
        #Hides the Spike train window
        if not self.ui.pushButton_4.isChecked():
            self.ui.widget_11.setVisible(True)
            self.ui.widget_hold.setVisible(True)
            self.ui.pushButton_4.setText("Hide time window")
        else:
            self.ui.widget_11.setVisible(False)
            self.ui.widget_hold.setVisible(False)
            self.ui.pushButton_4.setText("Show time window")
            self.ui.groupBox_2.update()
            self.ui.groupBox_2.repaint()

    def isNumber(self, x):
        isnum = True
        try:
            num = float(x)
        except ValueError:
            isnum = False
        return isnum


    def linestartevent(self, event):
        #When the start edit line is pressed the histogram will start listening for pressing button events.
        self.linestart = self.ui.widget_4.canvas.mpl_connect('button_press_event', self.handleStart)

    def lineendevent(self, event):
        #When the end edit line is pressed the histogram will start listening for pressing button events.
        self.lineend = self.ui.widget_4.canvas.mpl_connect('button_press_event', self.handleEnd)


    def fillstartline(self, event):
        #Fills out the start line edit with the start time of the recording session.
        self.printStartLine(self.ct.mint)
        self.ui.line_start.setText(str(self.ct.mint/self.ct.sampling))


    def fillendline(self, event):
        #Fills out the end line edit with the end time of the recording session.
        self.printEndLine(self.ct.maxt)
        self.ui.line_end.setText(str(self.ct.maxt/self.ct.sampling))


    def miEvent(self):
        #When the user presses the enter button in the start line edit, the application will draw a green line in
        #the histogram
        if self.checkingData(self.ui.line_start):
            x = float(self.ui.line_start.text())*self.ct.getSampling()
        else:
            return

        self.printStartLine(x)
        self.ui.widget_4.canvas.mpl_disconnect(self.linestart)
        self.linestart = None


    def miEvent2(self):
        #When the user presses the enter button in the end line edit, the application will draw a red line in
        #the histogram
        if self.checkingData(self.ui.line_end):
            x = float(self.ui.line_end.text())*self.ct.getSampling()
        else:
            return

        self.printEndLine(x)
        self.ui.widget_4.canvas.mpl_disconnect(self.lineend)
        self.lineend = None



    def printStartLine(self, x):
        #Prints a green line in the time histogram at the position specified by x
        if self.previousstartx != x:
            self.previousstartx = x
            if self.previousline is not None:
                self.previousline.remove()
            self.previousline = self.ui.widget_4.canvas.ax.vlines(x,[0],self.ct.getMaxTimeHistogram() + 10,lw=1, colors='green')
            self.ui.widget_4.canvas.draw()

    def printEndLine(self, x):
        #Prints a red line in the time histogram at the position specified by x
        if self.previousendx != x:
            self.previousendx = x
            if self.previouslineend is not None:
                self.previouslineend.remove()
            self.previouslineend = self.ui.widget_4.canvas.ax.vlines(x,[0],self.ct.getMaxTimeHistogram() + 10,lw=1, colors='red')
            self.ui.widget_4.canvas.draw()


    def handleStart(self,event):
        #When the user clicks on the time histogram and previously clicked on the start line, then it will print
        #a green line on the plot.
        if self.ui.widget_4.ntb._active is None:
            if self.previousline is not None:
                self.previousline.remove()
                self.previousline = None


            self.printStartLine(event.xdata)
            self.ui.line_start.setText(str(event.xdata/self.ct.sampling))


    def handleEnd(self,event):
        #When the user clicks on the time histogram and previously clicked on the end line, then it will print
        #a red line on the plot.
        if self.ui.widget_4.ntb._active is None:
            if self.previouslineend is not None:
                self.previouslineend.remove()
                self.previouslineend = None


            self.printEndLine(event.xdata)
            self.ui.line_end.setText(str(event.xdata/self.ct.sampling))



    def checkingData(self, textedit):
        #Validates that the value given in either edit lines is valid
        if self.isNumber(textedit.text()) & (textedit.text() != ""):
            x = float(textedit.text())*self.ct.getSampling()
        else:
            self.printLog("A numeric value must be entered")
            textedit.setStyleSheet("border: 1px solid red")
            return False


        if (x>self.ct.getMaxTime()) | (x<self.ct.getMinTime()):
            minv = self.ct.getMinTime()/self.ct.getSampling()
            maxv = self.ct.getMaxTime()/self.ct.getSampling()
            self.printLog("Values must be between[" + str(minv) + "," + str(maxv) + "]")
            textedit.setStyleSheet("border: 1px solid red")
            return False

        textedit.setStyleSheet("")
        return True



    def checkingDataFinal(self, textedit1, textedit2):
        #Validates the start and end time values
        if self.linestart is not None:
            self.ui.widget_4.canvas.mpl_disconnect(self.linestart)
            self.linestart = None
        if self.lineend is not None:
            self.ui.widget_4.canvas.mpl_disconnect(self.lineend)
            self.lineend = None

        if self.checkingData(textedit1):
            x1 = float(textedit1.text())*self.ct.getSampling()
        else:
            return False

        if self.checkingData(textedit2):
            x2 = float(textedit2.text())*self.ct.getSampling()
        else:
            return False

        if x2 < x1:
            self.printLog("Start time must be less than end time")
            return False

        return True


    def changecolor(self):
        #Changes the current spike color
        colors = ['QWidget { background-color : blue}', 'QWidget { background-color : red}',
                  'QWidget { background-color : green}', 'QWidget { background-color : black}',
                  'QWidget { background-color : yellow}']
        colors2 = ['blue', 'red', 'green', 'black', 'yellow']
        c= str(self.ui.pushButton.styleSheet())
        ind = colors.index(c)
        ind = ind +1
        if ind == len(colors):
            ind = 0
        self.colorpoint = colors2[ind]
        self.ui.pushButton.setStyleSheet(colors[ind])
        self.ct.updateSpikeColor()

    def getSpikeColor(self):
        #Gets the current spike color
        return self.colorpoint


    def enableClusterControls(self, mode):
        #Enables/Disables the cluster/spike view radio buttons in the Cluster map
        self.ui.radioButton_2.setEnabled(mode)
        self.ui.radioButton.setEnabled(mode)

    def enablePCAControls(self,mode):
        #Enables/Disables the complete/filtered view radio buttons in the PCA window
        self.ui.radioButton_5.setEnabled(mode)
        self.ui.radioButton_6.setEnabled(mode)



    def configureAxesPCA(self):
        self.ui.widget_5.canvas.ax.set_xticks((-4, 0,4))
        self.ui.widget_5.canvas.fig.tight_layout(rect=[0.15, 0.17, .9, .85])
        self.ui.widget_5.canvas.ax.set_title('PCA, Spatial Clusters',fontsize=10)
        self.ui.widget_5.canvas.ax.set_xlabel('PC 1', fontsize=10)
        self.ui.widget_5.canvas.ax.set_ylabel('PC 2', fontsize=10)


    def onpickwave(self, event):
        #Changes the waveform color of the selected spike
        colors = ['red', 'blue', 'green', 'black', 'yellow']
        c= event.artist.get_color()
        ind = colors.index(c)
        ind = ind +1
        if ind == len(colors):
            ind = 0
        self.colorspike = colors[ind]
        event.artist.set_color(self.colorspike)
        self.ui.widget_9.canvas.draw()


    def getScale(self):
        #Returns the current scale of the ISI histogram
        if self.ui.radioButton_3.isChecked():
            return 'linear'
        else:
            return 'log'

    def linearview(self, enabled):
        self.ct.handleChangeHistogram()

    def changeBinsEvent(self):
        self.ct.handleChangeHistogram()

    def getBinsNumber(self):
        return int(self.ui.lineEdit.text())

    def printPCA(self):
        #Lauches the thread associated to the PCA computation
        self.ui.textBrowser.setText("Computing PCA started...\n")
        self.thread.start()

    def updateUi(self):
        #When the PCA computation has finished, it receives such signal to print in the PCA window the clusters that
        #are in hold on state or selected
        self.ui.textBrowser.setText("Computing PCA finished...\n")
        self.ui.pushButton_7.setDisabled(True)
        self.ct.printPCAClusters()





    def getTimeStart(self):
        return self.timestart

    def getTimeEnd(self):
        return self.timeend


    def getTimeStartSec(self):
        return float(self.ui.line_start.text())

    def getTimeEndSec(self):
        return float(self.ui.line_end.text())


    def filterTime(self):
        #Starts filtering the clusters that are within the start and end times given by the user

        if not self.checkingDataFinal(self.ui.line_start, self.ui.line_end):
            return

        self.timestart = (float(self.ui.line_start.text())*self.ct.getSampling())
        self.timeend = (float(self.ui.line_end.text())*self.ct.getSampling())


        self.printLog("")
        self.printLog("Filtering started...\n")
        self.ct.filterIndexes(self.timestart, self.timeend)
        self.printLog("Filtering finished...\n")
        self.ui.pushButton_2.setDisabled(True)
        self.ui.pushButton_10.setDisabled(True)
        self.ui.pushButton_7.setDisabled(False)

        self.disableFilterControls(False)


    def completepcaview(self, enabled):
        if enabled:
            self.ui.widget_5.ntb.desactiveToolbar()
            self.ct.completepcaviewcontrol()
        else:
            self.ui.widget_5.ntb.exploratoryPCAMode()
            self.ct.filteredpcaviewcontrol()

        self.ui.widget_5.canvas.draw()






    def spikesview(self, enabled):
        if enabled:
            self.ct.spikesviewcontrol()
            self.ui.tableView.setDisabled(False)

            if self.ui.lineEdit_2.text() != "":
                self.textc = self.ui.lineEdit_2.text()
                self.ui.lineEdit_2.setText("")
            if self.texts != "":
                self.ui.lineEdit_2.setText(self.texts)

            self.ui.mpl.ntb.setEnableCursorCluster(False)
            self.ui.mpl.ntb.setEnableClusterShape(False)
            self.ui.mpl.ntb.setEnableSelectedCluster(False)
            self.ui.mpl.ntb.setEnableHeatCluster(False)

            self.ui.mpl.ntb.exploratorySpikeMode()

            self.ui.mpl.canvas.draw()


    def heatview(self, enabled):
        if enabled:
            self.ct.heatviewcontrol()
            #self.ui.tableView.setDisabled(True)
            #self.ui.radioButton.setDisabled(True)
            #self.ui.radioButton_2.setDisabled(True)

            if self.ui.lineEdit_2.text() != "":
                self.textc = self.ui.lineEdit_2.text()
                self.ui.lineEdit_2.setText("")

            #self.ui.mpl.ntb.setEnableCursorCluster(False)
            #self.ui.mpl.ntb.setEnableClusterShape(False)
            #self.ui.mpl.ntb.setEnableSelectedCluster(False)
            self.ui.mpl.ntb.setEnableHeatCluster(True)

            self.ui.mpl.canvas.draw()



    def exploratoryMode(self):
        self.ui.mpl.ntb.exploratoryMode()

    def exploratoryPCAMode(self):
        self.ui.widget_5.ntb.exploratoryPCAMode()

    def clusterview(self, enabled):
        if enabled:

            self.ct.clusterviewcontrol()
            self.ui.tableView.setDisabled(False)
            self.ui.radioButton.setDisabled(False)
            self.ui.radioButton_2.setDisabled(False)

            if self.ui.lineEdit_2.text() != "":
                self.texts = self.ui.lineEdit_2.text()
                self.ui.lineEdit_2.setText("")
            if self.textc != "":
                self.ui.lineEdit_2.setText(self.textc)

            self.ui.mpl.ntb.setEnableCursorCluster(True)
            self.ui.mpl.ntb.setEnableClusterShape(True)
            self.ui.mpl.ntb.setEnableSelectedCluster(True)
            self.ui.mpl.ntb.setEnableHeatCluster(True)


            if self.ui.mpl.ntb.isToolbarActive() == False:
                self.ui.mpl.ntb.exploratoryMode()

            self.ui.mpl.canvas.draw()

    def setPointVisible(self, ref, mode):
        if ref is not None:
            ref.set_visible(mode)

    def setPCAClusterVisible(self, ref, mode):
        if ref is not None:
            ref[0].set_visible(mode)

    def setRadialClustersVisible(self, ref, mode):
        if ref is not None:
            ref.set_visible(mode)

    def setClusterShapesVisible(self, ref, mode):
        if ref is not None:
            ref.set_visible(mode)

    def listhandler2(self, index):

        if (index.isValid()):
            i = index.row()
            c = index.column()

        if c == 0:
            model = self.ui.tableView.model()
            index = model.index(i,1)
            index2 = model.index(i,0)
            cl = int(model.data(index))
            ind = self.ct.getMapClusters()[cl]

            if (model.data(index2) == "U") & (self.ct.getLastClusterSel() != ind):
                model.setData(index2, 1)
            elif (model.data(index2) == "H") & (self.ct.getLastClusterSel() != ind):
                model.setData(index2, 0)
            elif (model.data(index2) == "U") & (self.ct.getLastClusterSel() == ind):
                model.setData(index2, 3)
            elif (model.data(index2) == "H") & (self.ct.getLastClusterSel() == ind):
                model.setData(index2, 2)


            self.ct.handleChange(ind, i)
            self.ct.filterExtraPoint(ind)





        if c in [1,2,3]:
            model = self.ui.tableView.model()
            index = model.index(i,1)
            cl = int(model.data(index)) #int(model.data(index))
            currentitem = self.ct.getMapClusters()[cl]
            self.ct.handleSelection(currentitem)






        if c == 4:
            model = self.ui.tableView.model()
            index = model.index(i,1)
            index2 = model.index(i,4)

            ind = int(model.data(index))
            ind = self.ct.getMapClusters()[ind]

            if model.data(index2) == "U":
                self.ct.flagClusters([ind])
                model.setData(index2, 1)
            else:
                self.ct.unflagClusters([ind])
                model.setData(index2, 0)

            self.ct.filterExtraPoint(ind)

    def handleTableSelection(self,row, select):
        model = self.ui.tableView.model()
        index = model.index(row,1)
        index2 = model.index(row,0)
        cl = int(model.data(index))
        ind = self.ct.getMapClusters()[cl]

        if (model.data(index2) == "U") & (not select):
            model.setData(index2, 0)
        elif (model.data(index2) == "H") & (not select):
            model.setData(index2, 1)
        elif (model.data(index2) == "U") & (select):
            model.setData(index2, 2)
        elif (model.data(index2) == "H") & (select):
            model.setData(index2, 3)



    def getZoomLimits(self):
        #Gets the axes limits of the Cluster map
        xmin =  self.ui.mpl.canvas.ax.get_xlim()[0]
        xmax =  self.ui.mpl.canvas.ax.get_xlim()[1]

        ymin =  self.ui.mpl.canvas.ax.get_ylim()[0]
        ymax =  self.ui.mpl.canvas.ax.get_ylim()[1]
        return xmin, xmax, ymin, ymax


    def scrollTo(self, index):
        #Scrolls the cluster table to the given index position
        ind = self.ui.tableView.model().index(index,0)
        self.ui.tableView.scrollTo(ind)
        self.ui.tableView.setFocus()
        self.ui.tableView.setCurrentIndex(ind)



    def getCheckedStatus(self, row):
        #Gets the hold on status of a row
        model = self.ui.tableView.model()
        index = model.index(row,0)
        data = model.data(index)
        if data == 'U':
            return Qt.Unchecked
        else:
            return Qt.Checked

    def changeCheckedStatus(self, row):
        #Changes the hold on status of a row, if it is hold on it will turn it off and the other way around

        model = self.ui.tableView.model()
        index = model.index(row,1)
        index2 = model.index(row,0)
        cl = int(model.data(index))
        ind = self.ct.getMapClusters()[cl]

        if (model.data(index2) == "U") & (self.ct.getLastClusterSel() != ind):
            model.setData(index2, 1)
        elif (model.data(index2) == "H") & (self.ct.getLastClusterSel() != ind):
            model.setData(index2, 0)
        elif (model.data(index2) == "U") & (self.ct.getLastClusterSel() == ind):
            model.setData(index2, 3)
        elif (model.data(index2) == "H") & (self.ct.getLastClusterSel() == ind):
            model.setData(index2, 2)


    def getFlaggedStatus(self, row):
        model = self.ui.tableView.model()
        index = model.index(row,4)
        data = model.data(index)
        if data == 'U':
            return Qt.Unchecked
        else:
            return Qt.Checked


    def setFlaggedStatus(self, row, status):
        model = self.ui.tableView.model()
        index = model.index(row,4)
        if status == Qt.Checked:
            model.setData(index, 1)
        else:
            model.setData(index, 0)


    def setCheckedStatus(self, row, status):
        model = self.ui.tableView.model()
        index = model.index(row,0)
        if status == Qt.Checked:
            model.setData(index, 1)
        else:
            model.setData(index, 0)


    def listAllClusters(self, clusters, countspikes):
        #Creates the cluster table

        self.tabledata = []
        timedif = (self.getTimeEnd()-self.getTimeStart())/self.ct.getSampling()

        self.ui.tableView.clicked.connect(self.listhandler2)
        # set the table model
        header = ['Hold', 'Cluster', 'Spikes/s', 'PCA', 'Flag']


        for c in range(len(clusters)):

            row = ["U", str(clusters[c]), "{:.3f}".format(countspikes[c]), "", "U" ]
            self.tabledata.append(row)


        self.tm = MyTableModel(self.tabledata, header, self)


        self.tm.setColors(self.ct.getColours())

        self.proxyModel = MySortFilterProxyModel()
        self.proxyModel.setValidRows([True]*(len(clusters)))
        self.proxyModel.setSortRole(Qt.UserRole)
        self.proxyModel.setHHeader(header)

        self.proxyModel.setSourceModel(self.tm)
        self.proxyModel.setSource(self.tm)
        self.ui.tableView.setModel(self.proxyModel)
        self.ui.tableView.resizeColumnsToContents()

        self.ui.tableView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.ui.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.ui.tableView.setSortingEnabled(True)
        self.ui.tableView.horizontalHeader().sortIndicatorChanged.connect(self.ct.indicator)
        self.ui.tableView.keyPressEvent = self.keyPress


    def keyPress(self,e):
        #Manages cluster selection by clicking on the up and down arrows
        index = self.ui.tableView.selectionModel().currentIndex()
        row = index.row()

        index = self.ui.tableView.model().index(row, 3)
        self.ui.tableView.selectionModel().select(index, QItemSelectionModel.Select);

        column = -1

        if e.key() == QtCore.Qt.Key_Up:
            row = row - 1
            column = 1
        elif e.key() == QtCore.Qt.Key_Down:
            row = row + 1
            column = 1
        elif ((e.key() == QtCore.Qt.Key_F) and (QApplication.keyboardModifiers() & QtCore.Qt.ControlModifier)):
            column = 4
        elif ((e.key() == QtCore.Qt.Key_X) and (QApplication.keyboardModifiers() & QtCore.Qt.ControlModifier)):
            column = 0
        elif ((e.key() == QtCore.Qt.Key_S) and (QApplication.keyboardModifiers() & QtCore.Qt.ControlModifier)):
            column = 1
        else:
            return

        if row == -1:
            row = self.ct.gettableLength() - 1
        elif row == self.ct.gettableLength():
            row = 0


        self.scrollTo(row)
        ind = self.ui.tableView.model().index(row,column)
        self.listhandler2(ind)


    def setSpikeTrainLims(self, xmin, xmax):
        self.ui.widget_11.canvas.ax.set_xlim([xmin, xmax])


    def herenothing(self, other):
        return False


    def onpick3(self, event):
        self.ct.handleHistogramPick(event.xdata)




    def select_file_raw(self):
        self.spikefile = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF5 data files (*.hdf5)")
        if self.spikefile:
            self.ui.line_raw.setText(self.spikefile)
            self.ui.line_raw.setToolTip(self.spikefile)

    def openfileraw(self):
        #self.spikefile = 'C:\\Users\\user\\Documents\\MsCfiles\\data\\P29_16_05_14_retina02_left_stim2_smallarray_fullfield_SpkD45_v18_clustered.hdf5'


        #self.spikefile = 'C:\\Users\\user\\Documents\\MsCfiles\\data\\P44_04_02_15_ret2_right_stim1_class1_SpkD45_v18_clustered_0.3_0.16.hdf5'


        try:
            if self.spikefile != "":
                clustered = self.ct.loadData(str(self.spikefile))
                if clustered == "false":
                    self.printLog("Try press the Cluster button")
                else:
                    self.ui.line_start.setDisabled(False)
                    self.ui.line_end.setDisabled(False)
                    self.ui.pushButton_2.setDisabled(False)
                    self.ui.pushButton_10.setDisabled(False)
                    self.printLog("Data file has been loaded correctly")
            else:
                self.printLog("An .hdf5 file must be provided using the Open File dialog.")
        except ValueError:
            self.printLog(str(ValueError.message))

    def run_cluster(self):
        self.printLog("Clustering")
        result, fileName = self.ct.runCluster(str(self.spikefile))
        if result and fileName: 
            self.spikefile = fileName
            self.printLog(result)
        else:
            self.printLog("Cluster unsuccessful")

    def uncheckClusters(self):
        self.ct.uncheckClusters(self.ct.getFilteredRows())

    def checkClusters(self):
        self.ct.checkClusters(self.ct.getFilteredRows())



    def printInfo(self, text):
        self.ui.lineEdit_2.setText(text)




    def setAxisTimeWidget(self, minx, maxx, miny, maxy):
        self.ui.widget_4.canvas.ax.set_xlim([minx, maxx])
        self.ui.widget_4.canvas.ax.set_ylim([miny, maxy])

    def setTicksTimeWidget(self, xticks, yticks, xlabels = None, ylabels = None):
        self.ui.widget_4.canvas.ax.set_xticks(xticks)
        self.ui.widget_4.canvas.ax.set_yticks(yticks)

        if xlabels is not None:
            self.ui.widget_4.canvas.ax.set_xticklabels(xlabels)
        if ylabels is not None:
            self.ui.widget_4.canvas.ax.set_yticklabels(ylabels)


    def plotTimeWidget(self, x, y):
        self.ui.widget_4.canvas.ax.plot(x,y)
        self.ui.widget_4.canvas.draw()


    def setAxisMainWidget(self, minx, maxx, miny, maxy):
        self.ui.mpl.canvas.ax.set_xlim(minx, maxx)
        self.ui.mpl.canvas.ax.set_ylim(miny, maxy)
        self.ui.mpl.canvas.ax.set_aspect('equal')

    def plotMainWidget(self,x,y, s,c ,marker,edgecolors,alpha = None,picker=None, zorder = None):
        return self.ui.mpl.canvas.ax.scatter(x,y,s = s,c=c,marker=marker,edgecolors=edgecolors,alpha=alpha,picker=picker, zorder = zorder)


    def drawMainWidget(self):
        self.ui.mpl.canvas.draw()

    def drawPCAWidget(self):
        self.ui.widget_5.canvas.draw()

    def drawWaveFormWidget(self):
        self.ui.widget_9.canvas.draw()

    def plotWaveComplete(self, shapes, color, index):
        plInds=range(np.min([30,shapes.shape[1]]))
        sl = 1.0*np.shape(shapes)[0]

        xdata = np.arange(sl)/self.ct.getSampling()
        x = [self.ui.widget_9.canvas.ax.plot(xdata,shapes[:,i],color=color,alpha=0.2) for i in plInds]
        plot = self.ui.widget_9.canvas.ax.plot(xdata,np.mean(shapes,axis=1),color=color,lw=2.5, alpha = 1)


        labels = ["{:.2f}".format(xdata[int(i * (sl-1))]*1000) for i in [0, 0.2, 0.4, 0.6, 0.8,1]]
        ticks =  [xdata[int(i * (sl-1))] for i in [0, 0.2, 0.4, 0.6, 0.8,1]]


        self.ui.widget_9.canvas.ax.set_xticks(ticks)
        self.ui.widget_9.canvas.ax.set_xticklabels(labels)

        return x,plot

    def drawHistogramWidget(self):
        self.ui.widget_12.canvas.draw()
        #self.ui.widget_12.canvas.show()


    def setAxisHistogramScale(self, scaletype, minx, maxx, miny, maxy,  nbins):
        #self.ui.widget_12.canvas.ax.set_title('ISI Histogram',fontsize=10)

        self.ui.widget_12.canvas.ax.set_xscale(scaletype)

        ticks = [int(i * maxy) for i in [0, 0.25,0.5 ,0.75, 1]]
        self.ui.widget_12.canvas.ax.set_yticks(ticks)
        self.ui.widget_12.canvas.ax.set_ylim([0, maxy])

        if scaletype == 'linear':
            data = maxx-minx
            ticks = [(i * data) for i in [0, 0.25,0.5 ,0.75, 1]]
            labels = ["{:.2f}".format((i * data)) for i in [0, 0.25, 0.5, 0.75, 1]]
            self.ui.widget_12.canvas.ax.set_xlim([minx, maxx])
            self.ui.widget_12.canvas.ax.set_xticks(ticks)
            self.ui.widget_12.canvas.ax.set_xticklabels(labels)

        else:

            self.ui.widget_12.canvas.ax.set_xlim([min(nbins), max(nbins)])



    def plotHistogram(self, x ,nbins, scaletype):

        n, bins, hist = self.ui.widget_12.canvas.ax.hist(x= x, bins=nbins, normed=0, histtype='bar', color= 'blue', rwidth=0.9, picker = 30)
        self.setAxisHistogramScale(scaletype, min(x), max(x), 0, max(n), nbins)
        return n,bins,hist

    def drawLineSpikesWidget(self):
        self.ui.widget_11.canvas.draw()


    def removeHistogramPoints(self, ref):
        if ref is not None:
            ref.remove()
        return None

    def removeSelectedSpikes(self, ref):
        if ref is not None:
            ref.remove()
        return None

    def deleteHistogram(self, ref):
        if ref is not None:
            for i in ref:
                i.remove()
        return None

    def removeClusterVTLines(self, ref):
        if ref is not None:
            ref[0].remove()
        return None

    def deleteLine(self, art):

        if (art is not None):
            if (len(art[0].get_axes().lines) != 0):
                art.pop(0).remove()
        return None

    def deleteLines(self, arts):
        if (arts is not None):
            if(len(arts[0][0].get_axes().lines) != 0):
                for ind in range(len(arts)):
                    arts[ind].pop(0).remove()
        return None

    def deleteVerticalSpike(self, art):
        if (art is not None):
            art.remove()
        return None

    def deletePoint(self, ref):
        if (ref is not None):
            ref.remove()
        return None

    def deleteSpikes(self, ref):
        if ref is not None:
            ref.remove()
        return None

    def deleteProjection(self, ref):
        if ref is not None:
            ref.pop(0).remove()
        return None


    def printClusterProjection(self,x ,y ,c, picker = None):
        return self.ui.widget_5.canvas.ax.plot(x,y,'k.',c=c,ms=6,rasterized=True,lw=0, picker = picker)

    def printClusterPoint(self,x,y, color, edgecolor):
        return self.ui.mpl.canvas.ax.scatter(x,y, s = 25,c=color,marker='o',edgecolors=edgecolor, zorder =10000)

    def printSpikePoint(self, x, y, c, s =25, marker='+', edgecolors = None, alpha= None, picker=None):
        ref = self.ui.mpl.canvas.ax.scatter(x, y, s=s, c=c, marker = marker, edgecolors=edgecolors, picker = picker, zorder =10000)
        return ref


    def printProjectionPoint(self, x, y, c, s =25, marker='o', alpha= None, edgecolors = None, picker=None):
        return self.ui.widget_5.canvas.ax.scatter(x, y, s=s, c=c, marker = marker, edgecolors = edgecolors, picker = picker, zorder =10000)



    def addPCAClusterInfo(self, i, text):
        model = self.ui.tableView.model()
        index = model.index(i,3)
        model.setData(index, "{:.3f}".format(text))


    def setPCAValid(self):
        self.tm.setPCAValid()



    def setCurrentRowList1(self, index):
        return index

    def getSelectedCluster(self, index):
        model = self.ui.tableView.model()
        index = model.index(index,1)
        cl = int(model.data(index))
        return cl






    def setRowHidden(self,i, boolean):
        if boolean:
            rows = self.proxyModel.getValidRows()
            rows[i] = False
            self.proxyModel.setValidRows(rows)
            self.proxyModel.invalidateFilter()
        else:
            rows = self.proxyModel.getValidRows()
            rows[i] = True
            self.proxyModel.setValidRows(rows)
            self.proxyModel.invalidateFilter()



    def deselectList1(self):
        self.ui.tableView.clearSelection()

    def confWaveFormWidget(self):
        self.ui.widget_9.canvas.ax.grid(which='major', axis='x', linewidth=1, linestyle='-', color='0.75')
        self.ui.widget_9.canvas.ax.grid(which='major', axis='y', linewidth=1, linestyle='-', color='0.75')
        #self.ui.mpl.canvas.ax.grid(which='major', axis='x', linewidth=1, linestyle='-', color='0.75')
        #self.ui.mpl.canvas.ax.grid(which='major', axis='y', linewidth=1, linestyle='-', color='0.75')

    def printVerticalSpikes(self, t):

        if len(t) != 0:

            nbins,tb = np.histogram(t,1000)
            tbs = np.ceil(tb[2]-tb[1])
            self.maxv = max(nbins)
            return self.ui.widget_11.canvas.ax.plot(tb[1:]-tbs/2,nbins, color = 'black')

    def printVerticalSelectedSpikes(self, t):
        return self.ui.widget_11.canvas.ax.vlines(t,[0],[self.maxv], colors='red')

    def printVerticalSpike(self,t):
        return self.ui.widget_11.canvas.ax.vlines(t,[0],[self.maxv], colors='magenta', zorder =10)

    def isClusterView(self):
        return self.ui.radioButton_2.isChecked()

    def setPointsVisible(self, ref, mode):
        if ref is not None:
            ref.set_visible(mode)

    def printSpikeWave(self, shape):
        sl = 1.0*len(shape)
        a = self.ui.widget_9.canvas.ax.plot(np.arange(sl)/self.ct.getSampling(),shape,color=self.colorspike,alpha=0.9, zorder =1000, picker = 5)
        return a

    def printLog(self, text):
        self.ui.textBrowser.setText(text)

    def setTextItem(self, item, text):
        item.setText(text)

    def disableSpikeControls(self, mode):
        self.ui.pushButton_5.setDisabled(mode)
        self.ui.pushButton_6.setDisabled(mode)


    def disableFilterControls(self, mode):
        self.ui.pushButton_8.setDisabled(mode)
        self.ui.pushButton_9.setDisabled(mode)
        self.ui.checkbutt.setDisabled(mode)
        self.ui.uncheckbutt.setDisabled(mode)
        self.ui.pushButton_3.setDisabled(mode)
        self.ui.pushButton_11.setDisabled(mode)
        self.ui.pushButton_12.setDisabled(mode)


    def isCompletePCAView(self):
        return self.ui.radioButton_5.isChecked()

class MyThread(QThread):

    def __init__(self, mw, ct, parent=None):
        ''' Constructor. '''
        QThread.__init__(self, parent)
        self.mw = mw
        self.ct = ct

    def run(self):
        self.ct.computePCA()
