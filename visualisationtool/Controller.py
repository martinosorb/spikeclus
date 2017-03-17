from collections import Counter
from PyQt4.QtCore import Qt
import numpy as np
import time
from DataBase import DataBase
from matplotlib.backends.qt_compat import QtCore
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
#from matplotlib.collections import LineCollection
#import scipy.special._ufuncs_cxx
#import sklearn.utils.lgamma
#import sklearn.utils.weight_vector
import math
from numpy.linalg import norm
#from PyQt4.QtGui import QColor
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon, Ellipse, CirclePolygon #, Circle
from matplotlib.collections import PatchCollection
import os
import sys
CLUSTERING_PATH = os.getcwd() + '/../clustering'
sys.path.insert(0, CLUSTERING_PATH)
from herdingspikes import spikeclass, ImportInterpolated


class ActionController():


    def __init__(self, parent=None):
        self.db = None
        self.clf = None
        self.abc = None
        self.times = None
        self.centres = None
        self.clusterID = None
        self.data = None
        self.sampling = None
        self.shapes = None
        self.nClusters = None
        #self._points = None
        self._points = np.empty([0, 4], dtype=float)
        self.pointspca = np.empty([0, 4], dtype=float)
        self.pointsclust = None
        self.lastClusterSelected = None
        self.clusterpoint = None
        self.swave = None
        self.wave = None
        self.n = None
        self.bins = None
        self.patchesh = None
        self.ii = None
        self.timesh = None
        self.selectedh= None
        self.selectedspk = None
        self.vspikes = None
        self.spikepoint = None
        self.spikewave = None
        self.vspike = None
        self.lastSpikeSelected = None
        self.listRowSelected = np.array([],dtype=int)
        self.listSelected = []
        self.mw = None
        self.listSaved = []
        self.lastxleftbinSelected = -1
        self.previousbutton = None
        self.listWaveSelected = []
        self.listSpikesSelected = []
        self.listProjecSelected = []
        self.listProjecFiltSelected = []
        self.nbins = 0
        self.checkAll = False
        self.fit = None
        self.clustproj = None
        self.zord = 10
        self.spikeproj = None
        self.projfilt = None
        self.listSpikePlot = []
        self.spikespoint = None
        self.countspikes = np.array([],dtype=float)
        self.meanpca = np.array([],dtype=int)
        self.filtered = np.array([],dtype=int)
        self.flagged= np.array([],dtype=int)
        self.mint = None
        self.maxt = None
        self.p = None
        self.patches = []
        self.selpoints = np.empty([0, 4], dtype=float)
        self.listselpoints = []
        self.listselpatches =[]
        self.des = False
        self.des2 = False
        self.offsets = None
        self.circles = []
        self.now = 0
        self.past = -1
        self.last = 0
        self.lasti = -1
        self.lastsp = 0
        self.lastspi = -1
        self.lastcp = None


    def find_first(self, item, vec):
        """return the index of the first occurence of item in vec"""
        for i in xrange(len(vec)):
            if item == vec[i]:
                return i
        return -1

    def getPCAList(self):
        return self.fit



    def isPCALimits(self, index, mint, maxt):
        noise = False
        if self.fit is not None:
            mean = self.meanpca[index]
            if (mean > mint) & (mean < maxt):
                noise = True
        return noise

    def setView(self, mw):
        self.mw = mw



    def loadData(self, spikeFile):

        self.db = DataBase()

        clustered = self.db.setupDatabase(spikeFile)
        if clustered == "false":
            return "false"
        else:
            self.times = self.db.getTimes()
            self.centres = self.db.getCentres()
            self.clusterID = self.db.getClusterID()
            self.data = self.db.getData()
            self.sampling = self.db.getSampling()
            self.shapes = self.db.getShapes()
            self.nClusters = np.shape(self.centres)[1]

            self.setColours()

            self.setTimeHistogram()

            self.mw.confWaveFormWidget()

    def runCluster(self, spikeFile, h=0.3, alpha=0.28):
        Folder, Filename = os.path.split(spikeFile)
        O = ImportInterpolated(spikeFile) 
        scorePCA = O.ShapePCA(ncomp=2, white=True)
        O.CombinedMeanShift(h, alpha, PrincComp=scorePCA, mbf=10)
        result = 'Found ' + str(O.NClusters()) + ' clusters for ' + str(O.NData()) + ' spikes.'
        fileName = Folder + '/' + Filename.replace('.hdf5', '_clustered_' + str(h) + '_' + str(alpha) + '.hdf5')
        O.Save(fileName)
        return result, fileName

    def setColours(self):
        np.random.seed(1)
        coloursl=plt.cm.spectral(np.random.permutation(np.linspace(0,1,self.nClusters)))
        coloursl=np.append(np.array([0,0,0,0.5]),coloursl[:-1])
        self.colours=np.reshape(coloursl,(self.nClusters,4))



    def getColour(self, cluster):
        return self.colours[cluster]

    def getColours(self):
        return self.colours

    def indicator(self):
        self.sortclf = np.array(range(max(self.clf)+1))

        for i in range(self.gettableLength()):
            cluster = self.mw.getSelectedCluster(i)
            self.sortclf[cluster] = i



    def filterIndexes(self, start, end):

        timedif = (self.mw.getTimeEnd()-self.mw.getTimeStart())/self.sampling

        self.mw.setSpikeTrainLims(self.mw.getTimeStartSec() -4 , self.mw.getTimeEndSec() + 4)

        self.sInds = (self.times>=start)&(self.times<=end)
        cl = self.clusterID[self.sInds]
        self.clf = np.unique(cl)


        nana2 = np.argsort(self.clf)
        self.clf = self.clf[nana2]




        self.tablelengh = len(self.clf)
        self.sortclf = np.array(range(max(self.clf)+1))

        for i in range(len(self.clf)):
            self.sortclf[self.clf[i]] = i


        self.nData = np.shape(self.data)[1]



        self.xmin =  self.mw.ui.mpl.canvas.ax.get_xlim()[0]
        self.xmax =  self.mw.ui.mpl.canvas.ax.get_xlim()[1]

        self.ymin =  self.mw.ui.mpl.canvas.ax.get_ylim()[0]
        self.ymax =  self.mw.ui.mpl.canvas.ax.get_ylim()[1]

        self.clusterIDF = self.clusterID[self.sInds]
        self.clusterIDF = np.append(self.clusterIDF,100000)

        self.countspikes = Counter(self.clusterIDF)


        self.nana = np.argsort(self.clusterIDF)
        self.clusterIDF = np.append(self.clusterIDF[self.nana],100000)
        self.clf2 = np.append(self.clf,100000)





        self.printMapAll()
        self.mw.listAllClusters(self.clf, self.countspikes)
        self.filtered = np.array(range(len(self.getClusterList())),dtype=int)



        sl = 1.0*np.shape(self.shapes[:,0])[0]
        xdata = np.arange(sl)/self.getSampling()
        labels = ["{:.2f}".format(xdata[int(i * (sl-1))]*1000) for i in [0, 0.2, 0.4, 0.6, 0.8,1]]

        self.mw.ui.widget_9.canvas.ax.set_xticklabels(labels)
        self.mw.ui.widget_9.canvas.ax.set_ylim([-40, 20])
        self.mw.ui.widget_9.canvas.draw()


    def gettableLength(self):
        return self.tablelength

    def getMapClusters(self):
        return self.origindex


    def handleSelection(self, currentitem):
        if (self.previousbutton != None):
            if self.previousbutton == currentitem: #If I pressed the current cluster, then deactivate it
                self.mw.deselectList1()
            else:
                self.previousbutton = currentitem
        cx = self.centres[0,self.clf][currentitem]
        cy = self.centres[1,self.clf][currentitem]

        self.handleDataExploration(currentitem,cx,cy,self.clf[currentitem])


    def setTimeHistogram(self):
        times = self.getTimes()
        self.nbins,tb = np.histogram(self.times,1000)
        tbs = np.ceil(tb[2]-tb[1])

        self.mint = min(self.times)
        self.maxt = max(self.times)
        mint = self.mint
        maxt = self.maxt
        minb = min(self.nbins)
        maxb = max(self.nbins)

        #Time in seconds
        mints = mint/self.sampling
        maxts = maxt/self.sampling

        stepy =(maxb -minb)/3.0
        yticks = np.arange(int(np.floor(minb)), int(np.ceil(maxb)+1), stepy)

        step = (int(np.ceil(maxt)) - int(np.floor(mint)))/4.0
        xticks = np.arange(int(np.floor(mint)),int(np.ceil(maxt))+1, step)

        #Labels in seconds
        step2 = (int(np.ceil(maxts)) - int(np.floor(mints)))/4.0
        xtickl = np.arange(int(np.floor(mints)),int(np.ceil(maxts))+1,step2)
        strings = ["%.1f" % number for number in xtickl]

        self.mw.setAxisTimeWidget(mint, maxt, minb, maxb+50)
        self.mw.setTicksTimeWidget(xticks = xticks, xlabels= strings, yticks = yticks)
        self.mw.plotTimeWidget(tb[1:]-tbs/2,self.nbins)

    def printClusterSpikes(self, clindex):

        indexc = (self.clusterID == self.clf[clindex]) & (self.sInds)
        dx = self.data[0][indexc]
        dy = self.data[1][indexc]
        abc2 = self.mw.printSpikePoint(dx,dy,s=10, c=tuple(self.colours[self.clf[clindex]]),marker='D',edgecolors='none', alpha=0.5,picker=3 )
        abc2.set_visible(False)
        return abc2

    def getClusterInfo(self):
        indexc = (self.clusterID == self.clf[self.lastClusterSelected]) & (self.sInds)
        maxx = max(self.data[0][indexc])
        minx = min(self.data[0][indexc])
        miny = min(self.data[1][indexc])
        maxy = max(self.data[1][indexc])
        return minx, maxx, miny, maxy

    def printMapAll(self):
        #Set up the data
        self.sortclf = np.array(range(max(self.clf)+1))
        self.tablelength = len(self.clf)
        self.origindex = np.array(range(max(self.clf)+1))
        cx = self.centres[0,self.clf]
        cy = self.centres[1,self.clf]

        c1 = np.array([1, .0, .0, 1])
        c2 = np.array([1, .5, .0, 1])
        c3 = np.array([1, 1, .0, 1])
        c4 = np.array([.5, 1, .0, 1])
        c5 = np.array([0, 1, 1, 1])
        c6 = np.array([0, 0, 1, 1])
        c7 = np.array([0, 0, 0.5, 1])
        c8 = np.array([0, 0, 0, 1])

        listDataColors =[]
        self.listDataRadio = np.empty(len(self.clf), dtype= 'float')
        self.listZOrder= []
        self.listheatcolor = np.empty([len(self.clf),4], dtype= 'float')
        countspikes = np.empty(len(self.clf), dtype= 'float')

        timedif = (self.mw.getTimeEnd()-self.mw.getTimeStart())/self.getSampling()


        datax = self.data[0][self.sInds]
        datay = self.data[1][self.sInds]
        lasti = 0

        most  = self.countspikes.most_common()

        maxc = math.log1p(most[0][1])
        minc = math.log1p(most[len(self.clf)-1][1])



        maxr=100
        minr=10


        for index, c in enumerate(self.clf):
            radio = (math.log1p(self.countspikes[c]) - minc)/maxc
            countspikes[index] = float(self.countspikes[c]/timedif)

            dataradio = maxr*radio
            if radio > 0.9:
                self.listheatcolor[index] = c1

            elif radio > 0.8:
                self.listheatcolor[index] = c2

            elif radio > 0.7:
                self.listheatcolor[index] = c3

            elif radio > 0.6:
                self.listheatcolor[index] = c4

            elif radio > 0.5:
                self.listheatcolor[index] = c5

            elif radio > 0.4:
                self.listheatcolor[index] = c6

            elif radio >= 0.3:
                self.listheatcolor[index] = c7

            elif radio >= 0:
                self.listheatcolor[index] = c8
                dataradio = minr


            self.sortclf[c] = index
            self.origindex[c] = index
            self.listDataRadio[index] = dataradio
            listDataColors.append(self.colours[self.clf[index]])

            self.listZOrder.append(float(1/self.countspikes[c]))


            ind = self.countspikes[c]
            indxs = self.nana[lasti: lasti + ind]
            dx = datax[indxs]
            dy = datay[indxs]

            lasti = ind + lasti




            self.i = 0
            if len(dx) > 2:

                if dx[0] == dx[1]:
                    cir = CirclePolygon((dx[0], dy[0]), .1, facecolor=self.colours[c], alpha=0.2)
                    self.patches.append(cir)
                    continue

                hull = ConvexHull(np.column_stack((dx, dy)))
                ax = dx[hull.vertices]
                ay = dy[hull.vertices]
                polygon = Polygon(np.column_stack((ax, ay)), True, facecolor=self.colours[c],
                                  alpha=0.2)
                self.patches.append(polygon)
            elif len(dx) == 2:
                x = dx[1]- dx[0]
                y = dy[1]- dy[0]
                long = np.absolute(x + y * 1j)
                angle = np.angle([x+y*1j], deg=True)

                elp = Ellipse((x/2,y/2) , width=long, height= .1, angle = angle, facecolor=self.colours[c],
                              alpha=0.2)
                self.patches.append(elp)
            else:
                cir = CirclePolygon((dx[0], dy[0]), .1, facecolor=self.colours[c],
                                    alpha=0.2)
                self.patches.append(cir)

        self.countspikes = countspikes
        self.p = PatchCollection(self.patches, alpha=0.2, picker = None)
        self.mw.ui.mpl.canvas.ax.add_collection(self.p)


        self.abc = self.mw.plotMainWidget(cx,cy,s = 10,c=listDataColors,marker='o',edgecolors='none',alpha=0.5,picker=5)


        self.plotradio = self.mw.plotMainWidget(cx,cy,s = self.listDataRadio,c=self.listheatcolor,marker='o',edgecolors='none',
                                                alpha=0.5, zorder=self.listZOrder, picker =5)

        self.mw.setRadialClustersVisible(self.plotradio, False)
        self.mw.setClusterShapesVisible(self.p, False)

        enum = np.array(range(cx.size))
        self.pointsclust = np.column_stack((enum, cx, cy, self.clf))


        self.mw.setAxisMainWidget(min(self.data[0]),max(self.data[0]),min(self.data[1]),max(self.data[1]))
        self.mw.enableClusterControls(True)
        self.mw.drawMainWidget()

        self.mw.exploratoryMode()

    def enableShapesPicker(self,mode):
        if mode:
            self.p.set_picker(5)
        else:
            self.p.set_picker(None)

    def enableClustersPicker(self,mode):
        if mode:
            self.abc.set_picker(5)
        else:
            self.abc.set_picker(None)

    def setRadialClustersVisible(self, mode):
        self.mw.setRadialClustersVisible(self.plotradio, mode)


    def setClusterShapesVisible(self, mode):
        self.mw.setClusterShapesVisible(self.p, mode)
        for patch in self.listselpatches:
            if self.mw.isClusterView():
                self.mw.setPointVisible(patch,True)
                patch.set_picker(5)





    def sigmoid(self, x):
        return 50*(1 / (1 + math.exp(-0.1*x)))

    def getOriginalLimits(self):
        return self.xmin, self.xmax, self.ymin, self.ymax

    def snap(self, x, y):
        #Return the value in self._points closest to (x, y).
        idx = np.nanargmin(((self._points[:,1:3] - (x,y))**2).sum(axis = -1))
        return self._points[idx]

    def snap2(self, x, y):
        #Return the value in self._pointsclust closest to (x, y).
        idx = np.nanargmin(((self.pointsclust[:,1:3] - (x,y))**2).sum(axis = -1))
        return self.pointsclust[idx]


    def snapsel(self,x,y):
        idx = np.nanargmin(((self.selpoints[:,1:3] - (x,y))**2).sum(axis = -1))
        return self.selpoints[idx]

    def snap3(self, x, y):

        idx = np.nanargmin(((self.pointspca[:,1:3] - (x,y))**2).sum(axis = -1))
        return self.pointspca[idx]

    def handleChange(self, clindex, row):

        clindex = int(clindex)

        if self.mw.getCheckedStatus(row) == 2:      #If insert
            point = [clindex, self.centres[0][self.clf[clindex]], self.centres[1][self.clf[clindex]], self.clf[clindex]]                                        #Chec
            self.listSaved.append(point) #Helps to keep track of the last point selected            #Draw red cross
            #Prints point selected in red
            x = point[1]
            y = point[2]
            pointref = self.mw.plotMainWidget(x ,y ,s=25, c='blue',marker = '+', edgecolors ='red', picker = None)

            if not self.mw.isClusterView():
                self.mw.setPointVisible(pointref, False)


            self.listSelected.append(pointref)                                                                          #Add row selected
            self.listRowSelected = np.append(self.listRowSelected, clindex)


            if clindex != self.lastClusterSelected:                                                                     #Prints spikes if not
                #Prints spikes corresponding to the cluster selected                                                    #already selected
                spikes = self.printClusterSpikes(clindex)
                if not self.mw.isClusterView():
                    self.mw.setPointVisible(spikes, True)
                self.listSpikePlot.append(spikes)

                self.addPointstoMap(clindex)
                self.addSelToList(clindex)

            else:
                self.listselpatches.append(self.patchpoint)
                self.listSpikePlot.append(self.spikespoint)
            self.inserthandleChange2(clindex)

        else:                     #If delete
            index = np.where(self.listRowSelected == clindex)[0]
            self.listSaved.pop(index)
            self.listRowSelected = np.delete(self.listRowSelected, index)
            listsel = self.listSelected.pop(index)
            self.mw.deletePoint(listsel)




            if clindex != self.lastClusterSelected:
                #Prints spikes corresponding to the cluster selected
                listspi = self.listSpikePlot.pop(index)
                self.mw.deleteSpikes(listspi)
                self.delPointsfromMap(clindex)
                self.delSelFromList(clindex, index)
                #self.deleteLastClusterSel()

                if self.lastSpikeSelected is not None:
                    if self.clusterID[self.lastSpikeSelected] == self.clf[clindex]:
                        self.spikepoint = self.mw.deletePoint(self.spikepoint)
                        self.spikewave = self.mw.deleteLine(self.spikewave)
                        self.vspike = self.mw.deleteVerticalSpike(self.vspike)
                        self.spikeproj = self.mw.deletePoint(self.spikeproj)
                        self.lastSpikeSelected = None

            else:
                self.spikespoint = self.listSpikePlot.pop(index)
                self.patchpoint = self.listselpatches.pop(index)

            self.deletehandlechange2(clindex, index)


        if len(self.listSaved) != 0:
            point = self.listSaved[len(self.listSaved)-1]
            text = "x: " + str(point[1]) + " y: " + str(point[2]) + " cl ID: " + str(int(point[3]))
        else:
            text = "No clusters selected"

        if not self.checkAll:   ##To optimise time in a for loop
            self.mw.printInfo(text)
            self.mw.drawMainWidget()
            self.mw.drawWaveFormWidget()
            self.mw.drawPCAWidget()

    def filterExtraPoint(self, clindex):

        indexf = self.mw.getConfWindow().filterIndexes([clindex])
        ind = np.where(self.filtered == clindex)[0]
        if indexf == [clindex]:
            if len(ind) == 0:
                self.filtered = np.append(self.filtered, clindex)
                self.mw.setRowHidden(clindex, False)
        else:
            if len(ind) != 0:

                self.flagged = np.delete(self.filtered, ind)
                self.mw.setRowHidden(clindex, True)





    def setCheckAll(self, mode):
        self.checkAll = mode

    def getLastClusterSel(self):
        return self.lastClusterSelected

    def deleteLastClusterSel(self):
        if self.lastClusterSelected is not None:
            i = self.lastClusterSelected
            c = self.clf[i]
            x = self.centres[0][c]
            y = self.centres[1][c]
            self.handleDataExploration(i,x,y,c)

    def addPointstoMap(self,clindex):
        indexc = (self.clusterID == self.clf[clindex]) & (self.sInds)
        listi = np.where(indexc == True)[0]
        listx = self.data[0][indexc]
        listy = self.data[1][indexc]
        listc = self.clusterID[indexc]
        cs = np.column_stack((listi,listx,listy,listc))
        self._points = np.append(self._points,cs,axis=0)





    def delPointsfromMap(self,clindex):
        indexc = np.where(self._points[:,3] == self.clf[clindex])[0]
        self._points = np.delete(self._points,indexc,0)

    def addPCAPointstoMap(self,clindex):

        if self.fit is not None:
            csInds = (self.clusterID == self.clf[clindex]) & (self.sInds)
            listi = np.where(csInds == True)[0]
            listx = self.fit[0][csInds]
            listy = self.fit[1][csInds]
            listc = self.clusterID[csInds]
            cs = np.column_stack((listi,listx,listy,listc))
            self.pointspca = np.append(self.pointspca,cs,axis=0)

    def delPCAPointsfromMap(self,clindex):
        ind = np.where(self.pointspca[:,3] == self.clf[clindex])[0]
        self.pointspca = np.delete(self.pointspca,ind,0)



    def handleDataExploration(self,i,x,y,c):

        i = int(i)

        self.mw.printLog("")
        if self.lastClusterSelected != i:
            if self.lastClusterSelected is not None:

                rowlast = self.sortclf[self.clf[self.lastClusterSelected]]
                self.mw.handleTableSelection(rowlast, False)

                self.lastcp = self.clusterpoint
                self.clusterpoint = self.mw.deletePoint(self.clusterpoint)
                self.patchesh = self.mw.deleteHistogram(self.patchesh)  #Removes only the histogram. Fully clf
                self.selectedh = self.mw.removeHistogramPoints(self.selectedh)  #Deletes the scatterpoints
                self.selectedspk = self.mw.removeSelectedSpikes(self.selectedspk) # Deletes spikes in timeline of histogram selection
                self.vspikes = self.mw.removeClusterVTLines(self.vspikes) #Deletes spikes in timeline


                if not (self.lastClusterSelected in self.listRowSelected):
                    self.delPointsfromMap(self.lastClusterSelected)
                    self.delSelPointsfromMap(self.lastClusterSelected)
                    self.spikespoint= self.mw.deleteSpikes(self.spikespoint)
                    self.patchpoint = self.mw.deletePoint(self.patchpoint)
                    self.spikepoint = self.mw.deletePoint(self.spikepoint)
                    self.spikewave = self.mw.deleteLine(self.spikewave)
                    self.vspike = self.mw.deleteVerticalSpike(self.vspike)
                    self.spikeproj = self.mw.deletePoint(self.spikeproj)

                    self.wave = self.mw.deleteLine(self.wave)#Mean of the waveforms
                    self.swave = self.mw.deleteLines(self.swave) #Waveforms at the background
                    self.clustproj  = self.mw.deleteProjection(self.clustproj)
                    self.projfilt  = self.mw.deleteProjection(self.projfilt)
                    self.delPCAPointsfromMap(self.lastClusterSelected)

                    self.lastSpikeSelected = None
                    self.mw.disableSpikeControls(True)


            if self.lastSpikeSelected is not None:
                if self.clf[i] == self.clusterID[self.lastSpikeSelected]:
                    self.mw.disableSpikeControls(False)
                else:
                    self.mw.disableSpikeControls(True)

            self.lastClusterSelected = i


            rownow = self.sortclf[self.clf[i]]

            self.mw.handleTableSelection(rownow, True)


            self.clusterpoint = self.mw.printClusterPoint(x, y, 'blue', None) #Prints the cluster

            if not (i in self.listRowSelected):
                self.spikespoint = self.printClusterSpikes(i)
                self.patchpoint = self.mw.ui.mpl.canvas.ax.add_patch(self.patches[int(i)])
                self.patchpoint.set_picker(5)
                if self.mw.isClusterView() == False:
                    self.mw.setPointsVisible(self.spikespoint,True)
                    self.mw.setPointVisible(self.patchpoint, False)
                    self.patchpoint.set_picker(None)

                self.addPointstoMap(i)
                self.addSelPointsToMap(i)

                self.swave, self.wave = self.printWave(i) #Prints the waveform mean and waveforms
                self.clustproj, self.projfilt = self.printPCACluster(i)
                self.addPCAPointstoMap(i)

                self.mw.setPCAClusterVisible(self.clustproj, self.mw.isCompletePCAView())
                self.mw.setPCAClusterVisible(self.projfilt, not self.mw.isCompletePCAView())

            else:
                #print "ya estaba chequeado, solo actualizare como ultimo selected"
                index = np.where(self.listRowSelected == i)[0]
                self.swave = self.listSpikesSelected[index]
                self.wave = self.listWaveSelected[index]
                self.patchpoint = self.listselpatches[index]
                self.spikespoint = self.listSpikePlot[index]
                if self.fit is not None:
                    self.clustproj  = self.listProjecSelected[index]
                    self.projfilt = self.listProjecFiltSelected[index]

            if self.fit is not None:
                oi = self.lastSpikeSelected
                if oi is not None:
                    self.spikeproj = self.mw.deletePoint(self.spikeproj)
                    if self.isClusterPrinted(self.clusterID[oi]):
                        self.spikeproj=  self.mw.printProjectionPoint(self.fit[0,oi],self.fit[1,oi], c=self.mw.getSpikeColor() , marker='o')


            for k in range(len(self.swave)):
                self.swave[k][0].set_zorder(self.zord)


            self.wave[0].set_zorder(self.zord)
            if self.fit is not None:
                self.clustproj[0].set_zorder(self.zord)
                self.projfilt[0].set_zorder(self.zord)
            self.zord = self.zord + 1


            nbins = self.mw.getBinsNumber()
            scale = self.mw.getScale()
            self.n, self.bins, self.patchesh, self.ii, self.timesh= self.printHistogram(c,scale , nbins) #Prints the histogram

            self.vspikes = self.mw.printVerticalSpikes(self.timesh) #Prints the spike events of the cluster


            text = "Cluster x: " + str(x) + " y: " + str(y) + " cl ID: " + str(int(c))
            self.mw.printInfo(text)
            if self.mw.isClusterView() == False:
                self.mw.setPointsVisible(self.clusterpoint,False)

            index = np.where(self.listRowSelected == i)[0] #If the cluster also belongs to the selected ones
                                                           #then make the button change as currently active
            if index.size != 0:
                self.previousbutton = self.mw.setCurrentRowList1(i)

        else:

            self.mw.deselectList1()

            self.mw.disableSpikeControls(True)
            rowlast = self.sortclf[self.clf[self.lastClusterSelected]]
            self.mw.handleTableSelection(rowlast, False)
            self.clusterpoint = self.mw.deletePoint(self.clusterpoint)
            self.patchesh = self.mw.deleteHistogram(self.patchesh)  #Removes only the histogram. Fully clf
            self.selectedh = self.mw.removeHistogramPoints(self.selectedh)  #Deletes the scatterpoints
            self.selectedspk = self.mw.removeSelectedSpikes(self.selectedspk) #Deletes spikes in timeline
            self.vspikes = self.mw.removeClusterVTLines(self.vspikes)


            if (not (i in self.listRowSelected)): # | (self.mw.getItemList1(i).checkState() != 2):
                self.delPCAPointsfromMap(i)
                self.wave = self.mw.deleteLine(self.wave)#Mean of the waveforms
                self.swave = self.mw.deleteLines(self.swave) #Waveforms at the background
                self.clustproj  = self.mw.deleteProjection(self.clustproj)
                self.projfilt  = self.mw.deleteProjection(self.projfilt)

                self.delPointsfromMap(i)
                self.delSelPointsfromMap(i)
                self.spikespoint= self.mw.deleteSpikes(self.spikespoint)
                self.patchpoint = self.mw.deletePoint(self.patchpoint)

            self.lastClusterSelected = None

            if self.lastSpikeSelected is not None:
                cl = self.clusterID[self.lastSpikeSelected]
                if not self.isClusterPrinted(cl):
                    self.spikepoint = self.mw.deletePoint(self.spikepoint)
                    self.spikewave = self.mw.deleteLine(self.spikewave)
                    self.vspike = self.mw.deleteVerticalSpike(self.vspike)
                    self.spikeproj = self.mw.deletePoint(self.spikeproj)
                    self.lastSpikeSelected = None


            text = "No clusters selected for exploration"
            self.mw.printInfo(text)


        if not self.checkAll:

            self.mw.drawWaveFormWidget()

            self.mw.drawHistogramWidget()
            self.mw.drawLineSpikesWidget()

            self.mw.drawPCAWidget()



            #Main widget
            #print time.clock()
            self.mw.ui.mpl.canvas.ax.draw_artist(self.mw.ui.mpl.canvas.ax.patch)
            self.mw.ui.mpl.canvas.ax.draw_artist(self.abc)

            if self.mw.ui.mpl.ntb.isHeatButtonChecked:
                self.mw.ui.mpl.canvas.ax.draw_artist(self.plotradio)
            else:
                self.mw.ui.mpl.canvas.ax.draw_artist(self.abc)

            for i in range(len(self.listSelected)):
                self.mw.ui.mpl.canvas.ax.draw_artist(self.listSelected[i])
                self.mw.ui.mpl.canvas.ax.draw_artist(self.listselpatches[i])
                self.listselpatches[i].set_picker(5)
                if not self.mw.isClusterView:
                    self.mw.ui.mpl.canvas.ax.draw_artist(self.listSpikePlot[i])

            if self.lastClusterSelected is not None:
                self.mw.ui.mpl.canvas.ax.draw_artist(self.clusterpoint)
                self.mw.ui.mpl.canvas.ax.draw_artist(self.spikespoint)
                if (not (self.lastClusterSelected in self.listRowSelected)):
                    self.mw.ui.mpl.canvas.ax.draw_artist(self.patchpoint)
                    self.patchpoint.set_picker(5)

            if self.lastSpikeSelected is not None:
                self.mw.ui.mpl.canvas.ax.draw_artist(self.spikepoint)

            self.mw.ui.mpl.canvas.update()
            self.mw.ui.mpl.canvas.flush_events()



    def handleChangeHistogram(self):
        if self.lastClusterSelected is not None:
            self.patchesh = self.mw.deleteHistogram(self.patchesh)
            c = self.clf[self.lastClusterSelected]
            nbins = self.mw.getBinsNumber()
            scale = self.mw.getScale()
            self.n, self.bins, self.patchesh, self.ii, self.timesh= self.printHistogram(c,scale ,nbins) #Prints the histogram
            self.mw.drawHistogramWidget()

    def updateSpikeColor(self):
        if self.lastSpikeSelected is not None:
            self.spikepoint.set_color(self.mw.getSpikeColor())
            if self.fit is not None:
                self.spikeproj.set_color(self.mw.getSpikeColor())
        if self.selectedh is not None:
            self.selectedh.set_color(self.mw.getSpikeColor())
        self.mw.drawPCAWidget()
        self.mw.drawMainWidget()

    def updateBoxes(self):
        for c in range(len(self.clf)):
            csInds = (self.clusterIDPCA == self.clf[c])
            self.meanpca = np.append(self.meanpca, np.mean(self.fit[0][csInds]))
            self.mw.addPCAClusterInfo(c,self.meanpca[c])




    def printPCAClusters(self):
        self.listProjecSelected = []
        self.updateBoxes()
        self.mw.setPCAValid()

        for i in range(len(self.listRowSelected)):
            proj, projfilt = self.printPCACluster(self.listRowSelected[i])
            self.addPCAPointstoMap(self.listRowSelected[i])
            zord = self.listWaveSelected[i][0].get_zorder()  #Keep the same zorder as in the wave plot
            proj[0].set_zorder(zord)
            self.listProjecSelected.append(proj)
            self.listProjecFiltSelected.append(projfilt)
            self.mw.setPCAClusterVisible(proj,False) #Hide the complete pca cluster proj, only show the filtered ones


        if self.lastClusterSelected is not None:
            index = np.where(self.listRowSelected == self.lastClusterSelected)[0]
            if index.size != 0:
                self.clustproj = self.listProjecSelected[index]
                self.projfilt = self.listProjecFiltSelected[index]
            else:
                self.clustproj, self.projfilt = self.printPCACluster(self.lastClusterSelected)
                self.addPCAPointstoMap(self.lastClusterSelected)
                self.mw.setPCAClusterVisible(self.clustproj,False)
                zord = self.wave[0].get_zorder()
                self.clustproj[0].set_zorder(zord)
                self.projfilt[0].set_zorder(zord)

        ind = self.lastSpikeSelected
        if (ind is not None):
            if self.isClusterPrinted(self.clusterID[ind]):
                self.spikeproj=  self.mw.printProjectionPoint(self.fit[0,ind],self.fit[1,ind],
                                                              c=self.mw.getSpikeColor(),
                                                              marker='o')

        self.mw.exploratoryPCAMode()
        self.mw.drawPCAWidget()



    def initSelectedClusters(self):
        self.mw.setPointVisible(self.abc, False)
        self.mw.drawMainWidget()

    def finishSelectedClusters(self):
        self.mw.setPointVisible(self.abc, True)
        self.mw.drawMainWidget()

    def addSelToList(self, clindex ):
        patch = self.mw.ui.mpl.canvas.ax.add_patch(self.patches[clindex])
        self.listselpatches.append(patch)
        patch.set_picker(5)
        if not self.mw.isClusterView():
            self.mw.setPointVisible(patch,False)
            patch.set_picker(None)
        self.addSelPointsToMap(clindex)


    def delSelFromList(self, clindex , index):
        patch = self.listselpatches.pop(index)
        self.mw.deletePoint(patch)
        self.delSelPointsfromMap(clindex)



    def addSelPointsToMap(self, clindex):
        x = self.centres[0][self.clf[clindex]]
        y = self.centres[1][self.clf[clindex]]
        c = self.clf[clindex]
        cs = np.column_stack((clindex,x,y,c))
        self.selpoints = np.append(self.selpoints,cs ,axis=0)

    def delSelPointsfromMap(self,clindex):
        indexc = np.where(self.selpoints[:,3] == self.clf[clindex])[0]
        self.selpoints = np.delete(self.selpoints,indexc,0)


    def addPointstoMap(self,clindex):
        indexc = (self.clusterID == self.clf[clindex]) & (self.sInds)
        listi = np.where(indexc == True)[0]
        listx = self.data[0][indexc]
        listy = self.data[1][indexc]
        listc = self.clusterID[indexc]
        cs = np.column_stack((listi,listx,listy,listc))
        self._points = np.append(self._points,cs,axis=0)


    def handlePCAExploration(self,oi,x,y,c):
        oi = int(oi)
        self.handleSpikeExploration(oi, self.data[0][oi], self.data[1][oi], self.clusterID[oi])

    def handleSpikeExploration(self,oi,x,y,c):
        oi = int(oi)
        self.mw.printLog("")
        if self.spikepoint is not None:
            self.spikepoint = self.mw.deletePoint(self.spikepoint)
            self.spikewave = self.mw.deleteLine(self.spikewave)
            self.vspike = self.mw.deleteVerticalSpike(self.vspike)
            self.spikeproj = self.mw.deletePoint(self.spikeproj)

        if self.lastSpikeSelected != oi:
            self.lastSpikeSelected = oi


            if self.lastClusterSelected is not None:
                if c == self.clf[self.lastClusterSelected]:
                    self.mw.disableSpikeControls(False)
                else:
                    self.mw.disableSpikeControls(True)

            self.spikepoint = self.mw.printSpikePoint(self.data[0][oi], self.data[1][oi],c=self.mw.getSpikeColor() , marker= 'o')
            if self.mw.isClusterView() == True:
                self.spikepoint.set_visible(False)
            self.spikewave = self.mw.printSpikeWave(self.shapes[:,oi])
            if (self.fit is not None):

                if self.isClusterPrinted(self.clusterID[oi]):
                    self.spikeproj=  self.mw.printProjectionPoint(self.fit[0,oi],self.fit[1,oi], c=self.mw.getSpikeColor(),
                                                                  marker='o')

            time = (self.times[oi]/self.sampling)
            self.vspike = self.mw.printVerticalSpike(time)
            text = "Spike x: " + str(x) + " y: " + str(y) + " cl ID: " + str(int(c))
            self.mw.printInfo(text)

        else:

            self.lastSpikeSelected = None
            self.mw.disableSpikeControls(True)
            text = "No spikes selected for exploration"
            self.mw.printInfo(text)
        self.mw.drawMainWidget()
        self.mw.drawWaveFormWidget()
        self.mw.drawLineSpikesWidget()
        self.mw.drawPCAWidget()

        #Holds on the waveforms checked

    def inserthandleChange2(self, clindex):

        if self.lastClusterSelected == clindex:
            if self.fit is not None:
                proj = self.clustproj
                projfilt = self.projfilt
            x = self.swave
            plot = self.wave

        else:
            if self.fit is not None:
                proj, projfilt = self.printPCACluster(clindex)
                self.addPCAPointstoMap(clindex)

                if self.mw.isCompletePCAView():
                    self.mw.setPCAClusterVisible(proj, True)
                    self.mw.setPCAClusterVisible(projfilt, False)
                else:
                    self.mw.setPCAClusterVisible(proj, False)
                    self.mw.setPCAClusterVisible(projfilt, True)

            x, plot = self.printWave(clindex)
            for k in range(len(x)):
                x[k][0].set_zorder(self.zord)
            if self.fit is not None:
                proj[0].set_zorder(self.zord)
                projfilt[0].set_zorder(self.zord)
            plot[0].set_zorder(self.zord)
            self.zord = self.zord + 1




        #Reference to the waveforms mean
        self.listWaveSelected.append(plot)
        #Reference to the individual waveforms
        self.listSpikesSelected.append(x)
        #Reference to the last projection to hold on


        if self.fit is not None:
            self.listProjecSelected.append(proj)
            self.listProjecFiltSelected.append(projfilt)
            oi = self.lastSpikeSelected
            if oi is not None:
                if self.isClusterPrinted(self.clusterID[oi]):
                    self.spikeproj=  self.mw.printProjectionPoint(self.fit[0,oi],self.fit[1,oi], c=self.mw.getSpikeColor(),
                                                                  marker='o')


    def deletehandlechange2(self, clindex, index):

        if self.fit is not None:
            proj = self.listProjecSelected.pop(index)
            projfilt = self.listProjecFiltSelected.pop(index)
        item =self.listWaveSelected.pop(index)
        itemspike = self.listSpikesSelected.pop(index)


        if self.lastClusterSelected != clindex:
            if self.fit is not None:
                self.mw.deleteProjection(proj)
                self.mw.deleteProjection(projfilt)
                self.delPCAPointsfromMap(clindex)
            self.mw.deleteLine(item)
            self.mw.deleteLines(itemspike)

        if (self.lastSpikeSelected is not None):
            cl = self.clusterID[self.lastSpikeSelected]
            if not self.isClusterPrinted(cl):
                #self.spikepoint = self.mw.deletePoint(self.spikepoint)
                #self.spikewave = self.mw.deleteLine(self.spikewave)
                #self.vspike = self.mw.deleteVerticalSpike(self.vspike)
                self.spikeproj = self.mw.deletePoint(self.spikeproj)





    def getList1Checked(self):
        return self.listRowSelected



    def isClusterPrinted(self,cluster):
        printed = False

        if self.lastClusterSelected is not None:
            if self.clf[self.lastClusterSelected] == cluster:
                printed = True

        if cluster in self.clf[self.listRowSelected]:
            printed = True
        return printed





    def printPCACluster(self, index):
        complete = None
        filtered = None
        if self.fit is not None:
            csInds=   (self.clusterID == self.clf[index])
            csIndsF = (self.clusterID == self.clf[index]) & self.sInds
            complete = self.mw.printClusterProjection(self.fit[0, csInds], self.fit[1, csInds], self.colours[self.clf[index]])
            filtered = self.mw.printClusterProjection(self.fit[0, csIndsF], self.fit[1, csIndsF], self.colours[self.clf[index]],picker = 5)
        return complete, filtered

    def gobackwards(self):

        self.mw.printLog("")

        index = -1

        timeref = self.lastSpikeSelected
        for i in range(len(self.ii)):
            if  timeref == self.ii[i]:
                index = i
                break
        if index ==  -1:
            return
        #Decrement time
        if index == 0:
            index = (len(self.ii))
        index = index -1
        index2= np.where(self._points[:,0] == self.ii[index])
        add = self._points[int(index2[0])]
        self.handleSpikeExploration(add[0], add[1], add[2], add[3])

    def goafterwards(self):
        self.mw.printLog("")

        index = -1

        timeref = self.lastSpikeSelected
        for i in range(len(self.ii)):
            if timeref == self.ii[i]:
                index = i
                break

        if index ==  -1:
            return
        #Increment time
        if index == (len(self.ii)-1):
            index = -1
        index = index +1

        index2= np.where(self._points[:,0] == self.ii[index])
        add = self._points[int(index2[0])]
        self.handleSpikeExploration(add[0], add[1], add[2], add[3])

    def printSelectedBin(self, x):
        index = -1
        indexes = []
        indexes2 = []

        for i in range(len(self.bins)):
            if x >= self.bins[i]:
                index = i
            else:
                break

        startbin = self.bins[index]
        endbin= self.bins[index+1]

        if index == (len(self.bins)-2):
            endbin= endbin + 0.1

        for i in range(len(self.hx)):
                if (self.hx[i] >= startbin) & (self.hx[i] < endbin):
                    indexes.append(i)

        for i in range(len(indexes)):
            indexes2.append(indexes[i])
            indexes2.append(indexes[i] + 1)

        indexes2 = np.unique(indexes2)

        if len(indexes2) != 0:
            dx = self.data[0][self.ii[indexes2]]
            dy = self.data[1][self.ii[indexes2]]
            ref = self.mw.printSpikePoint(dx,dy,s = 25,c=self.mw.getSpikeColor(),marker='x',alpha=0.5,
                                          edgecolors= self.mw.getSpikeColor())
        else:
            ref = None


        return ref, indexes2

    def handleHistogramPick(self,x):
        #print len(self.pointalone)
        if  self.selectedh is not None:
            self.selectedh = self.mw.removeHistogramPoints(self.selectedh)
            self.selectedspk = self.mw.removeSelectedSpikes(self.selectedspk)

        if self.lastxleftbinSelected != x:
            self.lastxleftbinSelected = x
            self.selectedh, self.timesselected = self.printSelectedBin(x)
            if len(self.timesselected) != 0:
                self.selectedspk = self.mw.printVerticalSelectedSpikes(self.timesh[self.timesselected])
            if self.mw.isClusterView():
                self.mw.setPointsVisible(self.selectedh, False)
        else:
            self.lastxleftbinSelected = -1
            text = "No bins selected for exploration"
        self.mw.drawMainWidget()
        self.mw.drawLineSpikesWidget()

    def printHistogram(self,c, type, nbins = 10):
        c = int(c)
        n = 0
        b = 0
        patches = []
        ii=[]
        timesspike =[]

        hind = (self.clusterID == c) & (self.sInds)
        #Positions in the original array
        ii = np.where(hind == True)

        times2 = self.times[hind]
        timesspike = times2/self.sampling
        self.hx = np.diff(timesspike)

        if  type == 'linear':
            bins = nbins
        else:
            if len(self.hx)>=2:
                edges = nbins + 1
                minb = np.min(self.hx)
                maxb = np.max(self.hx)

                if minb == 0.0:
                    minb = 0.000001
                bins = np.logspace(np.log10(minb), np.log10(maxb), edges)
            else:
                bins = nbins

        if len(self.hx)>=2:
             n, bins, patches = self.mw.plotHistogram(self.hx, bins, scaletype = type)
        else:
            self.patchesh = self.mw.deleteHistogram(self.patchesh)



        return n,bins,patches,ii[0], timesspike

    def printWave(self, clindex):
        clindex = int(clindex)
        clShapes = lambda cid: self.shapes[:,((self.clusterID == cid)&(self.sInds))]
        myShapes = clShapes(self.clf[clindex])
        color = self.colours[self.clf[clindex]]
        return self.mw.plotWaveComplete(myShapes, color, self.clf[clindex])

    def completepcaviewcontrol(self):

        for i in range(len(self.listRowSelected)):
            self.mw.setPCAClusterVisible(self.listProjecSelected[i],True)
            self.mw.setPCAClusterVisible(self.listProjecFiltSelected[i], False)

        if self.lastClusterSelected not in self.listRowSelected:
            self.mw.setPCAClusterVisible(self.clustproj,True)
            self.mw.setPCAClusterVisible(self.projfilt, False)

    def filteredpcaviewcontrol(self):

        for i in range(len(self.listRowSelected)):
            self.mw.setPCAClusterVisible(self.listProjecSelected[i],False)
            self.mw.setPCAClusterVisible(self.listProjecFiltSelected[i], True)

        if self.lastClusterSelected not in self.listRowSelected:
            self.mw.setPCAClusterVisible(self.clustproj,False)
            self.mw.setPCAClusterVisible(self.projfilt, True)

    def clusterviewcontrol(self):

        for i in range(len(self.listSelected)):
            self.mw.setPointVisible(self.listSelected[i],True)
            self.mw.setPointVisible(self.listSpikePlot[i], False)

        self.mw.setPointVisible(self.clusterpoint, True)
        self.mw.setPointVisible(self.spikepoint,False)
        self.mw.setPointVisible(self.spikespoint,False)

        if self.lastClusterSelected is not None:
            if (not (self.lastClusterSelected in self.listRowSelected)):
                self.mw.setPointVisible(self.patchpoint, True)
                self.patchpoint.set_picker(5)

        for i in range(len(self.listSelected)):
            self.mw.setPointVisible(self.listselpatches[i],True)
            self.listselpatches[i].set_picker(5)

        if self.mw.ui.mpl.ntb.isSelectButtonChecked():
            self.mw.setPointVisible(self.abc, False)
        else:
            self.mw.setPointVisible(self.abc, True)



        if self.mw.ui.mpl.ntb.isShapeButtonChecked():
            self.mw.setPointVisible(self.p, True)

        if self.mw.ui.mpl.ntb.isHeatButtonChecked():
            self.mw.setPointVisible(self.plotradio, True)

        self.mw.setPointVisible(self.selectedh, False)

    def spikesviewcontrol(self):

        for i in range(len(self.listSelected)):
            self.mw.setPointVisible(self.listSelected[i],False)
            self.mw.setPointVisible(self.listSpikePlot[i], True)

        self.mw.setPointVisible(self.clusterpoint, False)
        self.mw.setPointVisible(self.spikepoint, True)
        self.mw.setPointVisible(self.spikespoint, True)

        if self.lastClusterSelected is not None:
            if (not (self.lastClusterSelected in self.listRowSelected)):
                self.mw.setPointVisible(self.patchpoint, False)
                self.patchpoint.set_picker(None)

        for i in range(len(self.listSelected)):
            self.mw.setPointVisible(self.listselpatches[i],False)
            self.listselpatches[i].set_picker(None)
        self.mw.setPointVisible(self.abc, False)


        #self.mw.setPointVisible(self.p, False)
        if self.mw.ui.mpl.ntb.isShapeButtonChecked():
            self.mw.setPointVisible(self.p, False)

        if self.mw.ui.mpl.ntb.isHeatButtonChecked():
            self.mw.setPointVisible(self.plotradio, False)

        self.mw.setPointVisible(self.selectedh, True)

    def heatviewcontrol(self):

        #for i in range(len(self.listSelected)):
        #    self.mw.setPointVisible(self.listSelected[i],False)


        #self.mw.setPointVisible(self.clusterpoint, False)


        #if self.lastClusterSelected is not None:
        #    if (not (self.lastClusterSelected in self.listRowSelected)):
        #        self.mw.setPointVisible(self.patchpoint, False)
        #        self.patchpoint.set_picker(None)

        #for i in range(len(self.listSelected)):
        #    self.mw.setPointVisible(self.listselpatches[i],False)
        #    self.listselpatches[i].set_picker(None)

        self.mw.setPointVisible(self.abc, False)


        if self.mw.ui.mpl.ntb.isShapeButtonChecked():
            self.mw.setPointVisible(self.p, False)

        self.mw.setPointVisible(self.plotradio, True)





    def checkClusters(self, indexes):

        self.setCheckAll(True)

        for c in indexes:
            row = self.sortclf[self.clf[c]]
            if self.mw.getCheckedStatus(row) == QtCore.Qt.Unchecked:
                self.mw.setCheckedStatus(row, QtCore.Qt.Checked)
                self.handleChange(c, row)

        self.setCheckAll(False)

        self.mw.printLog("Checking clusters finished...\n")
        self.mw.drawMainWidget()
        self.mw.drawPCAWidget()
        self.mw.drawWaveFormWidget()

    def getClusterList(self):
        return self.clf


    def uncheckClusters(self, indexes):

        self.setCheckAll(True)

        for c in indexes:
            row = self.sortclf[self.clf[c]]
            if self.mw.getCheckedStatus(row) == QtCore.Qt.Checked:
                self.mw.setCheckedStatus(row, QtCore.Qt.Unchecked)
                self.handleChange(c, row)

        self.setCheckAll(False)

        self.mw.printLog("Unselecting all clusters finished...\n")

        self.mw.drawMainWidget()
        self.mw.drawWaveFormWidget()
        self.mw.drawPCAWidget()




    def filterClusters(self, inds):

        self.filtered = inds

        listDataColors= []
        listFiltered = [False]*len(self.clf)
        sortclf = np.empty([(max(self.clf)+1),1], dtype = 'int')
        self.tablelength = len(inds)
        listZOrder = []


        for count, i in enumerate(inds):
            listDataColors.append(self.colours[self.clf[i]])
            listZOrder.append(self.listZOrder[i])
            row = self.origindex[self.clf[i]]
            listFiltered[row]= True
            sortclf[self.clf[i]] = count

        self.sortclf = sortclf
        self.mw.proxyModel.setValidRows(listFiltered)
        #self.mw.proxyModel.setTableLength(self.gettableLength())

        #self.mw.ui.tableView.setSortingEnabled(False)
        self.mw.proxyModel.invalidateFilter()
        #self.mw.ui.tableView.setSortingEnabled(True)
        #self.mw.tm.sort(self.mw.ui.tableView.horizontalHeader().sortIndicatorSection())
        self.indicator()




        cx = self.centres[0,self.clf[inds]]
        cy = self.centres[1,self.clf[inds]]

        self.abc = self.mw.deletePoint(self.abc)
        self.abc = self.mw.plotMainWidget(cx,cy,s = 10,c=self.colours[self.clf[inds]],marker='o',edgecolors='none',alpha=0.5,picker=5)

        self.plotradio = self.mw.plotMainWidget(cx,cy,s = self.listDataRadio[inds],c=self.listheatcolor[inds],marker='o',edgecolors='none',
                                                alpha=0.5, zorder=listZOrder, picker=5)

        patches = [self.patches[i] for i in inds]
        self.p = self.mw.deletePoint(self.p)
        self.p = PatchCollection(patches, alpha=0.2, picker = None)
        self.mw.ui.mpl.canvas.ax.add_collection(self.p)
        self.enableShapesPicker(True)

        if not self.mw.ui.mpl.ntb.isShapeButtonChecked():
            self.mw.setClusterShapesVisible(self.p, False)
            self.enableShapesPicker(False)

        if not self.mw.ui.mpl.ntb.isHeatButtonChecked():
            self.mw.setRadialClustersVisible(self.plotradio, False)


        self.pointsclust = np.column_stack((inds, cx, cy, self.clf[inds]))

        self.mw.drawMainWidget()

        self.mw.printLog("Filtering finished...\n")



    def getFilteredRows(self):
        return self.filtered

    def flagClusters(self, indexes):


        for i in indexes:
            if not i in self.flagged:
                ind = self.sortclf[self.clf[i]]
                self.mw.setFlaggedStatus(ind,Qt.Checked)
                self.flagged = np.append(self.flagged, i)



        self.mw.printLog("Flagging finished...\n")


    def getFlaggedList(self):
        return self.clf[self.flagged]


    def unflagClusters(self, indexes):


        for i in indexes:
            if i in self.flagged:
                ind = self.sortclf[self.clf[i]]
                self.mw.setFlaggedStatus(ind,Qt.Unchecked)
                ind = np.where(self.flagged == i)[0]
                self.flagged = np.delete(self.flagged, ind)



        self.mw.printLog("Unflagging finished...\n")


    def getIndexesVisible(self, inds):
        xmin, xmax, ymin, ymax =  self.mw.getZoomLimits()
        return np.where((self.centres[0, self.clf[inds]] < xmax) & (self.centres[0, self.clf[inds]] > xmin) & \
                (self.centres[1, self.clf[inds]] < ymax) & (self.centres[1, self.clf[inds]] > ymin))[0]

    def getCurrentCluster(self):
        if self.lastClusterSelected is not None:
            c =  self.clf[self.lastClusterSelected]
        else:
            c = None
        return c


    def onpickchange(self, event):

        x, y = event.mouseevent.xdata, event.mouseevent.ydata
        if self.mw.ui.mpl.ntb.isSelectButtonChecked():
            i, x, y, c = self.snapsel(x, y)
        else:
            i, x, y, c = self.snap2(x,y)
            self.now = time.clock()

            if (self.now - self.past)< 0.5:
                self.past = self.now
                return
            self.past = self.now
        row = self.sortclf[self.clf[i]]
        self.mw.changeCheckedStatus(row)
        self.handleChange(i,row)

    def onpickexpl(self, event):

        x, y = event.mouseevent.xdata, event.mouseevent.ydata

        if self.mw.ui.mpl.ntb.isSelectButtonChecked():
            i, x, y, c = self.snapsel(x, y)
        else:
            i, x, y, c = self.snap2(x,y)

        c = int(c)

        # self.now = time.clock()
        # print(self.now - self.past)
        # if (self.now - self.past)< 0.5:
        #     self.past = self.now
        #     return
        # self.past = self.now

        row = self.sortclf[c]
        self.mw.scrollTo(row)
        self.handleDataExploration(i,x,y,c)






    def onpickspexpl(self, event):

        x, y = event.mouseevent.xdata, event.mouseevent.ydata
        oi, x, y, c = self.snap(x, y)

        self.handleSpikeExploration(oi,x,y,c)

    def onpicksexpl2(self, event):
        xo, yo = event.xdata, event.ydata
        if xo is None:
            return

        oi, x, y, c = self.snap(xo, yo)

        if norm(np.subtract([xo,yo] , [x,y])) > 0.05:
            if self.lastsp == 1:
                self.lastsp = 0
                self.handleSpikeExploration(oi,x,y,c)
        else:
            if (self.lastsp == 0) | (self.lastspi != oi):

                self.lastsp = 1
                self.lastspi = oi
                self.handleSpikeExploration(oi,x,y,c)





    def onpickpca(self, event):
        x, y = event.mouseevent.xdata, event.mouseevent.ydata

        oi, x, y, c = self.snap3(x, y)
        self.handlePCAExploration(oi,x,y,c)


    def onpickexpl2(self, event):
        xo, yo = event.xdata, event.ydata
        if xo is None:
            return

        oi, x, y, c = self.snap3(xo, yo)

        if norm(np.subtract([xo,yo] , [x,y])) > 0.05:
            if self.last == 1:
                self.last = 0
                self.handlePCAExploration(oi,x,y,c)
        else:
            if (self.last == 0) | (self.lasti != oi):
                self.last = 1
                self.lasti = oi
                self.handlePCAExploration(oi,x,y,c)



    def NData(self):
        #Returns the current number of datapoints.
        return np.shape(self.data)[1]

    def computePCA(self):

        p = PCA(n_components=3,whiten=True)

        # if >1Mio data points, use a subset to avoid memory errors
        if self.nData>1000000:
          # randomly select 1Mio shapes
          inds = np.random.choice(self.nData,1000000,replace=False)
          tf = p.fit(self.shapes[:,inds].T)
          # compute projections
          self.fit = p.transform(self.shapes.T).T
          self.clusterIDPCA = self.clusterID[inds]

        else:
          self.fit = p.fit_transform(self.shapes.T).T
          self.clusterIDPCA = self.clusterID

        self.mw.enablePCAControls(True)
        #plt.axis('scaled')



    def getMaxTimeHistogram(self):
        return np.max(self.nbins)

    def getMaxTime(self):
        return np.max(self.times)

    def getMinTime(self):
        return np.min(self.times)

    def getClusterID(self):
        return self.db.getclusterID()

    def getShapes(self):
        return self.db.getShapes()

    def getData(self):
        return self.db.getData()

    def getCentres(self):
        return self.db.getCentres()

    def getTimes(self):
        return self.db.getTimes()

    def getSampling(self):
        return self.db.getSampling()

    def getClusterNumber(self, index):
        return self.pointsclust[index][3]
