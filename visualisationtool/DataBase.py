import numpy as np
import h5py.h5ac
#import scipy.special._ufuncs_cxx


class DataBase():

    def setupDatabase(self, spikeFile):
        # RAW: 'Sampling', 'Locations', 'Shapes', 'Times'
        # CLUSTERED: 'Sampling', 'data', 'shapes', 'times', 'expinds', 'cluster_id', 'centres'
        f=h5py.File(spikeFile,'r')
        if f.keys() == ['Sampling','centres','cluster_id','data','expinds','shapes','times']:
            self.clusterID = np.array(f['cluster_id'].value,dtype=int)
            self.shapes = np.array(f['shapes'].value,dtype=float)
            self.data = np.array(f['data'].value,dtype=float)
            self.centres = np.array(f['centres'].value,dtype=float)
            self.times = np.array(f['times'].value,dtype=int)
            self.sampling = f['Sampling'].value
            return "clustered"

        if 'Sampling' and 'Locations' and 'Shapes' and 'Times' in f.keys():
            self.shapes = np.array(f['Shapes'].value,dtype=float).T
            self.data = np.array(f['Locations'].value,dtype=float).T
            self.sampling = f['Sampling'].value
            self.times = np.array(f['Times'].value,dtype=int)
            self.clusterID = None
            self.centres = None
            return "unclustered"

        else:
            return "false"
        f.close()

    def getClusterID(self):
        return self.clusterID

    def getShapes(self):
        return self.shapes

    def getData(self):
        return self.data

    def getCentres(self):
        return self.centres

    def getTimes(self):
        return self.times

    def getSampling(self):
        return self.sampling
