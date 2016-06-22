import h5py
import numpy as np

class hdf5reader_3B:
    def __init__(self, fileName):
        self.rf = h5py.File(fileName, 'r')
        
    def getConfig(self):
        '''Returns a tuple containing: nFrames, samplingRate, nChannels, indices'''

        # Read recording variables
        recVars = self.rf.require_group('3BRecInfo/3BRecVars/')
        bitDepth = recVars['BitDepth'].value[0]
        maxV = recVars['MaxVolt'].value[0]
        minV = recVars['MinVolt'].value[0]
        nFrames = recVars['NRecFrames'].value[0]
        samplingRate = recVars['SamplingRate'].value[0]
        signalInv = recVars['SignalInversion'].value[0]

        # Read chip variables
        chipVars = self.rf.require_group('3BRecInfo/3BMeaChip/')
        nRows = chipVars['NRows'].value[0]
        nCols = chipVars['NCols'].value[0]
        nChannels = nRows * nCols

        # Compute indices
        rawIndices = self.rf['3BRecInfo/3BMeaStreams/Raw/Chs'].value
        indices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices] # Swap X and Y

        self.nFrames = nFrames

        return (nFrames, samplingRate, nChannels, indices)

    def getData(self, timeStart, timeEnd):
        '''Returns a numpy array indexed by channel.'''
    
        if timeStart <= timeEnd: # Regular read
            return np.transpose(self.rf['3BData/Raw'][timeStart:timeEnd])
        else: # Reverse read
            return np.fliplr(np.transpose(self.rf['3BData/Raw'][timeEnd:timeStart]))