# 3Brain 3rd gen .brw (HDF5) 

import h5py

def openHDF5file(path):
    return h5py.File(path, 'r')

def getHDF5params(rf):
    # Read recording variables
    recVars = rf.require_group('3BRecInfo/3BRecVars/')
    bitDepth = recVars['BitDepth'].value[0]
    maxV = recVars['MaxVolt'].value[0]
    minV = recVars['MinVolt'].value[0]
    nFrames = recVars['NRecFrames'].value[0]
    samplingRate = recVars['SamplingRate'].value[0]
    signalInv = recVars['SignalInversion'].value[0]

    # Read chip variables
    chipVars = rf.require_group('3BRecInfo/3BMeaChip/')
    nRows = chipVars['NRows'].value[0]
    nCols = chipVars['NCols'].value[0]
    # nChipCh = nRows * nCols # Total number of channels

    # Get the actual number of channels used in the recording
    nRecCh = len(rf['3BData/Raw'][0])

    # Compute indices
    rawIndices = rf['3BRecInfo/3BMeaStreams/Raw/Chs'].value    

    # Name channels ([0..4095] for fullarray files) 
    chIndices = [(x-1) + (y-1)*nCols for (y,x) in rawIndices] 
    # chIndices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices] # Swap X and Y (old format)

    return (nFrames, samplingRate, nRecCh, chIndices)

def readHDF5(rf, t0, t1):
    ''' In order to use the algorithms designed for the old format, 
    the input data must be inverted.'''
    return 4095 - rf['3BData/Raw'][t0:t1].flatten() 

def readHDF5t(rf, t0, t1):
    ''' Transposed version for the interpolation method. '''
    if t0 <= t1:
        return 4095 - rf['3BData/Raw'][t0:t1].flatten('F')
    else: # Reversed read
        return 4095 - rf['3BData/Raw'][t1:t0].flatten('F') 