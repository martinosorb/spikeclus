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
    nRecCh = nRows * nCols

    # Compute indices
    rawIndices = rf['3BRecInfo/3BMeaStreams/Raw/Chs'].value
    # chIndices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices] # Swap X and Y (no longer required)
    chIndices = [(x-1) + (y-1)*nCols for (y,x) in rawIndices] # Name channels ([0..4095] for fullarray files) 

    return (nFrames, samplingRate, nRecCh, chIndices)

def readHDF5(rf, t0, t1):
    return 4095 - rf['3BData/Raw'][t0:t1].flatten() 