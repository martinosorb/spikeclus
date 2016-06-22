# distutils: language = c++
# distutils: sources = SpkDonline.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
from ctypes import CDLL
import ctypes
import h5py
import time
import multiprocessing
import os

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, 
                           long int * Indices, unsigned int nCPU)
        void SetInitialParams(int thres, int maa, int ahpthr, int maxsl, int minsl)
        void openSpikeFile(const char * name)
        void MedianVoltage(unsigned short * vm)
        void MeanVoltage(unsigned short * vm)
        void Iterate(unsigned short * vm, long t0)
        void IterateParallel(unsigned short * vm, long t0)
        void FinishDetection()


def detect(filePath, parallel = True):
 
    # Read data from a .brw (HDF5) file
    rf = h5py.File(filePath, 'r')

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
    chIndices = [(x-1) + (y-1)*nCols for (x,y) in rawIndices] # Swap X and Y
    
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    # Number of frames to read in one go (fit to the amount of available memory)
    tInc = max(nFrames, 10000)

    # Overlap across iterations
    tCut = int(0.001*int(samplingRate)) + int(0.001*int(samplingRate)) + 6

    # Number of processors available
    if parallel:
        nCPU = multiprocessing.cpu_count()
    else:
        nCPU = 1

    # Start detection
    print 'Starting detection...'
    print '# Using', nCPU, 'core(s)'
    print '# Sampling rate:', samplingRate
    print '# Number of recorded channels:', nRecCh
    print '# Analysing frames:', nFrames, 'Seconds:', nSec, 'tCut:', tCut, 'tInc:', tInc
    cdef Detection * det = new Detection()

    # Allocate indices and vm
    cdef np.ndarray[long, mode = "c"] Indices = np.asarray(chIndices, dtype=ctypes.c_long)
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros((nRecCh * tInc), dtype=ctypes.c_ushort)
    
    # Initialise detection algorithm
    det.InitDetection(nFrames, nSec, int(samplingRate), nRecCh, tInc, &Indices[0], nCPU)

    # Set the parameters: int thres, int maa, int ahpthr, int maxsl,  int minsl
    # det.SetInitialParams(9, 5, 0, 8, 3)
    
    # Open output file
    spikefilename = str.encode(os.path.splitext(filePath)[0] + "_Spikes.txt")
    det.openSpikeFile(spikefilename)

    # Setup timers  
    readT = medianT = iterateT = 0.0; 
    
    # For each chunk of data 
    for t0 in range(0, nFrames, tInc - tCut):
        t1 = t0 + tInc 

        if t1 > nFrames: # TODO - The last chunk is potentially smaller. 
            continue

        print 'Iteration', str(1 + t0/(tInc - tCut)) + '/' + str(1 + (nFrames - tInc)/(tInc - tCut)),':'
        
        # Read data
        print '# Reading', t1-t0, 'frames...'
        tic = time.time()
        vm = rf['3BData/Raw'][t0:t1].flatten() # !!! Indexed by time-step: vm[channel + time * NChannels]
        readT += time.time() - tic

        # Compute median voltage
        print '# Computing median voltage...'
        tic = time.time()
        det.MedianVoltage(&vm[0]) # det.MeanVoltage(&vm[0])
        medianT += time.time() - tic
        
        # Detect spikes
        print '# Detecting spikes...'
        tic = time.time()
        if parallel:
            det.IterateParallel(&vm[0], t0)
        else:
            det.Iterate(&vm[0], t0) # Deprecated
        iterateT += time.time() - tic

    det.FinishDetection()

    totalT = readT + medianT + iterateT
    print '# Elapsed time:', totalT, 's'
    print '# (Reading:', readT ,'s)'
    print '# (Median:', medianT,'s)'
    print '# (Detecting:', iterateT,'s)'
    print '# (Per frame:', 1000*totalT/nFrames, 'ms)'
