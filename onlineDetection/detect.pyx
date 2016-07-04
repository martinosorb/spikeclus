# distutils: language = c++
# distutils: sources = SpkDonline.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
from ctypes import CDLL
import ctypes
from readUtils import openHDF5file, getHDF5params, readHDF5
import time
import multiprocessing
import os

cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long tInc, 
                           long int * Indices, unsigned int nCPU)
        void SetInitialParams(int thres, int maa, int ahpthr, int maxsl, int minsl)
        void openSpikeFile(const char * name)
        # void MedianVoltage(unsigned short * vm)
        void MeanVoltage(unsigned short * vm, int tInc)
        void Iterate(unsigned short * vm, long t0, int tInc)
        void FinishDetection()


def detect(filePath):
 
    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)
   
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    # Number of frames to read in one go (fit to the amount of available memory)
    tInc = min(nFrames, 600000) # Capped at 5 GB of memory approx.

    # Overlap across iterations
    tCut = int(0.001*int(samplingRate)) + int(0.001*int(samplingRate)) + 6

    # Number of processors available
    nCPU = multiprocessing.cpu_count()

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
    MaxSl = int(samplingRate * 1 / 1000 + 0.5) + 1 # compute as in C# program
    MinSl = int(samplingRate * 0.3 / 1000 + 0.5)
    det.SetInitialParams(9, 5, 0, MaxSl, MinSl)
    
    # Open output file
    spikefilename = str.encode(os.path.splitext(filePath)[0] + "_Spikes.txt")
    det.openSpikeFile(spikefilename)

    # Setup timers  
    readT = medianT = iterateT = 0.0; 
    
    # For each chunk of data 
    t0 = 0
    while t0 + tInc <= nFrames:        
        t1 = t0 + tInc

        # Read data
        print '# Reading', t1-t0, 'frames...'
        tic = time.time()
        # Data is indexed by time-step: vm[channel + time * NChannels]
        # must be inverted in order to match the old .brw format.
        vm = readHDF5(rf, t0, t1)
        readT += time.time() - tic

        # Compute median voltage
        print '# Computing median voltage...'
        tic = time.time()
        det.MeanVoltage(&vm[0], tInc) #  det.MedianVoltage(&vm[0]) 
        medianT += time.time() - tic
        
        # Detect spikes
        print '# Detecting spikes...'
        tic = time.time()
        det.Iterate(&vm[0], t0, tInc)
        iterateT += time.time() - tic

        t0 += tInc - tCut # Advance across time
        if t0 < nFrames - tCut: # Last chunk can be smaller
            tInc = min(tInc, nFrames - t0)    

    det.FinishDetection()

    totalT = readT + medianT + iterateT
    print '# Elapsed time:', totalT, 's'
    print '# (Reading:', readT ,'s)'
    print '# (Median:', medianT,'s)'
    print '# (Detecting:', iterateT,'s)'
    print '# (Per frame:', 1000*totalT/nFrames, 'ms)'
