# distutils: language = c++
# distutils: sources = SpkDslowFilter.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
from ctypes import CDLL
import ctypes
from readUtils import openHDF5file, getHDF5params, readHDF5t
import time
import os

cdef extern from "SpkDslowFilter.h" namespace "SpkDslowFilter":
    cdef cppclass InterpDetection:
        InterpDetection() except +
        int* SetInitialParams (long nFrames, double nSec, int sf, double sfd, int NCh, int* Indices)
        void openFiles(const char *name)
        void AvgVoltageDefault(unsigned short* vm, long t0, int t)
        void InitialEstimation(unsigned short* vm, long t0)
        void StartDetection(unsigned short* vm, long t0, long nFrames, double nSec, double sfd, int* Indices)
        void skipLastReverse(int skipLast)
        void Iterate(unsigned short* vm, long t0, int tInc)
        void FinishDetection(unsigned short* vm, int skipLast, int tInc)

def interpDetect(filePath):
    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)
    
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    print "Number of channels:", nRecCh
    print "Sampling rate:", samplingRate
    print "Duration:", nSec

    # Allocate indices
    cdef np.ndarray[int, mode = "c"] Indices = np.asarray(chIndices, dtype=ctypes.c_int)
    
    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection() 
    
    dfTI = SpkD.SetInitialParams(nFrames, nSec, int(samplingRate), samplingRate, nRecCh, &Indices[0]);
    SpkD.openFiles( str.encode(os.path.splitext(filePath)[0]) );
    
    tInc = dfTI[2]
    
    # Allocate vm and vmx
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)

    # Setup timers  
    initialEstT = startT = loopT = finishT = 0.0; 

    if dfTI[0] > 0:
        t1 = 0

        print 'Initial estimations...'

        tic = time.time()
        for t0 in xrange(0, min(200*(dfTI[1]),nFrames-tInc), dfTI[1]):
            vm = readHDF5t(rf, t0, t0 + tInc) 
            SpkD.InitialEstimation(&vm[0], t0)
        initialEstT += time.time() - tic

        print 'Start detection...'
        tic = time.time()
        vm = readHDF5t(rf, 0, tInc)
        SpkD.StartDetection (&vm[0], 0, nFrames, nSec, samplingRate, &Indices[0])
        SpkD.Iterate (&vm[0], 0, tInc)
        t1 += dfTI[1]        
        startT += time.time() - tic

        tic = time.time()
        for t0 in xrange(dfTI[1], nFrames-tInc, dfTI[1]):
            vm = readHDF5t(rf, t0, t0 + tInc)
            SpkD.Iterate (&vm[0], t0, tInc)
            t1 += dfTI[1]          
        loopT += time.time() - tic

        tic = time.time()
        if t1 < nFrames - tInc + dfTI[1] - 1:
            vm = readHDF5t(rf, t1, nFrames)     
            SpkD.skipLastReverse ((int) (tInc - nFrames + t1))
            SpkD.Iterate (&vm[0], t1, nFrames - t1)

        vm = readHDF5t(rf, t1, nFrames)
        SpkD.FinishDetection (&vm[0], (int)(tInc - nFrames + t1), nFrames - t1)
        finishT += time.time() - tic

    else:
        t1 = nFrames
        
        print 'Initial estimations...'

        for t0 in xrange(nFrames, max(tInc,nFrames-200*dfTI[1]), -dfTI[1]):
            vm = readHDF5t(rf, t0 - tInc, t0)
            SpkD.InitialEstimation (&vm[0], t0 - tInc)

        print 'Start detection...'
        vm = readHDF5t(rf, nFrames - tInc, nFrames)
        SpkD.StartDetection (&vm[0], nFrames-tInc, nFrames, nSec, samplingRate, &Indices[0])
        SpkD.Iterate (&vm[0], nFrames-tInc, tInc)
        t1 -= dfTI[1]

        for t0 in xrange (nFrames-dfTI[1], tInc, -dfTI[1]):
            vm = readHDF5t(rf, t0 - tInc, t0)
            SpkD.Iterate (&vm[0], t0 - tInc, tInc)
            t1 -= dfTI[1]

        if t1 > tInc - dfTI[1] + 1:
            SpkD.skipLastReverse ((int)(tInc - t1))
            vm = readHDF5t(rf, 0, t1)
            SpkD.Iterate (&vm[0], 0, t1)
            
        vm = readHDF5t(rf, 0, t1)
        SpkD.FinishDetection (&vm[0], (int)(tInc - t1), t1)
    
    print '# Initialisation time:', initialEstT, 's'
    print '# Start time:', startT, 's'
    print '# Loop time:', loopT, 's'
    print '# Finish time:', finishT, 's'
