# distutils: language = c++
# distutils: sources = SpkDslowFilter.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes
from datetime import datetime
from readUtils import openHDF5file, getHDF5params, readHDF5t
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
        void Iterate(unsigned short* vm, long t0)
        void outNum(int i);
        void FinishDetection(unsigned short* vm, int skipLast)

def interpDetect(filePath):
    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)
    
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    # Allocate indices
    cdef np.ndarray[int, mode = "c"] Indices = np.asarray(chIndices, dtype=ctypes.c_int)
    
    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection() 
    
    dfTI = SpkD.SetInitialParams(nFrames, nSec, int(samplingRate), samplingRate, nRecCh, &Indices[0]);
    SpkD.openFiles( str.encode(os.path.splitext(filePath)[0]) );
    
    tInc = dfTI[2] 
    print 'dfTI', dfTI[0], dfTI[1], dfTI[2]
    # Allocate vm and vmx
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)
    cdef np.ndarray[unsigned short, mode = "c"] vmx = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)

    sw = datetime.now()
    # NChannels = 4096
    # df = 2 # has to be changed in Detection class as well
    if dfTI[0] > 0:

        t1 = 0
        tInc = dfTI [2] # has to be changed in Detection class as well
        # int tCut = (sf / 501 + sf / 1000) / dfTI[0] + 6;
        # int CutOffset = (sf / 1000) / df + 6;
        print 'Initial estimations...'
        for t0 in xrange(0, min(200*(dfTI[1]),nFrames-tInc), dfTI[1]):
            vm = readHDF5t(rf, t0, t0 + tInc) 
            SpkD.InitialEstimation(&vm[0], t0)
        print 'Start detection...'
        for t0 in xrange(0, dfTI[1], dfTI[1]):
            vm = readHDF5t(rf, t0, t0 + tInc)
            SpkD.StartDetection (&vm[0], t0, nFrames, nSec, samplingRate, &Indices[0])
            SpkD.Iterate (&vm[0], t0)
            t1 += dfTI [1]
        for t0 in xrange(dfTI[1], nFrames-tInc, dfTI[1]):
            vm = readHDF5t(rf, t0, t0 + tInc)
            SpkD.Iterate (&vm[0], t0)
            t1 += dfTI [1]
            
        if t1 < nFrames - tInc + dfTI [1] - 1:
            SpkD.outNum(2);
            vm = readHDF5t(rf, t1, nFrames)
            SpkD.skipLastReverse ((int) (tInc - nFrames + t1))
            SpkD.Iterate (&vm[0], t1)
            SpkD.outNum(3);
        vmx = readHDF5t(rf, t1, nFrames)
        SpkD.FinishDetection (&vmx[0], (int)(tInc - nFrames + t1))
        print 'ended'
    else:
        t1 = nFrames
        tInc = dfTI [2]
        # int tCut = -(sf / 501 + sf / 1000) / df + 6 +8;
        # int CutOffset = -(sf / 1000) / df + 6;
        print 'Initial estimations...'
        for t0 in xrange(nFrames, max(tInc,nFrames-200*dfTI[1]), -dfTI[1]):
            vm = readHDF5t(rf, t0 - tInc, t0)
            SpkD.InitialEstimation (&vm[0], t0 - tInc)
        print 'Start detection...'
        for t0 in xrange(nFrames, nFrames-dfTI[1], -dfTI[1]):
            vm = readHDF5t(rf, t0 - tInc, t0)
            SpkD.StartDetection (&vm[0], t0-tInc, nFrames, nSec, samplingRate, &Indices[0])
            SpkD.Iterate (&vm[0], t0-tInc)
            t1 -= dfTI [1]

        for t0 in xrange (nFrames-dfTI[1], tInc, -dfTI[1]):
            vm = readHDF5t(rf, t0 - tInc, t0)
            SpkD.Iterate (&vm[0], t0 - tInc)
            t1 -= dfTI [1]

        if t1 > tInc - dfTI [1] + 1:
            SpkD.skipLastReverse ((int)(tInc - t1))
            vm = readHDF5t(rf, 0, t1)
            SpkD.Iterate (&vm[0], 0)
            
        vmx = readHDF5t(rf, 0, t1)
        SpkD.FinishDetection (&vmx[0], (int)(tInc - t1))
    
    time = (datetime.now() - sw).microseconds;
    print 'Elapsed time:', time/1000
    print 'Milliseconds/frame:', time/nFrames
