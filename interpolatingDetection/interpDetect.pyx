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

cdef extern from "SpkDslowFilter.h" namespace "SpkDslowFilter":
    cdef cppclass Detection:
        Detection() except +
        int* SetInitialParams (long nFrames, double nSec, int sf, double sfd, int NCh, int* Indices)
        void openFiles(const char *name)
        void AvgVoltageDefault(short** vm, long t0, int t)
        void InitialEstimation(short** vm, long t0)
        void StartDetection(short** vm, long t0, long nFrames, double nSec, double sfd, int* Indices)
        void skipLastReverse(int skipLast)
        void Iterate(short** vm, long t0)
        void FinishDetection(short** vm, int skipLast)

def detect(rawfilename, sfd, nDumpFrames):
    # Read data and pipe it to the spike detector.
    rf = h5py.File(rawfilename + '.hdf5', 'r')
    cx, cy = rf['x'].value, rf['y'].value
    x1, x2, y1, y2 = np.min(cx), np.max(cx), np.min(cy), np.max(cy)
    nRecCh = len(cx)
    nFrames = len(rf['Ch0'])
    nSec = nFrames / sfd  # the duration in seconds of the recording
    sf = int(sfd)
    nSec = nDumpFrames / sfd

    print '# Number of frames:', nFrames
    print '# Duration (s):', nSec
    print '# Sampling rate: ', sfd
    print '# Number of recorded channels: ', nRecCh
    
    # Messy! To be consistent, X and Y have to be swappped
    cdef np.ndarray[int, mode = "c"] Indices = np.zeros(nRecCh, dtype=long)
    for i in range(nRecCh):
        Indices[i] = (cy[i] - 1) + 64 * (cx[i] - 1)

    cdef Detection* SpkD = new Detection() 
    
    dfTI = SpkD.SetInitialParams(nFrames, nSec, sf, sfd, nRecCh, &Indices[0]);
    SpkD.openFiles(rawfilename);

    sw = datetime.now()

    #cdef np.ndarray[short, ndim = 2, mode = "c"] vm = \
    #    np.zeros((nRecCh,  tInc), dtype=ctypes.c_short)
    tInc = dfTI [2]
    cdef short** vm[nRecCh, tInc] 

    # NChannels = 4096
    # df = 2 # has to be changed in Detection class as well
    if dfTI[0] > 0:
        t1 = 0
        tInc = dfTI [2] # has to be changed in Detection class as well
        # int tCut = (sf / 501 + sf / 1000) / dfTI[0] + 6;
        # int CutOffset = (sf / 1000) / df + 6;
        for t0 in xrange(0, min(200*(dfTI[1]),nFrames-tInc), dfTI[1]):
            ## short[][] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.InitialEstimation(&vm[0,0], t0)

        for t0 in xrange(0, dfTI[1], dfTI[1]):
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.StartDetection (&vm[0,0], t0, nFrames, nSec, sfd, Indices)
            SpkD.Iterate (&vm[0,0], t0)
            t1 += dfTI [1]
        
        for t0 in xrange(dfTI[1], nFrames-tInc, dfTI[1]):
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.Iterate (&vm[0,0], t0)
            t1 += dfTI [1]

        if t1 < nFrames - tInc + dfTI [1] - 1:
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
            SpkD.skipLastReverse (tInc - nFrames + t1)
            SpkD.Iterate (&vm[0,0], t1)

        ## short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
        SpkD.FinishDetection (&vmx[0][0], (int)(tInc - nFrames + t1))
    else:
        t1 = nFrames
        tInc = dfTI [2]
        # int tCut = -(sf / 501 + sf / 1000) / df + 6 +8;
        # int CutOffset = -(sf / 1000) / df + 6;
        for t0 in xrange(nFrames, max(tInc,nFrames-200*dfTI[1]), -dfTI[1]):
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            SpkD.InitialEstimation (&vm[0,0], t0 - tInc)

        for t0 in xrange(nFrames, nFrames-dfTI[1], -dfTI[1]):
            ## cshort[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0-tInc, tInc);
            SpkD.StartDetection (&vm[0,0], t0-tInc, nFrames, nSec, sfd, Indices)
            SpkD.Iterate (&vm[0,0], t0-tInc)
            t1 -= dfTI [1]

        for t0 in xrange (nFrames-dfTI[1], tInc, -dfTI[1]):
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            SpkD.Iterate (&vm[0,0], t0 - tInc)
            t1 -= dfTI [1]

        if t1 > tInc - dfTI [1] + 1:
            SpkD.skipLastReverse ((int)(tInc - t1))
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
            SpkD.Iterate (&vm[0,0], 0)

        ## short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
        SpkD.FinishDetection (&vmx[0][0], (int)(tInc - t1))
    
    time = (datetime.now() - sw).microseconds;
    print 'Elapsed time:', time/1000
    print 'Milliseconds/frame:', time/nFrames
