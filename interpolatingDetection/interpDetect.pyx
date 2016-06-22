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
from reader import hdf5reader_3B

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
        void FinishDetection(unsigned short* vm, int skipLast)

def interpDetect(filePath):
    # Read data from a .brw (HDF5) file
    reader = hdf5reader_3B(filePath)
    nFrames, samplingRate, nRecCh, chIndices = reader.getConfig()
    
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    # Allocate indices
    cdef np.ndarray[int, mode = "c"] Indices = np.asarray(chIndices, dtype=ctypes.c_int)
    
    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection() 
    
    dfTI = SpkD.SetInitialParams(nFrames, nSec, int(samplingRate), samplingRate, nRecCh, &Indices[0]);
    SpkD.openFiles(filePath);
    
    tInc = dfTI[2] 

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
            # short[][] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);    
            vm = reader.getData(t0, t0 + tInc).flatten('F') 
            SpkD.InitialEstimation(&vm[0], t0)
        print 'Start detection...'
        for t0 in xrange(0, dfTI[1], dfTI[1]):
            # short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            vm = reader.getData(t0, t0 + tInc).flatten('F')
            SpkD.StartDetection (&vm[0], t0, nFrames, nSec, samplingRate, &Indices[0])
            SpkD.Iterate (&vm[0], t0)
            t1 += dfTI [1]
        for t0 in xrange(dfTI[1], nFrames-tInc, dfTI[1]):
            # short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            vm = reader.getData(t0, t0 + tInc).flatten('F')
            SpkD.Iterate (&vm[0], t0)
            t1 += dfTI [1]
        if t1 < nFrames - tInc + dfTI [1] - 1:
            # short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
            vm = reader.getData(t1, nFrames - t1).flatten('F')
            SpkD.skipLastReverse (tInc - nFrames + t1)
            SpkD.Iterate (&vm[0], t1)

        ## short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
        vmx = reader.getData(t1, nFrames - t1).flatten('F')
        SpkD.FinishDetection (&vmx[0], (int)(tInc - nFrames + t1))
    else:
        t1 = nFrames
        tInc = dfTI [2]
        # int tCut = -(sf / 501 + sf / 1000) / df + 6 +8;
        # int CutOffset = -(sf / 1000) / df + 6;
        print 'Initial estimations...'
        for t0 in xrange(nFrames, max(tInc,nFrames-200*dfTI[1]), -dfTI[1]):
            # short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            vm = reader.getData(t0 - tInc, tInc).flatten('F')
            SpkD.InitialEstimation (&vm[0], t0 - tInc)
        print 'Start detection...'
        for t0 in xrange(nFrames, nFrames-dfTI[1], -dfTI[1]):
            # cshort[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0-tInc, tInc);
            vm = reader.getData(t0 - tInc, tInc).flatten('F')
            SpkD.StartDetection (&vm[0], t0-tInc, nFrames, nSec, samplingRate, &Indices[0])
            SpkD.Iterate (&vm[0], t0-tInc)
            t1 -= dfTI [1]

        for t0 in xrange (nFrames-dfTI[1], tInc, -dfTI[1]):
            # short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            vm = reader.getData(t0 - tInc, tInc).flatten('F')
            SpkD.Iterate (&vm[0], t0 - tInc)
            t1 -= dfTI [1]

        if t1 > tInc - dfTI [1] + 1:
            SpkD.skipLastReverse ((int)(tInc - t1))
            ## short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
            vm = reader.getData(0, t1).flatten('F')
            SpkD.Iterate (&vm[0], 0)
        ## short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
        vmx = reader.getData(0, t1).flatten('F')
        SpkD.FinishDetection (&vmx[0], (int)(tInc - t1))
    
    time = (datetime.now() - sw).microseconds;
    print 'Elapsed time:', time/1000
    print 'Milliseconds/frame:', time/nFrames
