# distutils: language = c++
# distutils: sources = interpolatingDetection.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes
from datetime import datetime

cdef extern from "iterpolatingDetection.h" namespace "SpkDslowFilter":
    cdef cppclass Detection:
        Detection() except +
        int* SetInitialParams (long nFrames, double nSec, int sf, double sfd, int NCh, int* Indices)
        void openSpikeFile(const std::string& name)
        void openShapeFile(const std::string& name)
        void openSpikeXFile(const std::string& name)
        void openShapeXFile(const std::string& name)
        void openInfoFile(const std::string& name)
        void openMeanFile(const std::string& name)
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

    print '# Number of frames:', nframes
    print '# Duration (s):', nSec
    print '# Sampling rate: ', sfd
    print '# Number of recorded channels: ', nRecCh
    
    # Messy! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=long)
    for i in range(nRecCh):
        Indices[i] = (cy[i] - 1) + 64 * (cx[i] - 1)

    cdef Detection* SpkD = new Detection() 
    
    cdef np.ndarray[int, mode = "c"] dfTI = SpkD.SetInitialParams (nFrames, nSec, sf, sfd, nRecCh, Indices);
    
    SpkD.openFiles(rawfilename);
    
    ############################
    
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, & Indices[0])
    det.SetInitialParams(9, 5, 0, 8, 3)

    # open output file
    spikefilename = str.encode(rawfilename + "_Spikes.txt")
    det.openSpikeFile(spikefilename)

    cdef np.ndarray[unsigned short, ndim = 1, mode = "c"] vm = np.zeros((nRecCh * tInc), dtype=ctypes.c_ushort)
    startTime = datetime.now()
    for t0 in range(0, nDumpFrames - tInc, tInc - tCut):
        if (t0 / tInc) % 100 == 0:
            print(str(t0 / sf) + " sec")
        for c in range(nRecCh):
            vm[c * tInc:c * tInc + tInc] = rf['Ch' +
                                              str(c)][t0:t0 + tInc].astype(dtype=ctypes.c_ushort)
        det.MedianVoltage(&vm[0])
        #det.MeanVoltage( & vm[0])  # a bit faster (maybe)
        det.Iterate( & vm[0], t0)
    det.FinishDetection()
    endTime = datetime.now()
    print('Time taken for detection: ' + str(endTime - startTime))
    print('Time per frame: ' + str(1000 * (endTime - startTime) / (nDumpFrames)))
    print('Time per sample: ' + str(1000 *
                                    (endTime - startTime) / (nRecCh * nDumpFrames)))
