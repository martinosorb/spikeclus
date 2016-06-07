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
    
    ''' TO ADAPT FROM C#
    
    // measure execution time
    var sw = new Stopwatch ();
    sw.Start ();
    //const int NChannels = 4096;
    //const int df = 2;//has to be changed in Detection class as well
    if (dfTI [0] > 0) {
        long t1 = 0;
        long tInc = dfTI [2];//has to be changed in Detection class as well
        //int tCut = (sf / 501 + sf / 1000) / dfTI[0] + 6;
        //int CutOffset = (sf / 1000) / df + 6;
        for (long t0=0; t0<Math.Min(200*(dfTI[1]),nFrames-tInc); t0+=dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.InitialEstimation (vm, t0);
        }
        for (long t0=0; t0<dfTI[1]; t0+=dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.StartDetection (vm, t0, nFrames, nSec, sfd, Indices);
            SpkD.Iterate (vm, t0);
            t1 += dfTI [1];
        }
        for (long t0=dfTI[1]; t0<nFrames-tInc; t0+=dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
            SpkD.Iterate (vm, t0);
            t1 += dfTI [1];
        }
        if (t1 < nFrames - tInc + dfTI [1] - 1) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
            SpkD.skipLastReverse ((int)(tInc - nFrames + t1));
            SpkD.Iterate (vm, t1);
        }
        short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, t1, nFrames - t1);
        SpkD.FinishDetection (vmx, (int)(tInc - nFrames + t1));
    } else {
        long t1 = nFrames;
        long tInc = dfTI [2];
        //int tCut = -(sf / 501 + sf / 1000) / df + 6 +8;
        //int CutOffset = -(sf / 1000) / df + 6;
        for (long t0=nFrames; t0>Math.Max(tInc,nFrames-200*dfTI[1]); t0-= dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            SpkD.InitialEstimation (vm, t0 - tInc);
        }
        for (long t0=nFrames; t0>nFrames-dfTI[1]; t0-=dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0-tInc, tInc);
            SpkD.StartDetection (vm, t0-tInc, nFrames, nSec, sfd, Indices);
            SpkD.Iterate (vm, t0-tInc);
            t1 -= dfTI [1];
        }
        for (long t0=nFrames-dfTI[1]; t0>tInc; t0-= dfTI[1]) {
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0 - tInc, tInc);
            SpkD.Iterate (vm, t0 - tInc);
            t1 -= dfTI [1];
        }
        if (t1 > tInc - dfTI [1] + 1) {
            SpkD.skipLastReverse ((int)(tInc - t1));
            short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
            SpkD.Iterate (vm, 0);
        }
        short[] [] vmx = brwRdr.GetRawDataADCCounts (Channels, 0, t1);
        SpkD.FinishDetection (vmx, (int)(tInc - t1));
    }
    sw.Stop ();
    Console.WriteLine ("Elapsed time: {0}", sw.Elapsed); // TimeSpan
    Console.WriteLine ("Milliseconds/frame: {0}", sw.Elapsed.TotalMilliseconds / nFrames); // TimeSpan 
    '''
