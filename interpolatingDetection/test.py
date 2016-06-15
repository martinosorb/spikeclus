#import pyximport
#pyximport.install()

from interpDetect import interpDetect
import h5py
import numpy as np
import os

# raw data file
# rawpath = 'data/'
# rawfile = rawpath+'P29_16_05_14_retina02_left_stim2_smallarray_fullfield_raw3'

rawpath = '/media/albert/OS/data/'
rawfile = rawpath + 'P29_16_05_14_retina02_left_stim3_fullarray_fullfield_raw'

print(rawfile)

sampling = 7022.05819854542
seconds = 1

# run detection
nDumpFrames = int(sampling * seconds)  # nFrames
interpDetect(rawfile, sampling, nDumpFrames)
