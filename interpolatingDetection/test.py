#import pyximport
#pyximport.install()

from interpDetect import interpDetect
import h5py
import numpy as np
import os

# raw data file
# rawpath = 'data/'
# rawfile = rawpath+'P29_16_05_14_retina02_left_stim2_smallarray_fullfield_raw3'

rawpath = os.environ['HOME'] + '/data/' # Testing in local folder
rawfile = rawpath + 'P29_16_05_14_retina02_left_stim2_smallarray_fullfield_raw3'

print(rawfile)

sampling = 23199.0903585

# run detection
nDumpFrames = int(sampling * 20)  # nFrames;  how many frames to analyze
interpDetect(rawfile, sampling, nDumpFrames)
