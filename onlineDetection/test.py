# import pyximport
# pyximport.install()

from detect import detect

path = '/media/albert/DATA/data/'
fileName = 'P29_16_05_14_retina01_right_stim1_smallarray_whitenoise100msHDF5.brw'

# Default call:
detect(path + fileName)

# Or change any parameter using its name:
# detect(path + fileName, Threshold = x, MinAvgAmp = x, AHPthr = x., MaxSl = x, MinSl = x)
