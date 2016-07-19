#import pyximport
#pyximport.install()

from interpDetect import interpDetect

rawpath = '/media/albert/DATA/data/'
fileName = 'P29_16_05_14_retina01_right_stim1_smallarray_whitenoise100msHDF5.brw'

interpDetect(rawpath + fileName)