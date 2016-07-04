#import pyximport
#pyximport.install()

from detect import detect

path = '/media/albert/OS/data/'
fileName = '30.brw'

detect(path + fileName)
