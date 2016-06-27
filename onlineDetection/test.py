#import pyximport
#pyximport.install()

from detect import detect

path = '/media/albert/OS/data/'
fileName = '10.brw'

detect(path + fileName)
