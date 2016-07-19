from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy
import os

if os.name == 'nt':
    print "OS: Windows"
    extraCompArgs = ['/DEBUG, /O2']
else:
    print "OS: Unix \ Others"
    extraCompArgs = ['-std=c++11', '-O3']

long_descr = """`See
<http://github.com/>`_.
"""

setup(
    version='1.0',
    author='Oliver Muthmann, Matthias H Hennig, Albert Puente Encinas',
    license='GPL3',
    description='Efficient spike detection for extracellular recordings.',
    long_description=long_descr,
    url='http://github.com/',
    ext_modules=cythonize(Extension(
           "interpDetect",
           sources=["interpDetect.pyx", "SpkDslowFilter.cpp"],
           language="c++",
           extra_compile_args= extraCompArgs,
           )),
    include_dirs=[numpy.get_include()], requires=['numpy', 'h5py']
)
