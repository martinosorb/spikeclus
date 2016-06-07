from Cython.Build import cythonize
from setuptools import setup, Extension


setup(
    version='...',
    author='...',
    license='...',
    description='...',
    long_description='...',
    url='...',
    ext_modules=cythonize(Extension(
           "interpDetect",
           sources=["interpDetect.pyx", "interpolatingDetection.cpp"],
           language="c++",
           extra_compile_args=['-std=c++11'],
           )),
    # include_dirs=[numpy.get_include()], requires=['numpy', 'h5py']
)
