# Online spike detection

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com). Cython port by [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html) and [Albert Puente](albertpuente93@gmail.com).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

This branch provides an alpha version which can process raw data the new HDF5 based .brw format (used by 3Brain to store Biocam4096 recordings with Brainwave X). Other formats can be easily integrated by appending their own methods in `readUtils.py` following the same scheme.

## Installation and use

Compile the code by running `python setup.py build_ext --inplace` and execute `python test.py` to detect spikes.

## Implementation

The interpolating detection creates 6 output files, one of which contains parameters used (_info.txt), one for global voltage fluctuations (_avg.txt), two files for spikes (real (_spikes.txt) and virtual (_spikesX.txt) channels) and two files for cut-outs (_shapes.txt,_shapesX.txt). It is not required to set a minimum width and averaged amplitude for spikes here, as the detection uses averages over 3 frames. Additionally, the trigger for recalibrations of the electrodes can be recorded if provided, to remove potential artifacts.

The output is then read and processed by the PYTHON tools provided in [postProcessing](../postProcessing).

## Requirements

### Python packages
* h5py
* cython

### Compilers
* Linux: C++11 compiler (gcc/g++)
* Windows: [Microsoft Visual C++ Compiler for Python 2.7 ](https://www.microsoft.com/en-us/download/details.aspx?id=44266) or any full-fledged Microsoft Visual C++ compiler (>= 2008).
