Interpolating spike detection
=============================

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

Note this implementation is currently only compatible with the old .brw file format used by 3Brain to store Biocam4096 recordings with Brainwave 2.0. A generic version with a python interface will be ready shortly.

## Installation and use

1. Download the current *BrwExtReader.dll* from [http://www.3brain.com/index.php/5/Downloads](http://www.3brain.com/index.php/5/Downloads) or use the version provided here
2. Compile Main.cs (use any C# compiler for that, e.g. Monodevelop on Linux), or use the executable provided
3. Run (e.g. using Mono)

## Implementation

The interpolating detection creates 6 output files, one of which contains parameters used (_info.txt), one for global voltage fluctuations (_avg.txt), two files for spikes (real (_spikes.txt) and virtual (_spikesX.txt) channels) and two files for cut-outs (_shapes.txt,_shapesX.txt). It is not required to set a minimum width and averaged amplitude for spikes here, as the detection uses averages over 3 frames. Additionally, the trigger for recalibrations of the electrodes can be recorded if provided, to remove potential artifacts.

The output is then read and processed by the PYTHON tools provided in [postProcessing](../postProcessing).
