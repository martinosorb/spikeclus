# Online spike detection

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com). Cython port by [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

This branch provides an alpha version which can process raw data from any source. The original implementation is only compatible with the old .brw file format used by 3Brain to store Biocam4096 recordings with Brainwave 2.0.

## Installation and use

### Generic code

Underdevelopment. Start by adapting `detect.pyx` to your needs, then compile the code by running `python setup.py build_ext --inplace`, and edit `test.py` and run this to detect spikes. The output is an ASCII file containing channel number, time stamp and (scaled) amplitude for each spike.

Note: As it is, the code expects spikes to be negative deflections, the sign of the parameter Ascale cen be changed for positive spikes.

### Original Brainwave support

This code is in the folder [brw_files](brw_files).

1. Download the current *BrwExtReader.dll* from [http://www.3brain.com/index.php/5/Downloads](http://www.3brain.com/index.php/5/Downloads) or use the version provided here
2. Compile Main.cs (use any C# compiler for that, e.g. [gmcs](http://www.mono-project.com/docs/about-mono/languages/csharp/) or [Monodevelop](http://www.monodevelop.com/) on Linux), or use the executable provided here.
3. Run (e.g. using [Mono](http://www.mono-project.com/)), a graphical user interface will appear.

## Implementation

The online detection will produce a single text file containing (integers arranged in three columns) channel ids, time stamps (in units of frames) and amplitudes (in units of v/64).

Some parameters can be adjusted by the experimenter. The repolarisation threshold ensures that the measured voltage goes back towards baseline. At least one frame during the maximum width of the spike should be larger than that value. Note, however, that this can be strongly distorted by noise and therefore one should not use a high value here.

The output from the spike detection is then read and processed by the PYTHON tools provided in [postProcessing](../postProcessing).
