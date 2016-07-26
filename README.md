Herding Spikes
==============

Software for high density electrophysiology.

development branch
------------------
Experimental changes will first appear here. Bound to be unstable...

## Sub-projects

 - [onlineDetection](onlineDetection): Online-capable spike detection, done independently on each recording channel.
 - [interpolatingDetection](interpolatingDetection): Spike detection with spatial interpolation. Returns cut-outs for detected events from multiple channels, which allows performing spike localisation.
 - [postProcessing](postProcessing): Programs for removing duplicate events in spikes detected with the online method, and localise spikes detected with the interpolation method.
 - [clustering](clustering): Perform spike sorting by location and PCA on interpolated data.
 - [visualisationtool](visualisationtool): A basic GUI tool for visualising and annotating sorted spikes.

## How to use this software

If you are interested in using this software and have questions or problems, get in touch with [us](http://homepages.inf.ed.ac.uk/mhennig/index.html).

### Fast and simple spike detection

The project *[onlineDetection](onlineDetection)* provides an efficient and reliable algorithm for detecting spikes in single channels.
Main features are a robust noise estimate based on signal percentiles rather than moments, and fast integer based computations allowing real-time performance even when recording 1000s of channels simultaneously. The current implementation only reads 3Brain [Brainwave](http://www.3brain.com/index.php/5/Downloads) files. We are working on a cython based version supporting any file format. While originally developed for high density multielectrode arrays, we expect this to perform well on most extracellular recordings.

### Spike sorting on large scale, high density multielectrode arrays

This process is fully automated, as manual inspection gets very time consuming for large-scale recordings. The following steps lead all the way from raw data to sorted units:

1. Spike detection and spatial signal interpolation - currently runs in 0.1 * real time for 4000 channels, scales linearly with recording duration). Code is in the sub-project [interpolatingDetection](interpolatingDetection).

2. Spatial event localisation, now has about real-time performance (scales linearly with event number).
The relevant code is in the sub-project  [postProcessing](postProcessing).

3. Clustering, code in [clustering](clustering). This step is parallelised and extremely efficient, typical experiments with millions of spikes are clustered in minutes. The code also includes functions for noise removal and quality control.

These methods have for far only been tested with data recorded with the [3Brain Biocam](http://www.3brain.com/biocam-system), we would be keen to other data sets. Steps 1 and 2 would have to be adjusted for other systems, step 3 is already rather generic.

# Contributors
- [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html): Spike sorting
- [Oliver Muthmann](mailto:ollimuh@googlemail.com): Spike detection and localisation
- [Martino Sorbaro](http://martinosorb.github.io): Spike sorting
- [Cesar Juarez Ramirez](mailto:cesaripn2@gmail.com): Visualisation toolkit

# References

J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.
