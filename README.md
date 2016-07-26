Branch with obsolete/old code
=============================

# Herding Spikes
Software for high density multi-electrode array recordings.

## Sub-projects

 - [onlineDetection](onlineDetection): Online-capable spike detection, done independently on each recording channel.
 - [interpolatingDetection](interpolatingDetection): Spike detection with spatial interpolation. Returns cut-outs for detected events from multiple channels, which allows performing spike localisation.
 - [InterpDetectionDev](InterpDetectionDev): Development version of the interpolation and spike localisation code.
 - [postProcessing](postProcessing): Programs for removing duplicate events in spikes detected with the online methot, and localise spikes detected with the interpolation method.
 - [clustering](clustering): Perform spike sorting by location and PCA on interpolated data.
 - [visualisationtool](visualisationtool): A GUI tool for visualising and annotating sorted spikes.


## Contributors
- [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html): Spike sorting
- [Oliver Muthmann](mailto:ollimuh@googlemail.com): Spike detection
- [Martino Sorbaro](http://martinosorb.github.io): Spike sorting
- [Cesar Juarez Ramirez](mailto:cesaripn2@gmail.com): Visualisation toolkit
