
#This script produces plots of HVSR data
#Written in January 2024 by Andrea Vergnano (andrea.vergnano@unito.it, cesare.comina@unito.it). You can use, share and modify it as long as you properly refer to the author (Creative Commons BY licence).

#Input: a series of .hv (+.log) and .grid files produced by Geopsy-hv software (with -hv and -rotate options) , and .spec files produced with the -spectrum option.


This script takes as input the results of the three Geopsy processing packages called H/V, H/V Rotate, and H/V Spectrum, respectively stored in files with extensions .hv , .grid  and .spec. Several files are processed in batch mode. The .hv files are processed to create an H/V ratio vs Frequency graph, the .grid files to create a rotation spectrum to assess the polarization of the signal, and the .spec files are graphed as the E, N and Z components of the spectrum. The rotation spectrum is also exported as a .tiff image with corresponding .tfw world file for fast import and visualization in a map in GIS environment (such as QGIS). The .hv files are also processed to determine if the various SESAME criteria are respected, and a spreadsheet report is automatically exported.
