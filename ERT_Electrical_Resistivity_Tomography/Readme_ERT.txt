Hello!

This folder collects the scripts I wrote about ERT (Electrical Resistivity Tomography)
I mainly work with R and RStudio. You can install R from this link https://cran.r-project.org/

Consider that I mainly work with Syscal Pro georesistivimeter manufactured by IRIS instruments, therefore the data format, the logic used for the scripts, is made to work with such instruments. It should be simple to adapt these scripts to the file formats used by other instruments.


Scripts in this folder:

1) CleanERT. A script that allows you to preprocess ERT data before inversion.

2) SequencERT. A collection of scripts to create 2D and 3D ERT quadrupole sequences

3) ReordERT. A script that takes any quadrupole sequence and reorders it in a way that electrodes that were just used as current electrodes are not used as potential electrodes, to avoid polarization errors. Also, it optimizes the sequence to be able to remove the first (and then also the second) cable after a certain measurement, while the acquisition is still running, to save time. The folder contains also a preliminary attempt in Fortran.


Each script has example input and output data.
