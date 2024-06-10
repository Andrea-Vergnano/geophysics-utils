Hi!

This script was written by me, Andrea Vergnano, in June 2024. Please feel free to use it for any purpose, citing my contribution in eventual publications.
You can contact me at andrea.vergnano@posteo.eu for asking any doubt or making contributions to the code.


This script basically creates a very optimized dipole dipole sequence for Electrical Resistivity Tomography survey. The optimization works for multichannel georesistivimeters, such as the Syscal Pro by IRIS instruments, that support up to N potential measurements at a time, given a fixed current dipole.  You can specify what is the number of measurements that your instrument can perform at the same time: in the case of Syscal pro this number is 10.

This script accepts as input the number of channels of your instrument, and how many channels you want to go on for doing roll-along measurements: it also creates an optimized roll-along sequence. 



The script can be run with R (and optionally Rstudio, which is an R interface easier to use). 

The output files are saved in the R working directory, which if not specified I think it is the home directory.



