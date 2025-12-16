Hello!

This folder collects the scripts I wrote about TDEM (Time Domain ElectroMagnetic method), a geophysical method.


Scripts in this folder:

1) Clean_TDEM_data. This scripts:
#a) Asks an input file to be selected interactively
#b) Allows the user to graphycally select which data to keep, to clean noisy parts of the signal, such as those due to superparamagnetic effects.
#c) Calculates apparent resistivity (Rhoapp), and transforms it to effective resistivity and depth, according to Meju (1998). A plot of two Meju's equations for Rho_eff is produced. A plot of Obukhov conditions is produced.
#d) Produces an output .csv file containing also the apparent and effective resistivity, and nicely formatted in columns to be imported in a custom MATLAB script to invert the TDEM data in 1D.




