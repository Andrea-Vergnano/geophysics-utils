# This script reads a raw data file from TEM-FAST instrument, a TDEM (Time-domain electromagnetic method) instrument.
# This raw data file contains some metadata lines, and then 4 columns, containing Channel, Time (in microseconds), Signal (in V/A), and Error (in V/A)
# This script runs in R version 4.3.2 and Rstudio 2023.09.1 Build 494.
# If you want to cite it, please cite the publication of mine (Andrea Vergnano) about TDEM and HVSR measurements in Jangany, Madagascar, 2024.

# What does this script do:
#0) Asks an input file to be selected interactively
#1) Allows the user to graphycally select which data to keep, to clean noisy parts of the signal, such as those due to superparamagnetic effects.
#2) Calculates apparent resistivity (Rhoapp), and transforms it to effective resistivity and depth, according to Meju (1998). A plot of two Meju's equations for Rho_eff is produced. A plot of Obukhov conditions is produced.
#3) Produces an output .csv file containing also the apparent and effective resistivity, and nicely formatted in columns to be imported in a custom MATLAB script to invert the TDEM data in 1D.


#An example is provided, L3_01. Raw data is the file without extension, processed data is the file with .csv extension.
