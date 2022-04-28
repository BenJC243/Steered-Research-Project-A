# Steered-Research-Project-A

This project is for partial fulfilment of MSc Bioinformatics at the University of Leicester. Group A consists of: Solomon, Megan, Ben, Rawan and Karthik. All work submitted here is our own.

For the Seurat analysis you will have to download all the files in the data folder and then run the R script I uploaded. To make sure the script runs smoothly you will have to change the "setwd()" command to the specific directory the files you previously downloaded. 

Ensure that these packages are installed before running the seurat analysis script: 
#install.packages("devtools")

#devtools::install_github("VCCRI/CIDR")

#install CIDR from https://github.com/VCCRI/CIDR

# Normalised Counts
FINAL_LOGNORM.csv.gz normalised using same method as Dermanis: "The counts of all genes for any given cell where converted to counts per million (CPM) by diving with the total number of reads and multiplying by 106 followed by conversion to a log base 10" --> see suppl Dermanis file. The counts were however
divided by the count total (column sum) of their respective sample.
FINAL_notLOG: same as above, but wasn't converted to log base 10
Count_matv2: Raw read counts
