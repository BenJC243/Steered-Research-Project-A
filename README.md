# Steered-Research-Project-A

Link to the website - https://bmfsab-solomon.shinyapps.io/SRP-A/?_ga=2.216594454.101771341.1651154579-973842359.1651154579
This project is for partial fulfilment of MSc Bioinformatics at the University of Leicester. Group A consists of: Solomon, Megan, Ben, Rawan and Karthik. All work submitted here is our own.

For the Seurat analysis you will have to download all the files in the 'GROUP A SUBMISSION FOLDER'r and then run the R script uploaded. To make sure the script runs smoothly you will have to change the "setwd()" command to the specific directory the files you previously downloaded. 

Ensure that these packages are installed before running the seurat analysis script: 
#install.packages("devtools")

#devtools::install_github("VCCRI/CIDR")

#install CIDR from https://github.com/VCCRI/CIDR

# Normalised Counts
FINAL_LOGNORM.csv.gz normalised using same method as Dermanis: "The counts of all genes for any given cell where converted to counts per million (CPM) by diving with the total number of reads and multiplying by 106 followed by conversion to a log base 10" --> see suppl Dermanis file. The counts were however
divided by the count total (column sum) of their respective sample; 

FINAL_notLOG: same as above, but wasn't converted to log base 10; 

count_matv2: Raw read counts
