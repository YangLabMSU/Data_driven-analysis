
library(dplyr)

# Load the script containing the analyze_data function
source("/Users/amirostadi/Desktop/1o4 image analysis/process_data_func.R")
data <- read.csv('/Users/amirostadi/Desktop/1o4 image analysis/zones_expanded_cell_array_modified.csv')
###

# Filter based on the input imaging_type
data_sel_dsg3 <- subset(data, imagingType == 'Dsg3')
data_sel_dsg3 <- data_sel_dsg3[, c("antibodyName","Correlation", "Entropy", "Contrast", "Homogeneity", "Energy", "Mean", "Std")]
# 
data_sel_Ecad <- subset(data, imagingType == 'E-cad')
data_sel_Ecad <- data_sel_Ecad[, c("antibodyName","Correlation", "Entropy", "Contrast", "Homogeneity", "Energy", "Mean", "Std")]
# 
data_sel_RhoA <- subset(data, imagingType == 'RhoA')
data_sel_RhoA <- data_sel_RhoA[, c("antibodyName","Correlation", "Entropy", "Contrast", "Homogeneity", "Energy", "Mean", "Std")]

data_sel_Factin <- subset(data, imagingType == 'F-actin')
data_sel_Factin <-  data_sel_Factin[, c("antibodyName","Correlation", "Entropy", "Contrast", "Homogeneity", "Energy", "CV", "Iso", "Mean", "Std")]

data_sel_IF <- subset(data, imagingType == 'IF')
data_sel_IF <-  data_sel_IF[, c("antibodyName","Correlation", "Entropy", "Contrast", "Homogeneity", "Energy", "CV", "Iso", "Mean", "Std")]



dsg3_data <- analyze_data(data_sel_dsg3)
Ecad_data <- analyze_data(data_sel_Ecad)
RhoA_data <- analyze_data(data_sel_RhoA)
Factin_data <- analyze_data(data_sel_Factin)
IF_data <- analyze_data(data_sel_IF)



