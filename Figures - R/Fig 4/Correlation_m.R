library(readr)
library(dplyr)
library(tidyr)
library(corrplot)
library(Hmisc)
library(stringr)


# Load the data
data <- read.csv('/Users/amirostadi/Desktop/1o4 image analysis/zones_expanded_cell_array.csv')

data <- data %>%
  mutate(
    treatmentTime = if_else(antibodyName == "Control", NA, treatmentTime),
    doseLevel = if_else(antibodyName == "Control", NA, doseLevel)
  )


data$identifier <- paste(data$antibodyName,data$imagingType, data$doseLevel, data$treatmentTime, sep = "_")


data <- data[, c("imagingType", "identifier", "antibodyName", "CorrelationN", "EnergyN", "HomogeneityN" , "EntropyN", "CV", "Iso", "MeanN", "StdN")]
colnames(data) <- c("imagingType", "identifier", "antibodyName",  "CorrelationN", "EnergyN", "HomogeneityN",  "EntropyN", "CV", "Iso", "MeanN", "StdN")

# Group the data by identifier
grouped_data <- data %>%
  group_by(identifier)

# Compute the average of each variable within each group while excluding NaN values
averaged_data <- grouped_data %>%
  summarize_at(vars( CorrelationN, EntropyN, EnergyN, HomogeneityN,  CV, Iso, MeanN, StdN), mean, na.rm = TRUE)


# Add a new column to averaged_data with the second part of the identifier
averaged_data <- averaged_data %>%
  mutate(imagingType = sapply(str_split(identifier, "_"), `[`, 2))

averaged_data <- averaged_data %>%
  mutate(antiBody = sapply(str_split(identifier, "_"), `[`, 1))

averaged_data <- averaged_data %>%
  mutate(dose = sapply(str_split(identifier, "_"), `[`, 3))

averaged_data <- averaged_data %>%
  mutate(time = sapply(str_split(identifier, "_"), `[`, 4))

# Pivot the data to have additional columns for each imaging type
averaged_data <- averaged_data %>%
  pivot_wider(names_from = imagingType, values_from = c(CorrelationN, EntropyN, EnergyN, HomogeneityN, CV, Iso, MeanN, StdN),  values_fill = NA)

averaged_data <- subset(averaged_data, select = -identifier)

# View the modified data
# averaged_data

# Group the data by antiBody, dose, and time
summarized_data <- averaged_data %>% 
  group_by(antiBody, dose, time) %>% 
  mutate(
    across(starts_with("CorrelationN"), mean, na.rm = TRUE),
    across(starts_with("EntropyN"), mean, na.rm = TRUE),
    across(starts_with("EnergyN"), mean, na.rm = TRUE),
    across(starts_with("HomogeneityN"), mean, na.rm = TRUE),
    
    
    across(starts_with("CV"), mean, na.rm = TRUE),
    across(starts_with("Iso"), mean, na.rm = TRUE),
    across(starts_with("MeanN"), mean, na.rm = TRUE),
    across(starts_with("StdN"), mean, na.rm = TRUE)
    
  ) %>% 
  distinct()  


summarized_data <- subset(summarized_data, select = -c(antiBody, dose, time))

# Remove columns with all missing values
summarized_data <- summarized_data[, !colSums(is.na(summarized_data)) == nrow(summarized_data)]

summarized_data <- summarized_data %>%
  select(ends_with("_Dsg3"), ends_with("_E-cad"), ends_with("_IF"), ends_with("_F-actin"), ends_with("_RhoA"))


colnames(summarized_data) <- gsub("N", "", colnames(summarized_data))

rcorr_results <- rcorr(as.matrix(summarized_data), type="pearson")
CorrelationN_matrix <- rcorr_results$r
p_values_matrix <- rcorr_results$P

# Define the significance level, for example 0.05
sig.level <- 0.05

# Plot the CorrelationN matrix with only significant CorrelationNs
corrplot(CorrelationN_matrix, 
         method="circle", 
         tl.col="black", 
         p.mat=p_values_matrix,
         sig.level=sig.level,
         insig="blank"
         # type="lower" 
         )


library(dplyr)
library(stringr)

# Convert the CorrelationN matrix into a dataframe for easier manipulation
cor_df <- as.data.frame(as.table(CorrelationN_matrix))

# Rename the columns for better readability
colnames(cor_df) <- c("Var1", "Var2", "Correlation")

# Filter out the diagonal elements (self-correlations) and duplicate pairs
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]
cor_df <- cor_df[!duplicated(t(apply(cor_df, 1, sort))), ]

# Extract imaging types from variable names
cor_df <- cor_df %>%
  mutate(
    Type1 = sapply(str_split(Var1, "_"), `[`, 2),
    Type2 = sapply(str_split(Var2, "_"), `[`, 2)
  )

# Ensure the correlations come from different imaging types
cor_df_diff_types <- cor_df %>%
  filter(Type1 != Type2)

# Get the top 4 highest absolute correlations from different imaging types
top_4_correlations_diff_types <- cor_df_diff_types %>% 
  arrange(desc(abs(Correlation))) %>% 
  head(10)

# Display the top 4 correlations from different imaging types
print(top_4_correlations_diff_types)
