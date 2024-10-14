library(readr)
library(dplyr)
library(tidyr)
library(corrplot)
library(Hmisc)
library(stringr)
library(ggplot2)
library(viridis)

# Load the data
# Load the data
data <- read.csv('/Users/amirostadi/Desktop/1o4 image analysis/zones_expanded_cell_array.csv')

data <- data %>%
  mutate(
    treatmentTime = if_else(antibodyName == "Control", NA, treatmentTime),
    doseLevel = if_else(antibodyName == "Control", NA, doseLevel)
  )


data$identifier <- paste(data$antibodyName,data$imagingType, data$doseLevel, data$treatmentTime,  sep = "_")


data <- data[, c("imagingType", "identifier", "antibodyName", "CorrelationN", "EntropyN", "EnergyN", "HomogeneityN" ,  "CV", "Iso", "MeanN", "StdN")]
colnames(data) <- c("imagingType", "identifier", "antibodyName", "CorrelationN",  "EntropyN",  "EnergyN", "HomogeneityN" , "CV", "Iso", "MeanN", "StdN")
data$CV <- as.numeric(data$CV)
data$Iso <- as.numeric(data$Iso)
# Group the data by identifier
grouped_data <- data %>%
  group_by(identifier)

# Compute the average of each variable within each group while excluding NaN values
averaged_data <- grouped_data %>%
  summarize_at(vars(CorrelationN, EntropyN, EnergyN, HomogeneityN,  CV, Iso, MeanN, StdN), mean, na.rm = TRUE)


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
  pivot_wider(names_from = imagingType, values_from = c( CorrelationN, EntropyN, EnergyN, HomogeneityN,  CV, Iso, MeanN, StdN),  values_fill = NA)

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
    across(starts_with("MeanN"), mean, na.rm = TRUE),
    across(starts_with("StdN"), mean, na.rm = TRUE),
    
    across(starts_with("CV"), mean, na.rm = TRUE),
    across(starts_with("Iso"), mean, na.rm = TRUE)
  ) %>% 
  distinct()  # Remove duplicates after summarizing

# View the summarized data
# summarized_data

summarized_data <- subset(summarized_data, select = -c( dose, time))

# Remove columns with all missing values
summarized_data <- summarized_data[, !colSums(is.na(summarized_data)) == nrow(summarized_data)]

summarized_data <- summarized_data %>%
  select(ends_with("_Dsg3"), ends_with("_E-cad"), ends_with("_IF"), ends_with("F-actin"), ends_with("_RhoA"))

antibody_order <- c("Control", "HLA", "AK23", "PX44", "PX43", "AntiTPO", "AtS13", "PVIgG")  # Replace with your desired order

# Convert the antibody column to a factor with the desired order
summarized_data$antiBody <- factor(summarized_data$antiBody, levels = antibody_order)

# Create the scatter plot with linear regression line
plot <- ggplot(summarized_data) +
  geom_smooth(aes(x = `HomogeneityN_RhoA`, y = `EntropyN_E-cad`), method = "lm", se = TRUE, alpha=.2,  linetype = "dashed") + # Add linear regression line
  geom_point(aes(x = `HomogeneityN_RhoA`, y = `EntropyN_E-cad`, color = antiBody), size=8, alpha=.8) +
  theme_light() +
  labs(
    y = "Entropy E-cad", # 
    x = "Homogeneity RhoA" 
  )+
  theme(
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 28),
    axis.text.x = element_text(angle = 0, size = 26),#, vjust = 0.5, hjust = 1, size = 18), # Rotate x-axis labels to 90 degrees
    axis.text.y = element_text(size = 26), 
    legend.position = "none",  # Remove the legend
    
    # legend.text = element_text(size = 24),
    aspect.ratio = 4/4 # Set aspect ratio to 3:4
  ) +
  scale_color_viridis_d()  # Use the Viridis color palette

plot


















# Create the scatter plot with linear regression line
plot <- ggplot(summarized_data) +
  geom_smooth(aes(x = `EntropyN_RhoA`, y = `EntropyN_Dsg3`), method = "lm", se = TRUE, alpha=.2,  linetype = "dashed") + # Add linear regression line
  geom_point(aes(x = `EntropyN_RhoA`, y = `EntropyN_Dsg3`, color = antiBody), size=8, alpha=.8) +
  theme_light() +
  labs(
    y = "Entropy Dsg3", # Replaced underscore with space
    x = "Entropy RhoA" # Replaced underscore with space
  )+
  theme(
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 28),
    axis.text.x = element_text(angle = 0, size = 26),#, vjust = 0.5, hjust = 1, size = 18), # Rotate x-axis labels to 90 degrees
    axis.text.y = element_text(size = 26), # Rotate x-axis labels to 90 degrees
    legend.position = "none",  # Remove the legend
    
    # legend.text = element_text(size = 24),
    aspect.ratio = 4/4 # Set aspect ratio to 3:4
  ) +
  scale_color_viridis_d()  # Use the Viridis color palette

plot

















# Create the scatter plot with linear regression line
plot <- ggplot(summarized_data) +
  geom_smooth(aes(x = `EnergyN_F-actin`, y = `CV_IF`), method = "lm", se = TRUE, alpha=.2,  linetype = "dashed") + # Add linear regression line
  geom_point(aes(x = `EnergyN_F-actin`, y = `CV_IF`, color = antiBody), size=8, alpha=.8) +
  theme_light() +
  labs(
    y = "CV IF", # Replaced underscore with space
    x = "Energy F-actin" # Replaced underscore with space
  )+
  theme(
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 28),
    axis.text.x = element_text(angle = 0, size = 26),#, vjust = 0.5, hjust = 1, size = 18), # Rotate x-axis labels to 90 degrees
    axis.text.y = element_text(size = 26), # Rotate x-axis labels to 90 degrees
    legend.position = "none",  # Remove the legend
    
    # legend.text = element_text(size = 24),
    aspect.ratio = 4/4 # Set aspect ratio to 3:4
  ) +
  scale_color_viridis_d()  # Use the Viridis color palette

plot