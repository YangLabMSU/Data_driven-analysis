# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

data <- read.csv('/Users/amirostadi/Desktop/Fig3/zones_expanded_cell_array.csv')
data$Treatment <- factor(data$Treatment)

for(i in 6:ncol(data)) {
  data[, i] <- as.numeric(as.character(data[, i]))
}


# Data preprocessing 
data <- subset(data, imagingType == 'Dsg3')
data$Treatment <- factor(data$Treatment, levels = c("Control", "CS", "SS", "AK23", "AK23+CS", "AK23+SS")) 


# Creating the boxplot with annotations for significant differences
ggplot(data, aes(x = Treatment, y = Entropy)) +
  geom_boxplot(size = 1.1, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = .4, aes(color = Treatment)) +
  # geom_text(data = annotations, aes(x = Treatment, y = MeanN, label = label), 
  #           nudge_y = 0.0, color = "black", vjust = 0, size = 10) + # Adjust `nudge_y` and `vjust` as necessary
  theme_light() +
  theme(
    text = element_text(size = 33),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none", # Remove legend
    aspect.ratio = 1
  ) +
  labs(x = "")
