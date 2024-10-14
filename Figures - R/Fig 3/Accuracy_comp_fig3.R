# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the CSV file
data_path <- "/Users/amirostadi/Desktop/1o4 image analysis/combined_data_accuracy.csv"
combined_data_accuracy <- read.csv(data_path)
combined_data_accuracy <- combined_data_accuracy %>%
  filter(Analysis != "+Dose/+Time")

# Convert Analysis to a factor with the specified order
combined_data_accuracy$Analysis <- factor(combined_data_accuracy$Analysis, levels = c("Original", "-10ug/ml", "-2ug/ml")) 

grouped_accuracy <- combined_data_accuracy %>%
  group_by(Analysis, Condition) %>%
  summarise(Accuracy = mean(Accuracy, na.rm = TRUE))

# View the grouped data frame
print(grouped_accuracy)

colors <- c(
  "Original" = "#E69F00",  # Golden Yellow
  "-10ug/ml" = "#56B4E9",  # Sky Blue
  "-2ug/ml" = "#009E73"   # Teal
)
# Create a bar plot of the grouped accuracy data
plot_accuracy <- ggplot(grouped_accuracy, aes(x = Condition, y = Accuracy * 100, fill = Analysis)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_light() +
  coord_cartesian(ylim = c(55, 85)) +  # Zoom in on the y-axis range
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate labels
    aspect.ratio = 1  # Set aspect ratio
  ) +
  scale_fill_manual(values = colors) +   # Specify the colors
  labs(x = "", y = "Accuracy (%)")   # Set axis labels

# Display the plot
print(plot_accuracy)
