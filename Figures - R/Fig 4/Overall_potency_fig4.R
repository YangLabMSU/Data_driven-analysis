library(ggplot2)

# Read the CSV file
df <- read.csv("/Users/amirostadi/Desktop/1o4 image analysis/potency_df_composit 2.csv")

# Create a new dataframe to store the weighted average values
new_df <- data.frame()

# Calculate the weighted average for each antibodyName
unique_antibodies <- unique(df$antibodyName)

for (antibody in unique_antibodies) {
  subset_df <- df[df$antibodyName == antibody, ]
  weighted_avg <- sum(subset_df$Potency * subset_df$Accuracy) / sum(subset_df$Accuracy)
  
  new_df <- rbind(new_df, data.frame(antibodyName = antibody, Potency = weighted_avg))
}

# Print the new dataframe
# print(new_df)

# Plot the potency scores
plot_potency <- ggplot(new_df, aes(x = reorder(antibodyName, Potency), y = Potency, fill = Potency)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light() +
  scale_fill_gradient(low = "dimgray", high = "dimgray") +
  
  theme(text = element_text(size = 20)) +
  theme(legend.position = "none")+
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1),
    aspect.ratio = 3/1.5
    # plot.margin = margin(2, 2, 2, 4, "cm")
  ) +  
  labs(x = "", y = "Pathogenicity")

plot_potency