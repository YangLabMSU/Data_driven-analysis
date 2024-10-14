library(ggplot2)

# Read the CSV file
df <- read.csv("/Users/amirostadi/Desktop/1o4 image analysis/imp_df_composit 2.csv")

new_df <- data.frame()

# Calculate the weighted average for each Par
unique_Par <- unique(df$Par)

for (antibody in unique_Par) {
  subset_df <- df[df$Par == antibody, ]
  weighted_avg <- sum(subset_df$Imp * subset_df$Accuracy) / sum(subset_df$Accuracy)
  
  new_df <- rbind(new_df, data.frame(Par = antibody, Imp = weighted_avg))
}

# Print the new dataframe
# print(new_df)

# Plot the Imp scores
plot_Imp <- ggplot(new_df, aes(x = reorder(Par, Imp), y = Imp, fill = Imp)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light() +
  scale_fill_gradient(low = "dimgray", high = "dimgray") +
  
  theme(text = element_text(size = 20)) +
  theme(legend.position = "none")+
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1),
    aspect.ratio = 3/1.5)+
    # plot.margin = margin(2, 2, 2, 4, "cm")
  labs(x = "", y = "Rel Importance (%)")

plot_Imp 