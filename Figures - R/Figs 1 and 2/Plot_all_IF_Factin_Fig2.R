source("/Users/amirostadi/Desktop/1o4 image analysis/process_data_func.R")
library(igraph)
library(ggplot2)

graph_IF <- IF_data$graph_data
graph_Factin <- Factin_data$graph_data
potency_df_IF <- IF_data$potency_data
importance_IF <- IF_data$importance_data
importance_Factin <- Factin_data$importance_data
cmTable_IF <- IF_data$cmTable_data
accuracy_IF <- IF_data$accuracy
cmTable_Factin <- Factin_data$cmTable_data
accuracy_Factin <- Factin_data$accuracy
potency_df_Factin <- Factin_data$potency_data

combined_df_Factin_IF_potency <- rbind(
  transform(potency_df_IF, IT = "IF"),
  transform(potency_df_Factin, IT = "Factin")
)


importance_IF$importance$par <- rownames(importance_IF$importance)
importance_IF$importance$IT <- "IF"
importance_Factin$importance$par <- rownames(importance_Factin$importance)
importance_Factin$importance$IT <- "Factin"
combined_df_Factin_IF_importance <- rbind(
  importance_IF$importance,
  importance_Factin$importance
)



new_order <- c(1,2,3,4,5,6,7,8)

graph_reordered <- permute.vertices(graph_Factin, new_order)

# Generate the new layout with the reordered graph
layout <- layout_in_circle(graph_reordered)

scale_factor <- 1

# Plot the graph with the new layout
plot(graph_reordered, layout = layout,
     edge.width = E(graph_reordered)$weight * scale_factor,
     vertex.label = V(graph_reordered)$name,
     vertex.label.color = "darkgreen",
     vertex.size = 65,
     vertex.label.cex = 1.4,
     edge.curved = .2,
     edge.color = 'darkgreen',
     asp = 1.2,
     frame = FALSE,
     vertex.color = "white"  # This line sets the node color to white
     
)


new_order <- c(1,2,3,4,5,6,7,8)

graph_reordered <- permute.vertices(graph_IF, new_order)

# Generate the new layout with the reordered graph
layout <- layout_in_circle(graph_reordered)

scale_factor <- 1

# Plot the graph with the new layout
plot(graph_reordered, layout = layout,
     edge.width = E(graph_reordered)$weight * scale_factor,
     vertex.label = V(graph_reordered)$name,
     vertex.label.color = "firebrick",
     vertex.size = 65,
     vertex.label.cex = 1.4,
     edge.curved = .2,
     edge.color = 'firebrick',
     asp = 1.2,
     frame = FALSE,
     vertex.color = "white"  # This line sets the node color to white
     
)


# Convert antibodyName to a factor with custom order
combined_df_Factin_IF_potency$antibodyName <- factor(combined_df_Factin_IF_potency$antibodyName,
                                                     levels = c("Control", "HLA", "AntiTPO", "PX43", "PX44", "AK23", "PVIgG", "AtS13"))
combined_df_Factin_IF_potency$IT <- factor(combined_df_Factin_IF_potency$IT, 
                                           levels = c("Factin", "IF"))
plot_potency <- ggplot(combined_df_Factin_IF_potency, 
                       aes(x = antibodyName, 
                           y = Potency, 
                           group = IT, 
                           fill = IT)) +    
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + 
  theme_light() + 
  scale_fill_manual(values = c("IF" = "firebrick", "Factin" = "darkgreen")) +   
  theme(text = element_text(size = 20)) + 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1), 
        aspect.ratio = 3/2.3) + 
  labs(x = "", y = "Pathogenicity")


# Convert antibodyName to a factor with custom order
combined_df_Factin_IF_importance$par <- factor(combined_df_Factin_IF_importance$par,
                                               levels = c('CV', 'Iso','Contrast', 'Homogeneity', 'Energy', 'Mean', 'Correlation', 'Std', 'Entropy'))

combined_df_Factin_IF_importance$IT <- factor(combined_df_Factin_IF_importance$IT, 
                                              levels = c("Factin", "IF"))
# Plot with the custom order
plot_importance <- ggplot(combined_df_Factin_IF_importance, 
                          aes(x = par, 
                              y = Overall, 
                              group = IT, 
                              fill = IT)) +    
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + 
  theme_light() + 
  scale_fill_manual(values = c("IF" = "firebrick", "Factin" = "darkgreen")) +   
  theme(text = element_text(size = 20)) + 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1), 
        aspect.ratio = 3/1.9) + 
  labs(x = "", y = "Rel Importance (%)")


plot_cm <- ggplot(data = cmTable_Factin, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Percentage), colour = "white") +
  # geom_text(aes(label = sprintf("%.0f", Percentage)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 24),
        aspect.ratio = 4/4,
        legend.position = "right") +
  labs(fill = "% data")+
  annotate("text", x = 0.5, y = 8.3, label = paste("accuracy:", round(accuracy_Factin*100,1), "%"), 
           size = 8, hjust = 0, vjust = 1)



plot_potency
# plot_importance
# plot_cm

# Plot_cond('Entropy', 'F-actin')













