source("/Users/amirostadi/Desktop/1o4 image analysis/process_data_func.R")


graph_dsg3 <- dsg3_data$graph_data
graph_ecad <- Ecad_data$graph_data

potency_df_dsg3 <- dsg3_data$potency_data
importance_dsg3 <- dsg3_data$importance_data
importance_ecad <- Ecad_data$importance_data
cmTable_dsg3 <- dsg3_data$cmTable_data
accuracy_dsg3 <- dsg3_data$accuracy
cmTable_ecad <- Ecad_data$cmTable_data
accuracy_ecad <- Ecad_data$accuracy
potency_df_ecad <- Ecad_data$potency_data

combined_df_ecad_dsg3_potency <- rbind(
  transform(potency_df_dsg3, IT = "Dsg3"),
  transform(potency_df_ecad, IT = "Ecad")
)


importance_dsg3$importance$par <- rownames(importance_dsg3$importance)
importance_dsg3$importance$IT <- "Dsg3"
importance_ecad$importance$par <- rownames(importance_ecad$importance)
importance_ecad$importance$IT <- "Ecad"
combined_df_ecad_dsg3_importance <- rbind(
  importance_dsg3$importance,
  importance_ecad$importance
)


layout <- layout_in_circle(graph_ecad)
scale_factor <- 1

# Plot the graph_dsg3
 plot(graph_ecad, layout = layout,
     edge.width = E(graph_ecad)$weight * scale_factor,
     vertex.label = V(graph_ecad)$name,
     vertex.label.color = "darkgreen",
     vertex.size = 70,
     vertex.label.cex = 1.5,
     edge.curved = .2,
     edge.color = 'darkgreen',
     asp = 1.2,
     frame = FALSE,
     vertex.color = "white"  # This line sets the node color to white
     
)
#
 
 layout <- layout_in_circle(graph_dsg3)
 scale_factor <- 1
 
 # Plot the graph_dsg3
 plot(graph_dsg3, layout = layout,
      edge.width = E(graph_dsg3)$weight * scale_factor,
      vertex.label = V(graph_dsg3)$name,
      vertex.label.color = "firebrick",
      vertex.size = 70,
      vertex.label.cex = 1.5,
      edge.curved = .2,
      edge.color = 'firebrick',
      asp = 1.2,
      frame = FALSE,
      vertex.color = "white"  # This line sets the node color to white
      
 )
 

# Convert antibodyName to a factor with custom order
combined_df_ecad_dsg3_potency$antibodyName <- factor(combined_df_ecad_dsg3_potency$antibodyName,
                                                     levels = c("Control", "HLA", "AntiTPO", "PX43", "PX44", "AK23", "PVIgG", "AtS13"))
combined_df_ecad_dsg3_potency$IT <- factor(combined_df_ecad_dsg3_potency$IT, 
                                                     levels = c("Ecad", "Dsg3"))
plot_potency <- ggplot(combined_df_ecad_dsg3_potency, 
                       aes(x = antibodyName, 
                           y = Potency, 
                           group = IT, 
                           fill = IT)) +    
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + 
  theme_light() + 
  scale_fill_manual(values = c("Dsg3" = "firebrick", "Ecad" = "darkgreen")) +   
  theme(text = element_text(size = 20)) + 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1), 
        aspect.ratio = 3/2.3) + 
  labs(x = "", y = "Pathogenicity")


# Convert antibodyName to a factor with custom order
combined_df_ecad_dsg3_importance$par <- factor(combined_df_ecad_dsg3_importance$par,
                                                     levels = c('Contrast', 'Homogeneity', 'Energy', 'Mean', 'Correlation', 'Std', 'Entropy'))

combined_df_ecad_dsg3_importance$IT <- factor(combined_df_ecad_dsg3_importance$IT, 
                                           levels = c("Ecad", "Dsg3"))
# Plot with the custom order
plot_importance <- ggplot(combined_df_ecad_dsg3_importance, 
                       aes(x = par, 
                           y = Overall, 
                           group = IT, 
                           fill = IT)) +    
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + 
  theme_light() + 
  scale_fill_manual(values = c("Dsg3" = "firebrick", "Ecad" = "darkgreen")) +   
  theme(text = element_text(size = 20)) + 
  theme(text = element_text(size = 24), 
        axis.text.x = element_text(angle = 0, hjust = .5, vjust = 1), 
        aspect.ratio = 3/1.9) + 
  labs(x = "", y = "Rel Importance (%)")

 
 plot_cm <- ggplot(data = cmTable_dsg3, aes(x = Reference, y = Prediction)) +
   geom_tile(aes(fill = Percentage), colour = "white") +
   # geom_text(aes(label = sprintf("%.0f", Percentage)), vjust = 1) +
   scale_fill_gradient(low = "white", high = "firebrick") +
   theme_light() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
         text = element_text(size = 24),
         aspect.ratio = 4/4,
         legend.position = "right") +
   labs(fill = "% data")+
   annotate("text", x = 0.5, y = 8.3, label = paste("accuracy:", round(accuracy_dsg3*100,1), "%"), 
            size = 8, hjust = 0, vjust = 1)
 
 
 
  plot_potency
 # plot_importance
  # plot_cm
 
 # Plot_cond('Entropy', 'E-cad')
 
 
 
 
 
 
 
 
 
 
 
 
 
