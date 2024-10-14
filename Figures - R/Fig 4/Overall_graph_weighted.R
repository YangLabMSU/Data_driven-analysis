library(dplyr)
library(igraph)


# Read all the CSV files into separate data frames
file_paths <- c("/Users/amirostadi/Desktop/1o4 image analysis/For composite graph weighted/misclassified_edges_F-actin.csv",
                "/Users/amirostadi/Desktop/1o4 image analysis/For composite graph weighted/misclassified_edges_RhoA.csv",
                "/Users/amirostadi/Desktop/1o4 image analysis/For composite graph weighted/misclassified_edges_Dsg3.csv",
                "/Users/amirostadi/Desktop/1o4 image analysis/For composite graph weighted/misclassified_edges_Ecad.csv",
                "/Users/amirostadi/Desktop/1o4 image analysis/For composite graph weighted/misclassified_edges_IF.csv")

file_names <- c("F-actin", "RhoA", "Dsg3", "Ecad", "IF") 

dfs <- lapply(file_paths, read.csv)

dfs[[1]]['Freq'] <- dfs[[1]]['Freq']*.552
dfs[[2]]['Freq'] <- dfs[[2]]['Freq']*.670
dfs[[3]]['Freq'] <- dfs[[3]]['Freq']*.736
dfs[[4]]['Freq'] <- dfs[[4]]['Freq']*.793
dfs[[5]]['Freq'] <- dfs[[5]]['Freq']*.706



# Add a column for identifying each file
for (i in 1:length(dfs)) {
  dfs[[i]]$File <- file_names[i]
}

# Merge the data frames on 'Reference' and 'Prediction' columns and sum 'Freq' values
misclassified_edges <- bind_rows(dfs) %>%
  group_by(Reference, Prediction) %>%
  summarise(Freq = sum(Freq))


misclassified_edges$Freq[misclassified_edges$Freq < 15] <- 0

# Ensure that your vertex names match the levels correctly
levels <- unique(c(as.character(misclassified_edges$Reference), as.character(misclassified_edges$Prediction)))
vertices <- data.frame(name=levels)

# Create the graph
graph <- graph_from_data_frame(d = misclassified_edges, directed = FALSE, vertices = vertices)

# Assign weights from your misclassification frequencies
E(graph)$weight <- misclassified_edges$Freq

# layout <- layout_as_tree(graph)
# layout <- layout_with_fr(graph)
# layout <- layout_with_kk(graph)
# layout <- layout_with_mds(graph)
# layout <- layout_as_star(graph)
# layout <- layout_on_grid(graph)
# layout <- layout_with_spring
layout <- layout_nicely(graph)
# layout <- layout_with_graphopt(graph)
# layout <- layout_in_circle(graph)

new_order <- c(1,2,3,4,5,6,7,8)

graph <- permute.vertices(graph, new_order)
scale_factor <- .5

# Plot the graph
plot(graph, layout = layout,
     edge.width = E(graph)$weight * scale_factor,
     vertex.label = V(graph)$name,
     vertex.label.color = "black",
     vertex.size = 40,
     vertex.label.cex = 1.1,
     edge.curved = .3,
     asp = 1,
     frame = FALSE,
     vertex.color = "white"  # This line sets the node color to white
     
)


