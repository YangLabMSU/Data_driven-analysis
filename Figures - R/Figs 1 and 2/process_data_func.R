analyze_data <- function(data) {
  library(randomForest)
  library(caret)
  library(ggplot2)
  library(gridExtra)
  library(tidyverse)
  library(dplyr)

  # Remove rows with missing values
  data <- na.omit(data)
  
  # Split data into training and testing sets
  set.seed(123)
  trainIndex <- createDataPartition(data$antibodyName, p = 0.8, list = FALSE)
  trainSet <- data[trainIndex, ]
  testSet <- data[-trainIndex, ]
  
  # Train a random forest model
  controlParams <- trainControl(method = "cv", number = 5)
  rfModel <- train(antibodyName ~ ., data = trainSet, method = "rf", 
                   trControl = controlParams, ntree = 200)
  
  # Model importance and predictions
  importance <- varImp(rfModel, scale = TRUE)
  predictions <- predict(rfModel, testSet)
  accuracy <- sum(predictions == testSet$antibodyName) / length(predictions)
  
  # Confusion matrix
  levels <- c("Control", "HLA", "AntiTPO", "PX44", "PX43", "AK23", "AtS13", "PVIgG")
  testSet$antibodyName <- factor(testSet$antibodyName, levels = levels)
  predictions <- factor(predictions, levels = levels)
  cm <- confusionMatrix(data = predictions, reference = testSet$antibodyName)
  cmTable <- as.data.frame(cm$table)
  
  # Calculate percentages for confusion matrix
  total_samples <- tapply(cmTable$Freq, cmTable$Reference, sum)
  cmTable$Percentage <- round(100 * cmTable$Freq / total_samples[cmTable$Reference], 2)
  
  # Create a new dataframe from importance scores
  new_df <- rownames_to_column(importance$importance, var = "Antibody")
  new_df <- new_df[order(new_df$Overall, decreasing = TRUE), ]
  
  # Calculate Potency for each Antibody
  calculate_potency <- function(data, control_data, importance_df) {
    potencies <- data %>%
      group_by(antibodyName) %>%
      summarise(across(everything(), median))
    
    control_medians <- control_data %>%
      summarise(across(everything(), median))
    
    potency_scores <- sapply(2:ncol(potencies), function(i) {
      param <- colnames(potencies)[i]
      weight <- importance_df$Overall[importance_df$Antibody == param]
      abs( weight * (potencies[[param]] - control_medians[[param]]) / control_medians[[param]]  )
    })
    
    potency_df <- potencies %>%
      mutate(Potency = rowSums(potency_scores, na.rm = TRUE))
    
    return(potency_df)
  }
  
  # Filter control data
  control_data <- data %>% filter(antibodyName == "Control")
  
  # Calculate potency scores
  potency_df <- calculate_potency(data, control_data, new_df)
  
  
  library(igraph)
  
  misclassified_edges <- subset(cmTable, Reference != Prediction)
  misclassified_edges <- subset(misclassified_edges, Freq > 2)
  
  misclassified_edges$pair <- apply(misclassified_edges, 1, function(row) {
    paste(sort(c(row['Reference'], row['Prediction'])), collapse = '-')
  })
  
  summarized_edges <- misclassified_edges %>%
    group_by(pair) %>%
    summarise(Total_Freq = sum(Freq))
  
  summarized_edges <- summarized_edges %>%
    separate(pair, into = c("Reference", "Prediction"), sep = "-")
  
  summarized_edges <- summarized_edges %>%
    rename(Freq = Total_Freq)
  
  misclassified_edges <- summarized_edges

  levels <- unique(c('PVIgG', 'AtS13', 'AK23', 'AntiTPO', 'PX44', 'PX43', 'HLA', 'Control'))
  vertices <- data.frame(name=levels)

  graph <- graph_from_data_frame(d = misclassified_edges, directed = FALSE, vertices = vertices)
  E(graph)$weight <- misclassified_edges$Freq
  
  
  
  # Return the plots and accuracy
  return(list(cmTable_data = cmTable, importance_data = importance, accuracy = accuracy, potency_data = potency_df, graph_data = graph, misclassified_edges = misclassified_edges, summarized_edges = summarized_edges))
}


Plot_cond <- function(variable_of_interest, imagingType1) {
  
  # Load required libraries
  library(ggplot2)
  library(RColorBrewer)
  library(tidyverse)
  library(dplyr)
  
  # Load the data
  data <- read.csv('/Users/amirostadi/Desktop/1o4 image analysis/zones_expanded_cell_array_modified.csv')
  
  # Filter data for F-actin imaging type
  data <- subset(data, imagingType == imagingType1)
  
  # Interaction and factor levels
  data$interaction <- interaction(data$doseLevel, data$treatmentTime)
  data$interaction <- factor(data$interaction, levels = c("2ug.4h", "2ug.24h", "10ug.4h", "10ug.24h"))
  order_levels <- c("Control", "PVIgG", "AK23", "AtS13", "HLA", "PX44", "PX43", "AntiTPO")
  data$antibodyName <- factor(data$antibodyName, levels = order_levels)
  
  # Conducting t-tests
  p_values <- list()
  antibody_names <- levels(data$antibodyName)
  antibody_names <- antibody_names[antibody_names != 'Control']
  
  for(ab_name in antibody_names) {
    test_res <- t.test(as.formula(paste(variable_of_interest, "~ antibodyName")),
                       data = data[data$antibodyName %in% c('Control', ab_name),])
    p_values[[ab_name]] <- test_res$p.value
  }
  
  # Prepare annotations
  annotations <- data.frame(antibodyName = names(p_values),
                            MeanN = sapply(names(p_values), function(ab) {
                              max(data[[variable_of_interest]]) 
                            }),
                            label = sapply(p_values, function(p) ifelse(p < 0.05, "*", ""))
  )
  
  annotations$MeanN <- as.numeric(annotations$MeanN)
  
  # Creating the boxplot with annotations for significant differences
  plot <- ggplot(data, aes_string(x = "antibodyName", y = variable_of_interest)) +
    geom_boxplot(size = 1.1, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = .4, aes(color = interaction)) +
    geom_text(data = annotations, aes_string(x = "antibodyName", y = "MeanN", label = "label"),
              nudge_y = 0.0, color = "black", vjust = 0, size = 10) +
    theme_light() +
    theme( 
      text = element_text(size = 24), 
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), 
      aspect.ratio = 1
    ) + 
    labs(y = variable_of_interest, x = "Antibody") +
    coord_cartesian(ylim = c(NA, max(data[[variable_of_interest]]) + 0.1 * max(data[[variable_of_interest]])))
  
  return(plot) 
  
}









