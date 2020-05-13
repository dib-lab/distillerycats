evaluate_model <- function(optimal_ranger, data, reference_class, set = "train", 
                           study_as_validation = "iHMP", accuracy_csv = NULL,
                           confusion_pdf = NULL) {
  library(ranger)
  library(readr)
  source("scripts/ggplotConfusionMatrix.R")
  
  # calculate prediction accuracy
  pred <- predict(optimal_ranger, data)
  pred_tab <- table(observed = reference_class, predicted = pred$predictions)
  print(pred_tab)
  
  # calculate performance
  if(dim(pred_tab)[1] == 1){
    performance <- pred_tab[2] / sum(pred_tab)
  } else {
    performance <- sum(diag(pred_tab)) / sum(pred_tab)
  }
  print(paste0("ACCURACY = ", performance))
  performance <- data.frame(accuracy = performance, set = set, study = study_as_validation)
  write_csv(performance, accuracy_csv)
  
  # plot pretty confusion matrix
  cm <- caret::confusionMatrix(data = pred$predictions, 
                               reference = factor(reference_class, 
                                                  levels = c("nonIBD", "CD", "UC")))
  plt <- ggplotConfusionMatrix(cm, plot_title = study_as_validation)
  ggsave(filename = confusion_pdf, plot = plt, scale = 1, width = 6, height = 4, dpi = 300)
}

