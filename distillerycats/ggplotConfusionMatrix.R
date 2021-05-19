ggplotConfusionMatrix <- function(m, plot_title = NULL){
  library(caret)
  library(ggplot2)
  library(scales)
  library(tidyr)
  mycaption <- paste("Accuracy", percent_format()(m$overall[1]),
                     "Kappa", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme_minimal() +
    theme(legend.position = "none", 
          text = element_text(size = 20),
          axis.text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5)) +
    labs(caption = mycaption, title = plot_title)
  return(p)
}

