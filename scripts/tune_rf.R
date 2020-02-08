library(ranger)
set.seed(1)

# inspiration for tuning organization:
# https://uc-r.github.io/random_forests#tune

# read in data ------------------------------------------------------------

ibd_novalidation_filt <- read.csv(snakemake@input[['ibd_novalidation']],
                                  row.names = 1)
diagosis <- read.table(snakemake@input[['ibd_novalidation_diagnosis']], 
                       header = T)
diagnosis <- diagnosis$x
ibd_novalidation_filt$diagnosis <- diagnosis

# make test and train -----------------------------------------------------

train <- sample(nrow(ibd_novalidation_filt), 
                0.7*nrow(ibd_novalidation_filt), 
                replace = FALSE)
train_set <- ibd_novalidation_filt[train, ]
test_set <- ibd_novalidation_filt[-train, ]
diagnosis_train <- diagnosis[train]
diagnosis_test <- diagnosis[-train]

# tune rf -----------------------------------------------------------------

# hyperparameter grid search
hyper_grid2 <- expand.grid(
  mtry       = seq(sqrt(ncol(ibd_novalidation_filt))/2, 
                   sqrt(ncol(ibd_novalidation_filt))*8, 
                   by = 20),
  node_size  = c(5, 7, 10),
  sampe_size = c(.70, .80),
  OOB_RMSE   = 0
)

for(i in 1:nrow(hyper_grid2)) { 
  # train model
  model <- ranger(
    formula         = diagnosis_train ~ ., 
    data            = train_set, 
    num.trees       = 10000,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sampe_size[i],
    seed            = 1
  )
  
  # add OOB error to grid
  hyper_grid2$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

# build optimal rf -------------------------------------------------------

optimal_ranger <- ranger(
  formula         = diagnosis_train ~ ., 
  data            = train_set, 
  num.trees       = 10000,
  mtry            = 960,
  min.node.size   = 5,
  sample.fraction = .7,
  seed            = 1,
  importance      = 'impurity'
)

saveRDS(optimal_ranger, snakemake@output[['optimal_rf']])

pred_test <- predict(optimal_ranger, test_set)
pred_test_tab <-table(observed = diagnosis_test, predicted = pred_test$predictions)
write.table(pred_test_tab, snakemake@output[['pred_test']])

pred_train <- predict(optimal_ranger, train_set)
pred_train_tab <- table(observed = diagnosis_train, predicted = pred_train$predictions)
write.table(pred_train_tab, snakemake@output[['pred_train']])
