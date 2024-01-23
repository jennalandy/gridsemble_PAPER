library(PRROC)
library(pROC)
library(dplyr)

calc_metrics <- function(
    test_statistics,
    fdr, 
    Fdr, 
    truth, 
    true_Fdr, 
    topq, 
    true_fdr = NULL, 
    cutoff = 0.2
) {
  pred = rep(1, length(test_statistics))
  pred[fdr > cutoff] = 0
  
  TP = sum(pred == 1 & truth == 1)
  FP = sum(pred == 1 & truth == 0)
  TN = sum(pred == 0 & truth == 0)
  FN = sum(pred == 0 & truth == 1)
  
  accuracy = (TP + TN)/(TP + TN + FP + FN)
  precision = TP/(TP+FP)
  recall = TP/(TP+FN) # fnr = 1-recall, same as sensitivity
  specificity = TN/(TN + FP) # fpr = 1-specificity
  f1 = (2 * precision * recall) / (precision+recall)
  
  metrics <- list(
    'pr' = get_prauc(fdr, truth),
    'roc' = get_roc(fdr, truth),
    'brier' = get_brier(fdr, truth),
    'Fdrerror' = get_MSE(Fdr, true_Fdr),
    'pr_topq' = get_prauc(fdr[topq], truth[topq]),
    'roc_topq' = get_roc(fdr[topq], truth[topq]),
    'brier_toq' = get_brier(fdr[topq], truth[topq]),
    'Fdrerror_topq' = get_MSE(Fdr[topq], true_Fdr[topq]),
    'pred_pos' = sum(pred == 1),
    'cutoff' = cutoff,
    'accuracy' = accuracy,
    'precision' = precision,
    'recall' = recall,
    'specificity' = specificity,
    'f1' = f1
  )
  
  if (!is.null(true_fdr)) {
    metrics$fdrerror <- get_MSE(fdr, true_fdr)
    metrics$fdrerror_topq <- get_MSE(fdr[topq], true_fdr[topq])
  }
  
  return(metrics)
}


get_prauc <- function(fdr, truth) {
  tryCatch({
    pr <- PRROC::pr.curve(
      scores.class0 = fdr[!as.logical(truth)],
      scores.class1 = fdr[as.logical(truth)]
    )
    return(as.numeric(pr$auc.integral))
  }, error = function(e) {
    return(NA)
  })
}

get_roc <- function(fdr, truth) {
  if (sum(truth) == 0 | sum(!truth) == 0) {
    return(0)
  }
  
  # direction = ">" accounts for inverse
  # relationship between fdr and y
  tryCatch({
    r <- pROC::roc(
      truth ~ fdr,
      direction = ">",
      levels = levels(as.factor(truth))
    )
    return(as.numeric(r$auc))
  }, error = function(e) {
    return(NA)
  })
}

get_brier <- function(fdr, truth) {
  prob_1 = 1-fdr
  mean((prob_1 - as.numeric(truth))**2)
}

get_MSE <- function(estimate, true) {
  mean((true - estimate)^2)
}