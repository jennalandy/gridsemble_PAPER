get_empirical_Fdr <- function(test_statistics, truth)  {
  Fdr <- rep(NA, length(test_statistics))
  # Pr(null | T <= t) = Pr(truth = 0 | T <= t)
  for (i in 1:length(test_statistics)) {
    t = test_statistics[i]
    Fdr[i] <- mean(1 - truth[abs(test_statistics) >= abs(t)])
  }
  return(Fdr)
}

Fdr_from_fdr <- function(fdr, test_statistics) {
  Fdr = vector(length = length(test_statistics))
  for (i in 1:length(test_statistics)) {
    Fdr[i] = mean(fdr[
      abs(test_statistics) >= abs(test_statistics[i]) # more extreme
    ])
  }
  return(Fdr)
}

p_from_t <- function(test_statistics, df = NULL, sides = 'two') {
  if (is.null(df)) {
    # assume standard normal
    one_sided <- unlist(lapply(test_statistics, function(z) {
      pnorm(-1*abs(z))
    }))
  } else {
    # assume t_df
    one_sided <-unlist(lapply(test_statistics, function(z) {
      pt(-1*abs(z), df = df)
    }))
  }
  
  if (sides == 'two') {
    return (2*one_sided)
  } else {
    return(one_sided)
  }
}

nrow_null0 <- function(dataframe) {
  ifelse(is.null(dataframe), 0, nrow(dataframe))
}