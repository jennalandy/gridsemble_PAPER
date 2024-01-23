source("./utilities.R")
library(dplyr)
library(BiocGenerics)
library(Biobase)
library(curatedOvarianData)

rowttests <- function(matrix, condition) {
  condition.first = unique(condition)[1]
  data.frame(
    statistic = sapply(1:nrow(matrix), function(i) {
      row = matrix[i,]
      row1 = row[condition == condition.first]
      row2 = row[condition != condition.first]
      row_test = t.test(row1, row2, var.equal=TRUE)
      unname(row_test$statistic)
    })
  )
}

sim_cod_based <- function(
    pi0 = 0.8, 
    minimum = 0.01, 
    maximum = 0.02, 
    mult = 1,
    skew = 0  # positive skews left (negative skew)
    # neg skews right (positive skew)
) {
  
  gc()
  
  if (!exists('TCGA_eset')) {
      data(TCGA_eset)
  }
  expression = exprs(TCGA_eset)
  labels = TCGA_eset$recurrence_status
  
  pi1 = 1-pi0
  
  
  # 1. deal with missing labels
  # to keep as much data as possible, add random labels for those with missing values
  # (we will permute labels anyways, so it doesn't matter)
  labels[is.na(labels)] <- sample(
    labels[!is.na(labels)],    # sample from existing labels to keep proportions constant
    size = sum(is.na(labels)),
    replace = T
  )
  
  expression = unname(expression)
  
  # 2. deal with missing expression data
  # fill NA with mean expression for the gene
  if (sum(is.na(expression)) > 0 ) {
    k <- which(is.na(expression), arr.ind=TRUE)
    expression[k] <- rowMeans(expression, na.rm = TRUE)[k[,1]]
  }
  
  # 3. Get rid of genes with no variation
  geneVars <- expression %>% apply(MARGIN = 1, FUN = function(x) {sd(x)})
  expression <- expression[geneVars != 0,]
  
  # 3. Get labels and permute them to break any association between expression
  # and label while preserving correlation structure of the Xs
  y <- ifelse(labels == unique(labels)[1], 1, 0)
  new_y <- sample(y, size = length(y))
  
  n = nrow(expression)
  m = ncol(expression)  
  
  # 4. Choose a subset of genes (proportion pi1 = 1- pi0)
  
  # shuffle expression rowsexpression <- expression[sample(1:n),]
  ndiff <- max(round(n * pi1), 1)
  
  # then first pi1*100% have DE
  # store truth = 1 if different, 0 otherwise
  truth <- c(rep(1, ndiff), rep(0, n - ndiff))
  
  
  
  # 5. Get and store initial t statistics for permuted data. 
  t_init <- rowttests(expression, factor(new_y))
  
  
  # Using the first gene, double check the direction in which
  # t-stats were calculated (group 0 - group 1 or vise versa).
  # store first = 0 or 1, where "first" is the first number
  # written in the subtraction statement
  group_means <- data.frame(e = expression[1,], new_y = new_y) %>% 
    group_by(new_y) %>% 
    summarize_all(mean)
  
  smaller <- group_means$new_y[group_means$e == min(group_means$e)]
  
  if (t_init$statistic[1] < 0) {
    first = smaller
  } else {
    first = as.numeric(!smaller)
  }
  
  # 6. Change expression level for "diff" genes where y = 1
  # strategically multiplying in direction of existing differences
  new_expression <- expression
  
  # diff means is the multiplicative difference we'll force upon 
  # the means within not-null tests. length = num not null.
  # these values will be multiplied by the expression values
  # where truth = 1 and group 1
  diff_means <- sapply(t_init$statistic[truth==1], function(t) {
    minimum_n <- minimum
    maximum_n <- maximum
    
    # subtracting group 0 - group 1
    if (first == 0) {
      if (t < 0) {
        # if t is less than 0, we want to inflate the expression values where y = 1
        # this will decrease mean0 - mean1 further negative
        return(runif(n = 1, 1 + skew*minimum + minimum_n, 1 + skew*minimum + maximum_n))
        
      } else {
        # else, increase mean0 - mean1 further positive
        return(runif(n = 1, 1 - maximum_n, 1 - minimum_n))
      }
      
      # subtracting group 1 - group 0
    } else {
      if (t > 0) {
        # if t is greater than 0, we want to deflate the expression values where y = 1
        # this will decrease mean1 - mean0 further negative
        return(runif(n = 1, 1 - maximum_n, 1 - minimum_n))
      } else {
        # otherwise increase mean1 - mean0 further positive
        return(runif(n = 1, 1 + skew*minimum + minimum_n, 1 + skew*minimum + maximum_n))
      }
    }
    
  })
  
  # append 1 values for the null tests
  diff_means <- c(diff_means, rep(1, n-ndiff))
  
  # create a matrix with all columns as diff means
  # for group 1, all 1s for group 0
  diff_mat <- matrix(
    rep(diff_means, m), 
    ncol = m
  )
  diff_mat[,!new_y] <- 1
  
  # multiply expression by diff matrix
  new_expression <- new_expression * diff_mat
  
  t <- rowttests(new_expression, factor(new_y))
  
  return(list(
    zz = t$statistic,
    truth = truth,
    true_Fdr = get_empirical_Fdr(t$statistic, truth),
    true_fdr = NULL,
    expression = expression,
    y = y,
    params = list(
      'pi0' = pi0,
      'minimum' = minimum,
      'maximum' = maximum
    )
  ))
}

sim_simple_asym <- function(
    m = 1000, 
    pi0 = 0.85, 
    pi1L = 1/3, 
    mu0 = 0, 
    sigma0 = 2,
    min1L = -10, 
    max1L = -4, 
    min1R = 3, 
    max1R = 6, 
    seed = NA
) {
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n_null = round(m*pi0)
  n_alt = m - n_null
  
  zz <- c(
    rnorm(n = n_null, mean = mu0, sd = sigma0), 
    runif(n = round(n_alt*pi1L), min1L, max1L), 
    runif(n = round(n_alt*(1 - pi1L)), min1R, max1R)
  )
  truth <- c(rep(0, round(m*pi0)), rep(1, m - round(m*pi0)))
  
  # null pdf
  null <- function(z) {
    (1/sqrt(2*pi*sigma0^2))*exp(-(z-mu0)^2/(2*sigma0^2))
  }
  
  # alt pdf
  alt <- function(zz) {
    # U(min, max)/2 + U(-max, -min)/2
    # U(min, max) ~ I(min <= z <= max)*1/(max-min)
    unlist(lapply(zz, function(z) {
      if ( min1L <= z & z <= max1L) {
        return(pi1L/(max1L - min1L))
      } else if (min1R <= z & z <= max1R) {
        return((1 - pi1L)/(max1R - min1R))
      } else {
        return(0)
      }
    }))
  }
  
  # mixture distribution
  mix <- function(zz) {
    pi0*null(zz) + (1-pi0)*alt(zz)
  }
  
  # fdr = pi0f0(z)/f(z)
  fdr <- pi0*null(zz)/mix(zz)
  
  Fdr <- Fdr_from_fdr(fdr, zz)
  
  
  shuffle <- sample(1:m)
  
  return(list(
    'zz' = zz[shuffle],
    'truth' = truth[shuffle],
    'true_Fdr' = Fdr[shuffle],
    'true_fdr' = fdr[shuffle],
    'params' = list(
      'pi0' = pi0, 
      'pi1L' = pi1L, 
      'mu0' = mu0, 
      'sigma0' = sigma0,
      'min1L' = min1L, 
      'max1L' = max1L, 
      'min1R' = min1R, 
      'max1R' = max1R
    )
  ))
}

sim_simple <- function(
    m = 1000,   # number of tests
    pi0 = 0.8,  # proportion null
    mu0 = 0,    # mean of null distribution
    sigma0 = 3, # sd of null distribution
    min1 = 2,   # min of +/- uniform alternative
    max1 = 6,  # max of +/- uniform alternative
    seed = NA   # random seed
) {
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n_null = round(m*pi0)
  n_alt = m - n_null
  
  zz <- c(
    rnorm(n = n_null, mean = mu0, sd = sigma0),   # null distribution
    runif(n = floor(n_alt/2), min1, max1),        # positive of alt
    runif(n = ceiling(n_alt/2), -1*max1, -1*min1) # negative of alt
  )
  truth <- c(rep(0, n_null), rep(1, n_alt))
  
  # null pdf
  null <- function(z) {
    (1/sqrt(2*pi*sigma0^2))*exp(-(z-mu0)^2/(2*sigma0^2))
  }
  
  # alt pdf
  alt <- function(zz) {
    # U(min, max)/2 + U(-max, -min)/2
    # U(min, max) ~ I(min <= z <= max)*1/(max-min)
    unlist(lapply(zz, function(z) {if (
      (-1*max1 <= z & z <= -1*min1) |
      (min1 <= z & z <= max1)
    ) {
      return(1/(2*(max1 - min1)))
    } else {
      return(0)
    }}))
  }
  
  # mixture distribution
  mix <- function(zz) {
    pi0*null(zz) + (1-pi0)*alt(zz)
  }
  
  # fdr = pi0f0(z)/f(z)
  fdr <- pi0*null(zz)/mix(zz)

  Fdr <- Fdr_from_fdr(fdr, zz)
  
  # shuffle the dataset
  shuffle <- sample(1:m)
  
  return(list(
    'zz' = zz[shuffle],
    'truth' = truth[shuffle],
    'true_fdr' = fdr[shuffle],
    'true_Fdr' = Fdr[shuffle],
    'params' = list(
      'sigma0' = sigma0,
      'pi0' = pi0,
      'mu0' = mu0,
      'min1' = min1,
      'max1' = max1
    )
  ))
}