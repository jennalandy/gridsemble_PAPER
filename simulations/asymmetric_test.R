source("./simulate.R")
source("./test.R")

run_test(
  data_generating_fn = function(i) {
    sim_simple_asym(
      m = 1000, 
      pi0 = 0.85, 
      pi1L = 1/3, 
      mu0 = 0, 
      sigma0 = 2,
      min1L = -12, 
      max1L = -5, 
      min1R = 3, 
      max1R = 9, 
      seed = NA
    )
  }, 
  N = 200, 
  path = "./asymmetric",
  overwrite = TRUE
)
