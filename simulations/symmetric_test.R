source("./simulate.R")
source("./test.R")

run_test(
  data_generating_fn = function(i) {
    sim_simple(
      m = 1000, # number of tests
      pi0 = 0.8,
      mu0 = 0,
      sigma0 = 3,
      min1 = 4, 
      max1 = 12
    )
  }, 
  N = 200, 
  path = "./symmetric",
  overwrite = TRUE
)
