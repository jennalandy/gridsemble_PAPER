source("./simulate.R")
source("./test.R")

run_inc_n_test(
  data_generating_fn = function() {
    sim_simple(
      m = 1000,
      pi0 = 0.8,
      mu0 = 0,
      sigma0 = 3,
      min1 = 4, 
      max1 = 12
    )
  }, 
  N = 200, 
  path = "inc_n_symmetric",
  overwrite = TRUE
)
