source("./simulate.R")
source("./test.R")

run_test(
  data_generating_fn = function(i) {
    sim_cod_based(
      minimum = 0.02,
      maximum = 0.04
    )
  }, 
  N = 200, 
  path = "cod_based",
  overwrite = TRUE,
  sim_n = 1000,
  df = 576
)
