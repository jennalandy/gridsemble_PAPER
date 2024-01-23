source("./simulate.R")
source("./test.R")
source("./get_random_grid.R")

load("./cod_based/all_datas.RData") # loads `all_datas` to environment

run_test(
  data_generating_fn = function(i) {
      all_datas[[i]] # use same datasets from cod-based for pairwise comparisons
  }, 
  N = length(all_datas), 
  path = "cod_based_random",
  overwrite = TRUE,
  sim_n = 1000,
  df = 576,
  grids = list(
      locfdr = locfdr_grid,
      fdrtool = fdrtool_grid,
      qvalue = qvalue_grid
  )
)
