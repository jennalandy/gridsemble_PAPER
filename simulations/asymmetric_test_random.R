source("./simulate.R")
source("./test.R")
source("./get_random_grid.R")

load("./asymmetric/all_datas.RData") # loads `all_datas` to environment
run_test(
  data_generating_fn = function(i) {
      all_datas[[i]] # use same datasets from asymmetric for pairwise comparisons
  }, 
  N = length(all_datas), 
  path = "asymmetric_random",
  overwrite = TRUE,
  grids = list(
      locfdr = locfdr_grid,
      fdrtool = fdrtool_grid,
      qvalue = qvalue_grid
  )
)
