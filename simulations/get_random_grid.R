# random grids don't have all combinatorial combinations
# each parameter is randomly sampled separately
# we match max grid sizes of gridsemble with grid search

fdrtool_grid_size = 22
fdrtool_grid <- data.frame(
  cutoff.method = c("fndr", "locfdr", rep("pct0", fdrtool_grid_size - 2)),
  pct0 = runif(fdrtool_grid_size, 0.4, 1)
)

locfdr_grid_size = 150
locfdr_grid <- data.frame(
  pct = runif(locfdr_grid_size, 0, 0.2),
  pct0 = runif(locfdr_grid_size, 0, 0.3),
  nulltype = sample(c(1,2,3), size = locfdr_grid_size, replace = TRUE),
  type = sample(c(0, 1), size = locfdr_grid_size, replace = TRUE)
)

qvalue_grid_size = 120
qvalue_grid <- data.frame(
  transf = sample(c("probit", "logit"), size = qvalue_grid_size, replace = TRUE),
  adj = runif(qvalue_grid_size, 0.5, 2),
  pi0.method = sample(c("bootstrap","smoother"), size = qvalue_grid_size, replace = TRUE),
  smooth.log.pi0 = sample(c(TRUE, FALSE), size = qvalue_grid_size, replace = TRUE)
)
