library(snow)
library(futile.logger)
library(BiocParallel)
library(gridsemblefdr)
library(cli)
library(reshape2)
library(fdrtool)
library(withr)

source("./evaluate.R")

run_test <- function(
    data_generating_fn, 
    N, 
    path,
    ensemble_size = 10, 
    n_synthetic = 10,
    overwrite = FALSE, 
    sim_n = NULL,
    parallel = FALSE,
    grids = NULL,
    df = NULL
) {
  if (!endsWith(path, '/')) {
    path = paste(path, '/', sep = '')
  }
  
  if (!dir.exists(path)) {
    dir.create(path)
  } else if (overwrite) {
    cat('\nDirectory', path, 'already exists')
    cat('\nOverwritting because given overwrite = TRUE')
  } else {
    cat('\nDirectory', path, 'already exists')
    cat('\nEnding tests because given overwrite = FALSE')
    cat('\nRerun with overwrite = TRUE or choose a different directory\n\n')
    return()
  }
  
  sink(file = paste(path, 'tests.log', sep = ''))
  
  start = Sys.time()
  cat(paste("Starting Test at", start))
  rep = 1
  
  all_datas <- list()
  while (rep <= N) {
    cat(paste("\n\nRep", rep, '/', N))
    tryCatch({
      dat <- data_generating_fn(rep)
      all_datas[[rep]] <- dat
      if (is.null(sim_n)) {
        sim_n = length(dat$zz)
      }
      
      if (is.null(grids)) {
        locfdr_grid <- build_locfdr_grid(dat$zz)
        fdrtool_grid <- build_fdrtool_grid(dat$zz)
        qvalue_grid <- build_qvalue_grid(dat$zz)
      } else {
        locfdr_grid <- reduce_locfdr_grid(dat$zz, grids$locfdr)
        fdrtool_grid <- reduce_fdrtool_grid(dat$zz, grids$fdrtool)
        qvalue_grid <- reduce_qvalue_grid(dat$zz, to_pval_function = p_from_t, grids$qvalue)
      }
     
      grid_size = nrow_null0(locfdr_grid) +
        nrow_null0(fdrtool_grid) +
        nrow_null0(qvalue_grid)
      cat(paste('\n\tgrid grid_size =', grid_size))
      
      if (ensemble_size > grid_size) {
        cat('\n\tensemble size > grid size, skipping')
        stop()
      }
      
      topq <- abs(dat$zz) > quantile(abs(dat$zz), 0.75)
      
      cat("\n\tsetup done")
      # ------- RUN ALL METHODS ------- 
      
      # default runs of other methods
      locfdr_run <- locfdr::locfdr(dat$zz, plot = F)
      locfdr_run$Fdr <- Fdr_from_fdr(locfdr_run$fdr, dat$zz)
     
      fdrtool_run <- fdrtool::fdrtool(dat$zz, plot = F, verbose = F)
      fdrtool_run$Fdr <- Fdr_from_fdr(fdrtool_run$lfdr, dat$zz)
      
      # if qvalue fails, rerun with lambda = 0
      qvalue_run <- NULL
      tryCatch({
        qvalue_run <- qvalue::qvalue(
          p_from_t(dat$zz, sides = 'two'), plot = 0, verbose = F
        )
      }, error = function(e) {})
      if (is.null(qvalue_run)) {
        qvalue_run <- qvalue::qvalue(
          p_from_t(dat$zz, sides = 'two'), plot = 0, lambda = 0, verbose = F
        )
      }
      qvalue_run$Fdr <- Fdr_from_fdr(qvalue_run$lfdr, dat$zz)
      
      cat("\n\t\tdefaults run")
      
      ensemble_run <- gridsemble(
        dat$zz, n_synthetic = 0, n_workers = 5,
        ensemble_size = ensemble_size, verbose = F,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df
      )
      cat("\n\t\tensemble run")
      
      ensemble_all_run <- gridsemble(
        dat$zz, n_synthetic = 0, n_workers = 5,
        ensemble_size = grid_size, verbose = F,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df
      )
      
      cat("\n\t\tensemble all run")
      
      grid_run <- gridsemble(
        dat$zz, n_synthetic = n_synthetic, 
        ensemble_size = 1,
        verbose = F, n_workers = 5,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df
      )
      
      cat("\n\t\tgrid run")
      
      gridsemble_run <- gridsemble(
        dat$zz, n_synthetic = n_synthetic,
        ensemble_size = ensemble_size,
        verbose = F, 
        n_workers = 5,
        locfdr_grid = locfdr_grid,
        fdrtool_grid = fdrtool_grid,
        qvalue_grid = qvalue_grid,
        parallel = parallel, 
        synthetic_size = sim_n,
        df = df
      )
      
      cat("\n\t\tgridsemble run")
      
      cat("\n\tall methods run")
      
      # ------- SAVE RESULTS ------- 
      
      generate_params_i <- gridsemble_run$generating_model$parameters
      
      pi0_estimates_i <- list(
        "locfdr" = locfdr_run$fp0["mlest", "p0"],
        "fdrtool" = fdrtool_run$param[1, 'eta0'],
        "qvalue" = qvalue_run$pi0,
        "ensemble" = ensemble_run$pi0,
        "ensemble_all" = ensemble_all_run$pi0,
        "grid" = grid_run$pi0,
        "gridsemble" = gridsemble_run$pi0
      )
      
      pi0_var_i <- list(
        "ensemble" = ensemble_run$pi0_var,
        "ensemble_all" = ensemble_all_run$pi0_var,
        "gridsemble" = gridsemble_run$pi0_var
      )
      
      fdr_var_i <- list(
        "ensemble" = median(ensemble_run$fdr_var),
        "ensemble_all" = median(ensemble_all_run$fdr_var),
        "gridsemble" = median(gridsemble_run$fdr_var)
      )
      
      locfdr_cutoff <- quantile(locfdr_run$fdr, 1 - min(pi0_estimates_i$locfdr, 1))
      locfdr_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = locfdr_run$fdr, 
        Fdr = locfdr_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = locfdr_cutoff
      )
      
      fdrtool_cutoff <- quantile(fdrtool_run$lfdr, 1 - min(pi0_estimates_i$fdrtool, 1))
      fdrtool_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = fdrtool_run$lfdr, 
        Fdr = fdrtool_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = fdrtool_cutoff
      )
      
      qvalue_cutoff <- quantile(qvalue_run$lfdr, 1 - min(pi0_estimates_i$qvalue, 1))
      qvalue_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = qvalue_run$lfdr, 
        Fdr = qvalue_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = qvalue_cutoff
      )
      
      ensemble_cutoff <- quantile(ensemble_run$fdr, 1 - min(pi0_estimates_i$ensemble))
      ensemble_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = ensemble_run$fdr, 
        Fdr = ensemble_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = ensemble_cutoff
      )
      
      ensemble_all_cutoff <- quantile(ensemble_all_run$fdr, 1 - min(pi0_estimates_i$ensemble_all))
      ensemble_all_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = ensemble_all_run$fdr, 
        Fdr = ensemble_all_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = ensemble_all_cutoff
      )
      
      grid_cutoff <- quantile(grid_run$fdr, 1 - min(pi0_estimates_i$grid))
      grid_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = grid_run$fdr, 
        Fdr = grid_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = grid_cutoff
      )
      
      gridsemble_cutoff <- quantile(gridsemble_run$fdr, 1 - min(pi0_estimates_i$gridsemble))
      gridsemble_metrics_i <- calc_metrics(
        test_statistics = dat$zz,
        fdr = gridsemble_run$fdr, 
        Fdr = gridsemble_run$Fdr, 
        truth = dat$truth, 
        true_Fdr = dat$true_Fdr,
        true_fdr = dat$true_fdr,
        topq = topq,
        cutoff = gridsemble_cutoff
      )
      
      gridsemble_package_counts_i = gridsemble_run$top_grid %>%
        mutate(method = factor(method, levels = c('qvalue','locfdr','fdrtool'))) %>%
        pull(method) %>%
        table()
      gridsemble_package_counts_i$i = rep

      if (rep == 1) {
        generate_params <- data.frame(generate_params_i)
        pi0_estimates <- data.frame(pi0_estimates_i)
        pi0_var <- data.frame(pi0_var_i)
        fdr_var <- data.frame(fdr_var_i)
        
        locfdr_metrics <- data.frame(locfdr_metrics_i)
        fdrtool_metrics <- data.frame(fdrtool_metrics_i)
        qvalue_metrics <- data.frame(qvalue_metrics_i)
        
        ensemble_metrics <- data.frame(ensemble_metrics_i)
        ensemble_all_metrics <- data.frame(ensemble_all_metrics_i)
        grid_metrics <- data.frame(grid_metrics_i)
        gridsemble_metrics <- data.frame(gridsemble_metrics_i)
        
        gridsemble_package_counts <- data.frame(gridsemble_package_counts_i)
      } else {
        generate_params <- rbind(
          generate_params,
          generate_params_i
        )
        pi0_estimates <- rbind(
          pi0_estimates, 
          pi0_estimates_i
        )
        pi0_var <- rbind(
          pi0_var,
          pi0_var_i
        )
        fdr_var <- rbind(
          fdr_var,
          fdr_var_i
        )
        
        locfdr_metrics <- rbind(
          locfdr_metrics, 
          locfdr_metrics_i
        )
        fdrtool_metrics <- rbind(
          fdrtool_metrics,
          fdrtool_metrics_i
        )
        qvalue_metrics <- rbind(
          qvalue_metrics,
          qvalue_metrics_i
        )
        
        ensemble_metrics <- rbind(
          ensemble_metrics,
          ensemble_metrics_i
        )
        ensemble_all_metrics <- rbind(
          ensemble_all_metrics,
          ensemble_all_metrics_i
        )
        grid_metrics <- rbind(
          grid_metrics,
          grid_metrics_i
        )
        gridsemble_metrics <- rbind(
          gridsemble_metrics,
          gridsemble_metrics_i
        )
        
        gridsemble_package_counts <- rbind(
          gridsemble_package_counts,
          gridsemble_package_counts_i
        )
      }
      
      cat("\n\tmetrics computed")
      
      # ------- SAVE ALL DATAFRAMES ON EACH ITER ------- 
      write.csv(generate_params, paste0(path, 'generate_params', '.csv'), row.names = FALSE)
      write.csv(pi0_estimates, paste0(path, 'pi0_estimates', '.csv'), row.names = FALSE)
      write.csv(pi0_var, paste0(path, 'pi0_var', '.csv'), row.names = FALSE)
      write.csv(fdr_var, paste0(path, 'fdr_var', '.csv'), row.names = FALSE)
      
      write.csv(locfdr_metrics, paste0(path, 'locfdr_metrics', '.csv'), row.names = FALSE)
      write.csv(fdrtool_metrics, paste0(path, 'fdrtool_metrics', '.csv'), row.names = FALSE)
      write.csv(qvalue_metrics, paste0(path, 'qvalue_metrics', '.csv'), row.names = FALSE)

      write.csv(ensemble_metrics, paste0(path, 'ensemble_metrics', '.csv'), row.names = FALSE)
      write.csv(ensemble_all_metrics, paste0(path, 'ensemble_all_metrics', '.csv'), row.names = FALSE)
      write.csv(grid_metrics, paste0(path, 'grid_metrics', '.csv'), row.names = FALSE)
      write.csv(gridsemble_metrics, paste0(path, 'gridsemble_metrics', '.csv'), row.names = FALSE)
      
      write.csv(gridsemble_package_counts, paste0(path, 'gridsemble_package_counts', '.csv'), row.names = FALSE)
      
      save(all_datas, file = paste0(path, 'all_datas', '.RData'))
      
      cat("\n\tdataframes saved")
      cat(paste("\nFinished", rep, '/', N, ': ', Sys.time()))
      
      rep = rep + 1
    },
    error = function(e) {
      cat(paste('\nError on', rep,'\n'))
      cat(paste("\n", e))
    })
  }
  
  now <- Sys.time()
  cat(paste('\n\nDone in', round(difftime(now, start, units = 'hours'), 2), 'hours'))
}



run_inc_n_test <- function(
    data_generating_fn, 
    N, 
    path,
    overwrite = FALSE,
    df = NULL
) {
  if (!endsWith(path, '/')) {
    path = paste(path, '/', sep = '')
  }
  
  if (!dir.exists(path)) {
    dir.create(path)
  } else if (overwrite) {
    cat('\nDirectory', path, 'already exists')
    cat('\nOverwritting because given overwrite = TRUE')
  } else {
    cat('\nDirectory', path, 'already exists')
    cat('\nEnding tests because given overwrite = FALSE')
    cat('\nRerun with overwrite = TRUE or choose a different directory\n\n')
    return()
  }
  
  sink(file = paste(path, 'tests.log', sep = ''))
  
  n_synthetic_values = c(0, 1, 5, 10, 20)
  ensemble_size_values = c(1, 5, 10, 20)
  
  start <- Sys.time()
  cat(paste("Starting test at", start))
  rep = 1
  while (rep <= N) {
    cat(paste('\n\nRep', rep, '/', N))
    
    dat <- data_generating_fn()
    
    locfdr_grid <- build_locfdr_grid(dat$zz)
    fdrtool_grid <- build_fdrtool_grid(dat$zz)
    qvalue_grid <- build_qvalue_grid(dat$zz)
    grid_size = nrow_null0(locfdr_grid) +
      nrow_null0(fdrtool_grid) +
      nrow_null0(qvalue_grid)
    cat(paste('\n\tGrid_size =', grid_size))
    
    if (grid_size < 20) {
      cat('\ngridsize too small (<20), skipping dataset')
      next
    }
    
    true_Fdr = dat$Fdr
    true_fdr = dat$fdr
    topq = abs(dat$zz) > quantile(abs(dat$zz), 0.75)
    
    for (n_synthetic in n_synthetic_values) {
      cat(paste('\n\tn_synthetic =', n_synthetic))
      
      for (ensemble_size in ensemble_size_values) {
        cat(paste('\n\t\tensemble_size =', ensemble_size))
        
        gridsemble_run <- gridsemble(
          dat$zz, 
          n_synthetic = n_synthetic,
          ensemble_size = ensemble_size,
          verbose = FALSE, 
          parallel = FALSE,
          locfdr_grid = locfdr_grid,
          fdrtool_grid = fdrtool_grid,
          qvalue_grid = qvalue_grid,
          df = df
        )
        cat('\n\t\t\tgridsemble done')
        
        ensemble_run <-  gridsemble(
          dat$zz, 
          n_synthetic = 0,
          ensemble_size = ensemble_size, 
          verbose = FALSE,
          parallel = FALSE,
          locfdr_grid = locfdr_grid,
          fdrtool_grid = fdrtool_grid,
          qvalue_grid = qvalue_grid,
          df = df
        )
        cat('\n\t\t\tensemble done')
        
        gridsemble_cutoff <- quantile(gridsemble_run$fdr, 1 - min(gridsemble_run$pi0))
        gridsemble_metrics_i <- calc_metrics(
          test_statistics = dat$zz,
          fdr = gridsemble_run$fdr, 
          Fdr = gridsemble_run$Fdr, 
          truth = dat$truth, 
          true_Fdr = dat$true_Fdr,
          true_fdr = dat$true_fdr,
          topq = topq,
          cutoff = gridsemble_cutoff
        )
        gridsemble_metrics_i$n_synthetic = n_synthetic
        gridsemble_metrics_i$ensemble_size = ensemble_size
        
        ensemble_cutoff <- quantile(ensemble_run$fdr, 1 - min(ensemble_run$pi0))
        ensemble_metrics_i <- calc_metrics(
          test_statistics = dat$zz,
          fdr = ensemble_run$fdr, 
          Fdr = ensemble_run$Fdr, 
          truth = dat$truth, 
          true_Fdr = dat$true_Fdr,
          true_fdr = dat$true_fdr,
          topq = topq,
          cutoff = ensemble_cutoff
        )
        ensemble_metrics_i$n_synthetic = n_synthetic
        ensemble_metrics_i$ensemble_size = ensemble_size
        
        if (
          rep == 1 &
          n_synthetic == n_synthetic_values[1] & 
          ensemble_size == ensemble_size_values[1]
        ) {
          gridsemble_metrics <- data.frame(gridsemble_metrics_i)
          ensemble_metrics <- data.frame(ensemble_metrics_i)
        } else {
          gridsemble_metrics <- rbind(
            gridsemble_metrics,
            gridsemble_metrics_i
          )
          ensemble_metrics <- rbind(
            ensemble_metrics, 
            ensemble_metrics_i
          )
        }
        
        write.csv(
          gridsemble_metrics, 
          paste0(path, "gridsemble_metrics", ".csv")
        )
        write.csv(
          ensemble_metrics, 
          paste0(path, "ensemble_metrics", ".csv")
        )
        cat("\n\t\t\tdataframes saved")
      }
    }
    cat(paste("\n\t\tFinished", rep, '/', N, ': ', Sys.time()))
    rep = rep + 1
  }
  now <- Sys.time()
  cat(paste('\n\nDone in', round(difftime(now, start, units = 'hours'), 2), 'hours'))
  sink()
}