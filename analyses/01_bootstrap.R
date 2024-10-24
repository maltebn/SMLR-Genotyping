library(tidyverse)
library(future.apply)
source("00_functions_and_global_definitions.R")

workers <- floor(availableCores()[[1]]*3/4)
plan(multisession, workers = workers)
reps <- 10^3

# dd_GT <- readRDS(file.path("..", "data", "00_dd_GT.rds"))
dd_GT_2nd <- readRDS(file.path("..", "data", "00_dd_GT_2nd.rds"))
# dd_GT_both <- dd_GT_2nd |> bind_rows(dd_GT[names(dd_GT_2nd)])
true_profiles_2nd <- readRDS(file.path("..", "data", "00_true_profiles_2nd.rds"))


dils <- list(50, c(50, 25), 25, 12.5, 6.25)
set.seed(1913)
global_seeds <- sample.int(10^6, size=length(dils))
# Note: you have to call plan(multisession, workers = workers) before
# future_lapply() will work correctly, and thus also before the generation of
# the future_seeds below (forgetting to set up a plan may result in an identical
# seed being used in all parallel sessions of future_lapply()).
future_seeds <- lapply(global_seeds, function(s) {
  future_lapply(seq_len(reps), FUN = function(x) {
    .Random.seed
  }, future.chunk.size = Inf, future.seed = s)
})

for (vst in c("sqrt", "log", "I")) {
  if (vst == "I") {
    initial_beta <- c(0, 0, -1)
  } else {
    initial_beta <- c(0, 1, -2)
  }
  for (intercept in c(T, F)) {
    if (intercept) {
      bootstrap_name <- paste0("01_obj_bootstrap_", vst)
    } else {
      initial_beta <- initial_beta[-1]
      bootstrap_name <- paste0("01_obj_bootstrap_", vst, "_NI")
    }
    
    dd_boot_fine <- lapply(seq_along(dils), function(i) {
      dd_GT_2nd |> filter(Dilution %in% dils[[i]]) |> 
        dd_boot_create_samplesize(reps, stratify = "Marker",
                                  f=vst, INT=intercept, b_int = initial_beta,
                                  method = "Nelder-Mead",
                                  control_list = list(maxit = 10^3, reltol=sqrt(.Machine$double.eps)), seeds = future_seeds[[i]])
    })
    
    dd_boot_fine <- dd_boot_fine |> 
      lapply(function(d){
        lapply(d, function(r){
          lapply(r, function(s){
            p <- s[[1]]
            s$beta0 <- p[1]
            s$beta1 <- p[2]
            s$beta2 <- p[3]
            s[-1]
          }) |> bind_rows()
        }) |> bind_rows()
      }) |> bind_rows()
    dd_boot_fine$Converged <- as.integer(dd_boot_fine$n_WC != 0)
    
    saveRDS(dd_boot_fine, file.path("..", "data", paste0(bootstrap_name, "_1K.rds")))
    gc()
  }
}