library(tidyverse)
library(future.apply)
source("00_functions_and_global_definitions.R")

workers <- floor(availableCores()[[1]]*3/4)
plan(multisession, workers = workers)

reps <- 10^3
train_proportion <- 9/12
set.seed(123)


dd_GT <- readRDS(file.path("..", "data", "00_dd_GT.rds"))
dd_GT_2nd <- readRDS(file.path("..", "data", "00_dd_GT_2nd.rds"))
true_profiles_2nd <- readRDS(file.path("..", "data", "00_true_profiles_2nd.rds"))
# We will treat the different runs in dd_GT as independent, i.e. simply consider
# a new run as new individuals which we encode by interaction(Run, ID, drop=T)
dd_GT$ID <- interaction(dd_GT$ID, dd_GT$Run, drop = TRUE)

dd_GT_both <- dd_GT_2nd |> bind_rows(dd_GT[names(dd_GT_2nd)])

# Test on relevant low dilutions
# (where the cases with missing heterozygous allele signals isn't too numerous)
dd_dils_FSIGEN_l <- list(list("dil_train" = c(250, 125, 62.5), "dil_test" = c(50, 31.25)),
                         list("dil_train" = c(50, 31.25), "dil_test" = c(50, 31.25)),
                         list("dil_train" = c(125, 25), "dil_test" = c(50, 31.25)),
                         list("dil_train" = c(50, 25), "dil_test" = c(50, 31.25)),
                         list("dil_train" = c(31.25, 25), "dil_test" = c(50, 31.25)),
                         list("dil_train" = 50, "dil_test" = c(50, 31.25)),
                         list("dil_train" = 31.25, "dil_test" = c(50, 31.25)),
                         list("dil_train" = 25, "dil_test" = c(50, 31.25)),
                         list("dil_train" = 1000, "dil_test" = c(50, 31.25)),
                         list("dil_train" = c(1000, 25), "dil_test" = c(50, 31.25)))

# dd_dils <- list(list("dil_train" = c(50, 25), "dil_test" = 1000),
#                 list("dil_train" = c(50, 25), "dil_test" = 500),
#                 list("dil_train" = c(50, 25), "dil_test" = 250),
#                 list("dil_train" = c(50, 25), "dil_test" = 125),
#                 list("dil_train" = c(50, 25), "dil_test" = 62.5),
#                 list("dil_train" = c(50, 25), "dil_test" = 50),
#                 list("dil_train" = c(50, 25), "dil_test" = 31.25),
#                 list("dil_train" = c(50, 25), "dil_test" = 25))

# dd_dils <- list(list("dil_train" = 1000, "dil_test" = 1000),
#                 list("dil_train" = 500, "dil_test" = 500),
#                 list("dil_train" = 250, "dil_test" = 250),
#                 list("dil_train" = 125, "dil_test" = 125),
#                 list("dil_train" = 62.5, "dil_test" = 62.5),
#                 list("dil_train" = 50, "dil_test" = 50),
#                 list("dil_train" = 31.25, "dil_test" = 31.25),
#                 list("dil_train" = 25, "dil_test" = 25))

# dd_dils <- list(list("dil_train" = c(50, 31.25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(31.25, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(1000, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(50, 25), "dil_test" = 62.5),
#                 list("dil_train" = c(50, 25), "dil_test" = 50),
#                 list("dil_train" = c(50, 25), "dil_test" = 31.25),
#                 list("dil_train" = c(50, 25), "dil_test" = 25))

# Test on relevant high dilutions
# (where the HSG makes WCs and the number of NCs for EQC increases rapidly)
dd_dils_FSIGEN_h <- list(list("dil_train" = c(250, 125, 62.5), "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = c(50, 31.25), "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = c(125, 25), "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = c(50, 25), "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = c(31.25, 25), "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = 50, "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = 31.25, "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = 25, "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = 1000, "dil_test" = c(250, 125, 62.5)),
                         list("dil_train" = c(1000, 25), "dil_test" = c(250, 125, 62.5)))

# dd_dils <- list(list("dil_train" = c(50, 31.25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(31.25, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(1000, 25), "dil_test" = c(250, 125, 62.5)))

# Test on extreme dilutions with the good performing fits from the above low and high dilutions
# dd_dils <- list(list("dil_train" = c(1000, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(500, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(250, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(125, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(62.5, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(1000, 500)),
#                 list("dil_train" = c(1000, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(500, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(250, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(125, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(62.5, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(25, 12.5, 6.25)),
#                 list("dil_train" = c(1000, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(500, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(250, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(125, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(62.5, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(250, 125, 62.5)),
#                 list("dil_train" = c(1000, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(500, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(250, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(125, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(62.5, 25), "dil_test" = c(50, 31.25)),
#                 list("dil_train" = c(50, 25), "dil_test" = c(50, 31.25)))

dd_dils_FSIGEN_lh <- list(list("dil_train" = c(50, 25), "dil_test" = c(50, 31.25)),
                     list("dil_train" = c(50, 25), "dil_test" = c(250, 125, 62.5)),
                     list("dil_train" = c(50, 25), "dil_test" = c(1000, 500)))


for (dil_range in c("low", "high", "all")) {
  cv_fun <- cross_val_create_custom
  
  for (vst in c("sqrt", "log", "I")) {
    initial_beta <- c(0, 1, -2)
    if (vst == "sqrt") {
      f_vst <- sqrt
    } else if (vst == "log") {
      f_vst <- function(x) log(x+1)
    } else if (vst == "I") {
      initial_beta <- c(0, 0, -1)
      f_vst <- function(x) x
    }
    for (INT in c(T, F)) {
      I_NI <- "_"
      if (!INT) {
        initial_beta <- initial_beta[-1]
        I_NI <- "_NI_"
      }
      
      if (dil_range == "low") {
        filename <- paste0("01_obj_crossval_", vst, I_NI, "1K_9individuals_low.rds")
        dd_dils_FSIGEN <- dd_dils_FSIGEN_l
      } else if (dil_range == "high") {
        filename <- paste0("01_obj_crossval_", vst, I_NI, "1K_9individuals_high.rds")
        dd_dils_FSIGEN <- dd_dils_FSIGEN_h
      } else if (dil_range == "all") {
        filename <- paste0("01_obj_crossval_", vst, I_NI, "1K_9individuals_low_high_all.rds")
        dd_dils_FSIGEN <- dd_dils_FSIGEN_lh
        cv_fun <- cross_val_create_custom_all
        # NOTE: the difference between cross_val_create_custom() and cross_val_create_custom_all()
        # is that the latter contains ALL probability estimates for the test data, not just
        # the probability estimates for the points with estimated values lower than the highest WC.
        # This makes the saved object MUCH larger in size, so when possible, it's preferable to use cross_val_create_custom().
      }
      # filename <- paste0("01_obj_crossval_", vst, I_NI, "1K_9individuals_extreme.rds")
      gc()
      
      
      dd_val_dils <- future_lapply(1:reps, function(x) {
        lapply(dd_dils_FSIGEN, function(d) {
          cv_fun(dd_GT_both,
                 dils_train = d$dil_train, dils_test = d$dil_test,
                 f=f_vst, intercept=INT, b_int=initial_beta, train_prop=train_proportion,
                 method="Nelder-Mead", hessian=F, control_list=list(maxit = 500))
        })
      }, future.seed = 1913, future.scheduling = 1)
      
      dd_val_dils <- sapply(1:length(dd_dils_FSIGEN), function(i) {
        lapply(dd_val_dils, function(x) x[[i]])
      }, simplify = FALSE)
      
      names(dd_val_dils) <- sapply(dd_dils_FSIGEN, function(d) {
        paste0("Train = ", paste0(d$dil_train, collapse = ", "), ", ", "Test = ", paste0(d$dil_test, collapse = ", "))
      })
      
      saveRDS(dd_val_dils, file.path("..", "data", filename))
    }
    
  }
}