#-----------------------#
# Global definitions ####
#-----------------------#
t_maf <- list("het" = c(0.35, 0.65), "hom" = 0.95)
t_ppc <- c(0.3, 0.7)
t_cov <- 100 # Threshold for the minimum coverage
t_q <- 0.85 # Probability threshold for symmetry model
no_call <- "NC" # Naming of no-calls (WARNING: do NOT change!!! I have used this
# object throughout the script for various logical operations
f_vst <- function(x) sqrt(x)


#------------------#
# Fitting model ####
#------------------#
# Fitting function 1:
fit_symmetry_model <- function(df, f = f_vst, intercept=TRUE,
                               b_int=c(0,1,-1),
                               method = c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN", "Brent"),
                               hessian = TRUE) {
  t1 <- f(df$A1.Reads)
  t2 <- f(df$A2.Reads)
  hom1 <- df$Genotype_true_AA == "A1A1"
  hom2 <- df$Genotype_true_AA == "A2A2"
  # sufficient statistics:
  ss0 <- sum(hom1+hom2)
  ss1 <- sum(t1*hom1+t2*hom2)
  ss2 <- sum(t2*hom1+t1*hom2)
  # minus log likelihood
  if (intercept) {
    if (length(b_int) != 3) {stop(paste0("The initial guess for the parameters has length ", length(b_int), ", but for a model with intercept, it should be a vector of length 3."))}
    ll <- function(b){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      sum(log(1+exp(b0+b1*t1+b2*t2)+exp(b0+b1*t2+b2*t1))) - b0*ss0 - b1*ss1 - b2*ss2
    }
  } else {
    if (length(b_int) != 2) {stop(paste0("The initial guess for the parameters has length ", length(b_int), ", but for a model without intercept, it should be a vector of length 2."))}
    ll <- function(b){
      b1 <- b[1]
      b2 <- b[2]
      sum(log(1+exp(b1*t1+b2*t2)+exp(b1*t2+b2*t1))) - b1*ss1 - b2*ss2
    }
  }
  
  # Fit and output
  fit <- NULL
  if (intercept) {
    par_vec <- c(beta0 = NA, beta1 = NA, beta2 = NA)
    ll_save <- function(b) {
      par_vec[1] <<- b[1]
      par_vec[2] <<- b[2]
      par_vec[3] <<- b[3]
      return(ll(b))
    }
    tryCatch({
      fit <- optim(b_int, fn = ll_save, method = method[1], hessian = hessian)
    }, error = function(e) {
      fit <- list(par = par_vec, convergence = 1)
    })
    # safe_optim <- purrr::safely(optim, otherwise = list("par"=par_vec, "convergence"=1))
    # fit <- safe_optim(b_int, fn = ll_save, method = method[1], hessian = hessian)
    if (fit$convergence == 0) {
      fisher_info <- solve(fit$hessian)
      se <- sqrt(diag(fisher_info))
      CI <- fit$par + 1.96*se %*% t(c(-1,1))
    } else {
      CI <- matrix(NA, nrow = 3, ncol = 2)
    }
    colnames(CI) <- c("Lower", "Upper")
    rownames(CI) <- c("beta0", "beta1", "beta2")
  } else {
    par_vec <- c(beta1 = NA, beta2 = NA)
    ll_save <- function(b) {
      par_vec[1] <<- b[1]
      par_vec[2] <<- b[2]
      return(ll(b))
    }
    tryCatch({
      fit <- optim(b_int, fn = ll_save, method = method[1], hessian = hessian)
    }, error = function(e) {
      fit <- list(par = par_vec, convergence = 1)
    })
    # safe_optim <- purrr::safely(optim, otherwise = list("par"=par_vec, "convergence"=1))
    # fit <- safe_optim(b_int[-1], fn = ll_save, method = method[1], hessian = hessian)
    if (fit$convergence == 0) {
      fisher_info <- solve(fit$hessian)
      se <- sqrt(diag(fisher_info))
      CI <- fit$par + 1.96*se %*% t(c(-1,1))
    } else {
      CI <- matrix(NA, nrow = 2, ncol = 2)
    }
    colnames(CI) <- c("Lower", "Upper")
    rownames(CI) <- c("beta1", "beta2")
  }
  
  return(list("par"=fit$par, "CI"=CI, "f"=f, "Converged" = fit$convergence))
}

# Fitting function 1:
fit_symmetry_model2 <- function(df, f = f_vst, intercept=TRUE, b_int=c(0,1,-2),
                                method="Nelder-Mead", control_list=list(maxit = 500),
                                hessian=FALSE){
  t1 <- f(df$A1.Reads)
  t2 <- f(df$A2.Reads)
  hom1 <- df$Genotype_true_AA == "A1A1"
  hom2 <- df$Genotype_true_AA == "A2A2"
  # sufficient statistics:
  ss0 <- sum(hom1+hom2)
  ss1 <- sum(t1*hom1+t2*hom2)
  ss2 <- sum(t2*hom1+t1*hom2)
  # minus log likelihood
  if (intercept) {
    if (length(b_int) != 3) {stop(paste0("The initial guess for the parameters has length ", length(b_int), ", but for a model with intercept, it should be a vector of length 3."))}
    ll <- function(b){
      b0 <- b[1]
      b1 <- b[2]
      b2 <- b[3]
      sum(log(1+exp(b0+b1*t1+b2*t2)+exp(b0+b1*t2+b2*t1))) - b0*ss0 - b1*ss1 - b2*ss2
    }
  } else {
    if (length(b_int) != 2) {stop(paste0("The initial guess for the parameters has length ", length(b_int), ", but for a model without intercept, it should be a vector of length 2."))}
    ll <- function(b){
      b1 <- b[1]
      b2 <- b[2]
      sum(log(1+exp(b1*t1+b2*t2)+exp(b1*t2+b2*t1))) - b1*ss1 - b2*ss2
    }
  }
  
  # Fit and output
  if (intercept) {
    fit <- optim(b_int, fn = ll, method = method, control=control_list, hessian=hessian)
  } else {
    fit <- optim(b_int, fn = ll, method = method, control=control_list, hessian=hessian)
  }
  
  if (hessian) {
    fisher_info <- solve(fit$hessian)
    se <- sqrt(diag(fisher_info))
    CI <- fit$par + 1.96*se %*% t(c(-1,1))
    return(list("par"=fit$par, "Converged" = fit$convergence, "f"=f, "Value" = fit$value, "CI"=CI))
  } else {
    return(list("par"=fit$par, "Converged" = fit$convergence, "f"=f, "Value" = fit$value))
  }
}


#--------------------#
# SMLR prediction ####
#--------------------#
predict_prob_threshold <- function(test_data, fit, threshold = FALSE, GT_nucleotide = TRUE) {
  t1 <- fit$f(test_data$A1.Reads)
  t2 <- fit$f(test_data$A2.Reads)
  if (length(fit$par) == 2) {
    b <- c(0, fit$par)
  } else {
    b <- fit$par
  }
  A1A1 <- exp(b[1] + b[2]*t1 + b[3]*t2) / (exp(b[1] + b[2]*t1 + b[3]*t2) + exp(b[1] + b[3]*t1 + b[2]*t2) + 1)
  A2A2 <- exp(b[1] + b[3]*t1 + b[2]*t2) / (exp(b[1] + b[2]*t1 + b[3]*t2) + exp(b[1] + b[3]*t1 + b[2]*t2) + 1)
  A1A2 <- 1 / (exp(b[1] + b[2]*t1 + b[3]*t2) + exp(b[1] + b[3]*t1 + b[2]*t2) + 1)
  
  # In some cases, the MLEs and transformed signals can make the both the
  # numerator and denominator result in Inf/Inf, leading to an NaN for one of
  # the three estimated probabilities. If this happens for the p_max, we get a
  # problem when predicting, but we can mitigate the problem by inferring the
  # value of Inf/Inf from the other two probabilities.
  nan11 <- which(is.nan(A1A1))
  nan22 <- which(is.nan(A2A2))
  nan12 <- which(is.nan(A1A2))
  if (length(nan11) > 0) {A1A1[nan11] <- 1 - A2A2[nan11] - A1A2[nan11]}
  if (length(nan22) > 0) {A2A2[nan22] <- 1 - A1A1[nan22] - A1A2[nan22]}
  if (length(nan12) > 0) {A1A2[nan12] <- 1 - A1A1[nan12] - A2A2[nan12]}
  
  pred_prob <- cbind(A1A1, A1A2, A2A2)
  p_max <- apply(pred_prob,
                 MARGIN = 1,
                 FUN = max) %>% unname
  Genotype_AA <- apply(pred_prob,
                       MARGIN = 1,
                       FUN = function(x) {names(which(x==max(x)))}) %>% unname
  
  # NOTE: I have experienced that Genotype_AA becomes a list if two (or more)
  # classes share the same probability that happens to equal p_max.
  # In such cases we currently decide to make no prediction:
  NN <- apply(cbind(p_max, pred_prob), MARGIN = 1, FUN = function(x) {sum(x[1] == x[-1]) > 1})
  Genotype_AA[NN] <- no_call
  # If we observe 0 reads in both channels we also choose to make no prediction:
  # (this can now be omitted as we decided to filter away double zeros for all
  # of our data in the beginning of the analysis, since all models and methods
  # should automatically make a no call or NA when both signals are zero)
  NN <- (test_data$A1.Reads == 0 & test_data$A2.Reads == 0)
  Genotype_AA[NN] <- no_call
  # Hopefully the overwriting with no_call means that we can unlist Genotype_AA
  # without getting a vector that is longer than the number of observations:
  if (length(Genotype_AA) == length(unlist(Genotype_AA))) {
    Genotype_AA <- unlist(Genotype_AA)
  } else {
    stop("Genotype_AA seems to contain more elements than the number of observations in the data frame.")
  }
  
  if (GT_nucleotide) {
    Genotype_predict <- if_else(Genotype_AA == "A1A1",
                                paste0(test_data$A1, test_data$A1),
                                if_else(Genotype_AA == "A1A2",
                                        paste0(test_data$A1, test_data$A2),
                                        if_else(Genotype_AA == "A2A2",
                                                paste0(test_data$A2, test_data$A2), no_call)))
  }
  
  if (threshold) {
    Genotype_AA[p_max < threshold] <- no_call
    if (GT_nucleotide) {Genotype_predict[p_max < threshold] <- no_call}
  }
  
  if (GT_nucleotide) {test_data <- test_data %>% mutate(Genotype_pred = Genotype_predict)}
  
  colnames(pred_prob) <- paste0("p_", colnames(pred_prob))
  
  test_data <- test_data %>% 
    mutate(Genotype_pred_AA = Genotype_AA,
           p_max = p_max) %>% 
    bind_cols(as.data.frame(pred_prob))
  
  return(test_data)
}


#---------------------#
# HSGQC prediction ####
#---------------------#
threshold_rule <- function(test_data, qc = list("MAF" = t_maf, "PPC" = t_ppc), gq = FALSE, gt_AA = FALSE){
  
  Genotype_predict <- test_data$Genotype
  # Genotype_predict <- if_else(is.na(test_data$QC), Genotype_predict, no_call)
  
  if (is.list(qc)) {
    if (is.list(qc$MAF)) {
      Genotype_predict <- case_when(Genotype_predict %in% c(no_call, "-") ~ no_call,
                                    Genotype_predict %in% c("AA", "CC", "GG", "TT") & (test_data$Maj.Allele.Freq < qc$MAF$hom*100) ~ no_call,
                                    !(Genotype_predict %in% c("AA", "CC", "GG", "TT")) & (test_data$Maj.Allele.Freq < qc$MAF$het[1]*100 | test_data$Maj.Allele.Freq > qc$MAF$het[2]*100) ~ no_call,
                                    .default = Genotype_predict)
    }
    if (is.numeric(qc$PPC)) {
      Genotype_predict[test_data$Perc.Pos.Cov < qc$PPC[1]*100 | test_data$Perc.Pos.Cov > qc$PPC[2]*100] <- no_call
    }
  }
  if (is.numeric(gq)) {
    Genotype_predict <- if_else(test_data$GQ < gq, no_call, Genotype_predict)
  }
  
  if (gt_AA) {
    Genotype_AA <- case_when(Genotype_predict %in% c(no_call, "-") ~ no_call,
                             Genotype_predict == paste0(test_data$A1, test_data$A1) ~ "A1A1",
                             Genotype_predict == paste0(test_data$A2, test_data$A2) ~ "A2A2",
                             .default = "A1A2")
    return(list("Predicted" = Genotype_predict, "AA" = Genotype_AA))
  } else {
    return(Genotype_predict)
  }
}


#------------------------------------------#
# Model from Mostad, Tillmar, and Kling ####
#------------------------------------------#
predict_MTK_model <- function(df, m, e=0.003, g_prior=c("uniform_bi", "uniform_tetra"), negloglik_me=F) {
  if (is.character(g_prior)) {
    if (g_prior[1] == "uniform_bi") {
      g_prior <- df |> select(Target.ID, A1, A2) |> group_by(Target.ID) |> 
        summarise(pA = if_else("A" %in% c(A1, A2), 0.5, 0),
                  pC = if_else("C" %in% c(A1, A2), 0.5, 0),
                  pG = if_else("G" %in% c(A1, A2), 0.5, 0),
                  pT = if_else("T" %in% c(A1, A2), 0.5, 0)) |> 
        mutate(AA = pA^2, CC = pC^2, GG = pG^2, TT = pT^2,
               AC = 2*pA*pC, AG = 2*pA*pG, AT = 2*pA*pT,
               CG = 2*pC*pG, CT = 2*pC*pT,
               GT = 2*pG*pT) |> select(!c(pA, pC, pG, pT))
    } else if (g_prior[1] == "uniform_tetra") {
      g_prior <- df |> select(Target.ID) |> group_by(Target.ID) |> 
        summarise(pA = 0.25, pC = 0.25, pG = 0.25, pT = 0.25) |> 
        mutate(AA = pA^2, CC = pC^2, GG = pG^2, TT = pT^2,
               AC = 2*pA*pC, AG = 2*pA*pG, AT = 2*pA*pT,
               CG = 2*pC*pG, CT = 2*pC*pT,
               GT = 2*pG*pT) |> select(!c(pA, pC, pG, pT))
    }}
  Pg <- df |> select(Target.ID) |> left_join(g_prior, by = "Target.ID") |> 
    mutate(across(!c("Target.ID", "AA", "CC", "GG", "TT"), ~ .x/2)) |> 
    mutate(CA = AC, GA = AG, TA = AT, GC = CG, TC = CT, TG = GT)
  gt <- Pg |> select(!Target.ID) |> names()
  
  k <- 0:m
  q <- k/m
  
  c1 <- df$A.Reads
  c2 <- df$C.Reads
  c3 <- df$G.Reads
  c4 <- df$T.Reads
  cm <- unname(cbind(c1,c2,c3,c4)) |> apply(MARGIN=1, function(x){sort(x, decreasing = F)}) |> t()
  
  # Note that lchoose(cm[,1], cm[,1]) always equals 0, which we exploit to get
  # better numerical stability: the above sorting makes sure that cm[,1] is the
  # largest count
  m_coef <- lchoose(cm[,1] + cm[,2],                   cm[,2]) + # The difficult coefficient (large factorial).
            lchoose(cm[,1] + cm[,2] + cm[,3],          cm[,3]) + # Since cm[,3] and cm[,4] are small numbers,
            lchoose(cm[,1] + cm[,2] + cm[,3] + cm[,4], cm[,4])   # the last two coefficients are easy to compute.
  
  lb <- function(alleles, acgt) {log((1-e)*(q*(alleles[1]==acgt) + (1-q)*(alleles[2]==acgt)) + e/4)}
  Pcg <- function(a) {apply(exp(t(m_coef +
                                  outer(c1, lb(a,"A")) +
                                  outer(c2, lb(a,"C")) +
                                  outer(c3, lb(a,"G")) +
                                  outer(c4, lb(a,"T"))) +
                                  lchoose(m, k) - m*log(2)),
                            MARGIN = 2, sum)}
  
  Pc <- lapply(gt, function(g) {
    a <- str_split_fixed(g, "", 2) |> c()
    Pcg(a)*Pg[,g]
  }) |> bind_cols() |> apply(MARGIN = 1, sum)
  
  if (!negloglik_me) {
    Pgc <- lapply(gt, function(g) {Pcg(c(str_split_fixed(g, "", 2)))*Pg[,g]/Pc}) |>
      bind_cols() |>
      mutate(AC = AC+CA, AG = AG+GA, AT = AT+TA,
             CG = CG+GC, CT = CT+TC,
             GT = GT+TG) |> select(AA, CC, GG, TT, AC, AG, AT, CG, CT, GT)
    
    MTK_pmax <- apply(Pgc, MARGIN = 1, max)
    MTK_pred <- apply(Pgc, MARGIN = 1, function(x) {names(Pgc)[which.max(x)]})
    return(df |> mutate("MTK_pred" = MTK_pred, "MTK_pmax" = MTK_pmax))
  } else {
    return(sum(-log(Pc)))
  }
}


#-----------------#
# Plot Machine ####
#-----------------#
plot_machine <- function(df, f=f_vst, legend=TRUE, title=NULL, x_lab=F, y_lab=F, x_lim=F, y_lim=F, complete_theme=theme_bw()){
  dd_missed <- df %>%
    filter(Genotype != Genotype_true) %>%
    mutate(Missed = if_else(Genotype == no_call, no_call, "WC")) |> 
    arrange(Missed)

  machineplot <- ggplot(df, aes(x = f(A1.Reads), y = f(A2.Reads))) +
    complete_theme +
    geom_point(aes(colour = Genotype_true_AA), size = 1, alpha = 0.15) +
    geom_point(aes(colour = Genotype_true_AA), data = dd_missed, size = 1) +
    scale_colour_manual(name="True genotype",
                        values=c("#F0E442", "#CC79A7", "#56B4E9"),
                        breaks=c("A1A1", "A1A2", "A2A2"),
                        labels=unname(c(TeX("$a_1 a_1$"), TeX("$a_1 a_2$"), TeX("$a_2 a_2$"))))
  if (is.character(x_lab) | is.expression(x_lab)) {
    machineplot <- machineplot + xlab(x_lab)
  } else {
    machineplot <- machineplot + xlab(TeX("$t_1$"))
  }
  if (is.character(y_lab) | is.expression(y_lab)) {
    machineplot <- machineplot + ylab(y_lab)
  } else {
    machineplot <- machineplot + ylab(TeX("$t_2$"))
  }
  if (is.numeric(x_lim) & is.numeric(y_lim)) {
    machineplot <- machineplot + coord_cartesian(xlim = x_lim, ylim = y_lim)
  } else if (is.numeric(x_lim)) {
    machineplot <- machineplot + coord_cartesian(xlim = x_lim, ylim = c(0, f(max(df$A2.Reads))))
  } else if (is.numeric(y_lim)) {
    machineplot <- machineplot + coord_cartesian(xlim = c(0, f(max(df$A1.Reads))), ylim = y_lim)
  }
  machineplot <- machineplot +
    # geoms below will use another colour scale
    new_scale_color() +
    geom_point(aes(shape = Missed, colour = Missed), data = dd_missed) +
    scale_shape_manual(values=c(3, 4)) +
    scale_colour_manual(name="Missed",
                        values=c("black", "red"),
                        breaks=c(no_call, "WC"),
                        labels=unname(c(no_call, "WC")))
  
  if (!legend) {machineplot <- machineplot + theme(legend.position="none")}
  if (is.character(title)) {machineplot <- machineplot + labs(title = title)}
  
  return(machineplot)
}


#---------------------------#
# Plot EQC (HSG with QC) ####
#---------------------------#
plot_EQC <- function(df, f=f_vst, threshold_maf = t_maf, legend=TRUE, title=NULL, complete_theme=theme_bw()){
  # The creation of the two lines below may look a bit weird, but we want the
  # Method-variable for the purpose of putting the lines in a legend.
  # We also want the lines to be as long as possible, but without extending them
  # beyond the maximal points.
  s1 <- c(0, max(df$A1.Reads))
  s2 <- c(0, max(df$A2.Reads))
  dd_poly1 <- data.frame(t1 = f(c(s1, rev(s1))),
                         t2 = f(c(s1*((1-threshold_maf$hom)/threshold_maf$hom), rev(s1)*(threshold_maf$het[1]/(1-threshold_maf$het[1])))),
                         Method = "HSGQC")
  dd_poly2 <- data.frame(t1 = f(c(s2*((1-threshold_maf$hom)/threshold_maf$hom), rev(s2)*(threshold_maf$het[1]/(1-threshold_maf$het[1])))),
                         t2 = f(c(s2, rev(s2))),
                         Method = "HSGQC")
  
  dd_missed <- df %>% 
    mutate(Genotype_pred = threshold_rule(df)) %>%
    filter(Genotype_pred != Genotype_true) %>%
    mutate(Missed = if_else(Genotype_pred == no_call, no_call, "WC")) |> 
    arrange(Missed)
  
  EQCplot <- ggplot(df, aes(x = f(A1.Reads), y = f(A2.Reads))) + complete_theme
  if (is.list(threshold_maf)) {
    EQCplot <- EQCplot +
      geom_polygon(aes(x=t1, y=t2, fill=Method), alpha=0.75, data = dd_poly1) +
      geom_polygon(aes(x=t1, y=t2, fill=Method), alpha=0.75, data = dd_poly2) +
      scale_fill_manual(name = "Allele balance",
                        values = "gray75",
                        labels = "Hom < 0.95 or\nHet < 0.35 or\nHet > 0.65")
  }
  
  EQCplot <- EQCplot +
    geom_point(aes(colour = Genotype_true_AA), size = 1, alpha = 0.15) +
    geom_point(aes(colour = Genotype_true_AA), data = dd_missed, size = 1) +
    scale_colour_manual(name="True genotype",
                        values=c("#F0E442", "#CC79A7", "#56B4E9"),
                        breaks=c("A1A1", "A1A2", "A2A2"),
                        labels=unname(c(TeX("$a_1 a_1$"), TeX("$a_1 a_2$"), TeX("$a_2 a_2$")))) +
    ylim(0, f(max(df$A2.Reads))) +
    xlab(TeX("$t_1$")) +
    ylab(TeX("$t_2$")) +
    # geoms below will use another colour scale
    new_scale_color() +
    geom_point(aes(shape = Missed, colour = Missed), data = dd_missed) +
    scale_shape_manual(values=c(3, 4)) +
    scale_colour_manual(name="Missed",
                        values=c("black", "red"),
                        breaks=c(no_call, "WC"),
                        labels=unname(c(no_call, "WC")))
  
  if (!legend) {EQCplot <- EQCplot + theme(legend.position="none")}
  if (is.character(title)) {EQCplot <- EQCplot + labs(title = title)}
  
  return(EQCplot)
}


#--------------#
# Plot SMLR ####
#--------------#
# If s1 and s2 are the raw signals from the possible base channels and f is a
# variance stabilizing transformation, then we define t1=f(s1) and t2=f(s2)
band_line <- function(df, model, q_threshold, f=f_vst){
  P11q <- function(t1, t2, q=q_threshold, b=model$par){
    if (length(b)==3) {
      exp(b[1] + b[2]*t1 + b[3]*t2) / (exp(b[1] + b[2]*t1 + b[3]*t2) + exp(b[1] + b[3]*t1 + b[2]*t2) + 1) - q
    } else {
      exp(b[1]*t1 + b[2]*t2) / (exp(b[1]*t1 + b[2]*t2) + exp(b[2]*t1 + b[1]*t2) + 1) - q
    }
  }
  P12q <- function(t1, t2, q=q_threshold, b=model$par){
    if (length(b)==3) {
      1 / (exp(b[1] + b[2]*t1 + b[3]*t2) + exp(b[1] + b[3]*t1 + b[2]*t2) + 1) - q
    } else {
      1 / (exp(b[1]*t1 + b[2]*t2) + exp(b[2]*t1 + b[1]*t2) + 1) - q
    }
  }
  # There will be a probability band around each separation line, but depending
  # on the value of q and whether the separation lines intersect each other in
  # the 1st quadrant, the probability band will either be one connected area or
  # two areas (when the inner limits of the band don't extend over the origin).
  # The inner limits of the probability band are where P(het|t1,t2) = q, i.e.
  # where P12q(t1,t2) = 0, and the outer limits are where P11q(t1,t2) = 0.
  # We want to search for the points (t1, t2) that meet these conditions.
  
  # We don't want to plot the probability bands further out than our data, so
  # the maximal points we will search for is simply determined by our data:
  t1_max <- max(f(df$A1.Reads))
  t2_max <- max(f(df$A2.Reads))
  t_max <- max(t1_max, t2_max)
  t_2nd.max <- min(t1_max, t2_max)
  # NOTE: These points will lie on the limit of the probability bands.
  
  # For one connected band area, we will consider the following:
  # We have that P12q(t1,t2) = P(het|t1,t2) - q, which is maximal when t1=t2.
  # Hence, if P(het|0,0) < q then P12q(0,0) < 0 and there will be a minimal
  # value of t1=t2 where P(het|t1=t2) = q which will mark the point where the
  # inner limits meet, i.e. a point where we don't need to search behind.
  # For two band areas, we shouldn't search behind the origin.
  
  # TODO: For all configurations of the separation lines, we can always increase
  # the q-threshold until P12q(0,0) < 0 becomes true. But for lower q-values, I
  # am not sure about the logical value of P12q(0,0) < 0 when the two separation
  # lines intersect in the 1st quadrant. That configuration seems difficult to
  # handle, so for article plots, one may then have to choose a large enough q.
  if (P12q(0,0) < 0) {
    min_P12_tx <- uniroot(function(x) P12q(x, x), lower=0, upper=t_max, maxiter=10^5, tol=10^(-9))$root
  } else {
    min_P12_tx <- uniroot(function(x) P12q(x, 0), lower=0, upper=t_max, maxiter=10^5, tol=10^(-9))$root
    if (min_P12_tx <= 0) {
      # When the two separation lines intersect in the 1st quadrant:
      min_P12_tx <- 0
    }
  }
  # Using if_else gives an error
  # NOTE: For one band area, this point will also lie on the limit, but
  #       for two band areas it will NOT lie on the limit!!!
  
  # Now imagine that we plot t2 against t1 (i.e. x=t1 and y=t2).
  # We define a sequence of t1-points for which we will find the t2-points that
  # lies on the inner limits, i.e. satisfies P12q(t1,t2) = 0.
  # The band limits will almost be straight lines for high t1- and t2-values,
  # but not necessarily for the low values, so we want a sequence where the
  # t1-points are closer in the lower end.
  # As explained later, it makes sense to start searching from the higher values
  P12_tx <- c(t_max,
              seq(from=t_2nd.max, to=min_P12_tx+0.55, length.out=10),
              seq(from=min_P12_tx+0.5,
                  to=min_P12_tx+0.15, by=-0.05),
              seq(from=min_P12_tx+0.1, to=min_P12_tx, by=-0.005))
  
  # We create a vector to fill with the t2-values that solves P12q(t1,t2) = 0.
  # Note that if P12q(0,0) < 0 then we already know the last of these values
  # which is min_P12_tx, which is sometimes more precisely determined when
  # searching along P12q(x,x) rather than alogn P12q(min_P12_tx,y), so we insert
  # this in an if-statement later.
  P12_ty <- rep(0, length(P12_tx))
  
  # We use uniroot in our search and this function needs to be told a lower and
  # an upper end for the search interval AND the end points should have opposite
  # signs when evaluated in the target function.
  # For this problem, zero will always be the lower limit, but the upper limit
  # depends on the data and on the function used to transform the data.
  # We already decided not to extend the probability bands further out than our
  # data, which gives us the upper end for the search interval.
  # We will always be able to find the t2-value satisfying P12q(t_max,t2) = 0,
  # but this is not always easy for P12q(min_t1t2, t2) = 0 (if two band areas).
  # A way to create a general code that can handle both one and two band areas
  # is to start our search for t2-values at t1=t_max and then move downwards
  # until t2=min_t1t2 (for two band areas, this happens before t1=min_t1t2).
  
  # In our search for the t2 in P12q(t_max,t2) = 0, we will search in [0; t_max]
  P12_ty[1] <- uniroot(function(y) P12q(t_max, y), lower=0, upper=t_max, maxiter = 10^5, tol = 10^(-9))$root
  # But now the upper search end should be chosen carefully: there are two
  # separation lines, hence in principle always two bands and inner limits,
  # i.e. there will be two solutions to the equation P12q(t1,t2) = 0.
  # If we plot something like
  # plot(seq(0,12,0.1), P12q(P12_tx[1], seq(0,12,0.1))); abline(h=0)
  # we will see a unimodal hill that intersects the x-axis (i.e. y=0) in two
  # places, so in order for the search interval's end points to have opposite
  # signs, our search interval shouldn't go from the one side (x=0) of the hill
  # to the other, but from x=0 to somewhere on the top of the hill.
  # If we draw a line perpendicular to the identity line in the (t1 x t2)-plane,
  # then the hill top of the profile along this line will be at t1=t2, but this
  # is not true for the vertical line at the t1-value we search at, where the
  # hill top will be a (very) little above t1=t2.
  # Therefore, the t1-value we search at is a good choice for the upper end of
  # the search interval.
  # Remember that we already know the last value of P12_ty when P12q(0,0) < 0
  if (P12q(0,0) < 0) {
    P12_ty[length(P12_tx)] <- min_P12_tx
    for (i in 2:(length(P12_tx)-1)) {
      P12_ty[i] <- uniroot(function(y) P12q(P12_tx[i], y), lower=-1, upper=P12_tx[i], maxiter = 10^5, tol = 10^(-9))$root
      # P12_ty[i] <- uniroot(function(y) P12q(P12_tx[i], y), lower=0, upper=P12_ty[i-1]+0.01, maxiter = 10^5, tol = 10^(-9))$root
    }
  } else if (min_P12_tx == 0) {
    # When the two separation lines intersect in the 1st quadrant then we may
    # have that P12q(0,0) < 0), but if not then min_P12_tx == 0 and we will have
    # to search for the corresponding P12_ty(0,t2):
    for (i in 2:length(P12_tx)) {
      P12_ty[i] <- uniroot(function(y) P12q(P12_tx[i], y), lower=-1, upper=P12_tx[i], maxiter = 10^5, tol = 10^(-9))$root
    }
  } else {
    # When P12q(0,0) >= 0 and the two separation lines don't intersect in the
    # 1st quadrant, then we know that the 'inner limit' of the lower band will
    # intersect the x-axis, so the last value of P12_ty should be a clear zero
    P12_ty[length(P12_tx)] <- 0
    for (i in 2:(length(P12_tx)-1)) {
      P12_ty[i] <- uniroot(function(y) P12q(P12_tx[i], y), lower=-1, upper=P12_tx[i], maxiter = 10^5, tol = 10^(-9))$root
    }
  }
  
  
  # We follow the same principles for the lower limits of the probability bands.
  # We make a t1-sequence and for each of its elements we want to find the t2
  # that solves P11q(t1,t2) = 0, but this equation only has one solution.
  
  # Consider the lower separation line: it can intersect the positive t1-axis or
  # the positive t2-axis. In the former case the lower end in our sequence of
  # t1-points should be zero, and in the latter case it should be some positive
  # number. Since the former case is a possibility, we should set the lower end
  # of the search interval to some "large enough" negative number when we search
  # for the lowest t1-point in our t1-sequence.
  min_P11_tx <- uniroot(function(x) P11q(x, 0), lower = -t_max/2, upper = t_max, maxiter = 10^5, tol = 10^(-9))$root
  if (min_P11_tx <= 0) {min_P11_tx <- 0}
  P11_tx <- c(t_max, seq(t_2nd.max, min_P11_tx, length.out=10))
  P11_ty <- rep(0, length(P11_tx))
  P11_ty[1] <- uniroot(function(y) P11q(t_max, y), lower=0, upper=t_max, maxiter = 10^5, tol = 10^(-9))$root
  for (i in 2:(length(P11_tx)-1)) {
    # OBS: we write 2:(length(P11_tx)-1)) since we want P11_ty[length(P11_tx)]
    # to be a clear zero when min_P11_tx > 0
    P11_ty[i] <- uniroot(function(y) P11q(P11_tx[i], y), lower=0, upper=P11_tx[i], maxiter = 10^5, tol = 10^(-9))$root
    # P11_ty[i] <- uniroot(function(y) P11q(P11_tx[i], y), lower=0, upper=P11_ty[i-1], maxiter = 10^5, tol = 10^(-9))$root
  }
  if (min_P11_tx == 0) {
    P11_ty[length(P11_tx)] <- uniroot(function(y) P11q(0, y), lower=0, upper=t_max, maxiter = 10^5, tol = 10^(-9))$root
  } else {
    # If min_P11_tx > 0 then P11_ty[length(P11_tx)] should be 0 as already defined
    P11_tx <- c(P11_tx, 0)
    P11_ty <- c(P11_ty, 0)
    P12_tx <- c(P12_tx, 0)
    P12_ty <- c(P12_ty, 0)
  }
  
  
  # if (min_P11_tx >= 0) {
  #   P11_tx <- c(t_max, seq(t_2nd.max, min_P11_tx, length.out=10))
  #   P11_ty <- rep(0, length(P11_tx))
  # } else {
  #   P11_tx <- c(t_max, seq(t_2nd.max, 0, length.out=10))
  #   P11_ty <- rep(0, length(P11_tx))
  #   P11_ty[length(P11_tx)] <- uniroot(function(y) P11q(0, y), lower=0, upper=t_max, maxiter = 10^5, tol = 10^(-9))$root
  # }
  # P11_ty[1] <- uniroot(function(y) P11q(t_max, y), lower=0, upper=t_max, maxiter = 10^5, tol = 10^(-9))$root
  # for (i in 2:(length(P11_tx)-1)) {
  #   P11_ty[i] <- uniroot(function(y) P11q(P11_tx[i], y), lower=0, upper=P11_tx[i], maxiter = 10^5, tol = 10^(-9))$root
  #   # P11_ty[i] <- uniroot(function(y) P11q(P11_tx[i], y), lower=0, upper=P11_ty[i-1], maxiter = 10^5, tol = 10^(-9))$root
  # }
  # P11_tx <- c(P11_tx, 0)
  # P11_ty <- c(P11_ty, 0)
  
  # TODO: These polygons probably look ugly when min_P11_tx == 0, but maybe not, if q is sufficiently high
  df_polygon1 <- tibble("x" = c(rev(P11_tx), P12_tx),
                        "y" = c(rev(P11_ty), P12_ty)) %>%
    filter(x >= 0, y >= 0, x <= t1_max)
  df_polygon2 <- tibble("x" = c(rev(P12_ty), P11_ty),
                        "y" = c(rev(P12_tx), P11_tx)) %>%
    filter(x >= 0, y >= 0, y <= t2_max)
  
  return(bind_rows(df_polygon1, df_polygon2))
}


# Vi skal plotte separationslinjen for den model, vi vælger til at prædiktere
# genotyper med, men det kan diskuteres, om vi skal plotte en linje for begge
# modeller (m/u intercept).
plot_model <- function(intercept_model=NULL, NI_model=NULL, df, q = t_q, predictions = c("Intercept", "No intercept"), q_digits = 2, legend=TRUE, title=NULL, complete_theme=theme_bw()){
  if (is.null(intercept_model) & is.null(NI_model)) {
    warning("You need to supply a model for at least one of the model arguments (intercept_model or NI_model)")
  } else if (is.null(intercept_model)) {
    predictions <- "No intercept"
  } else if (is.null(NI_model)) {
    predictions <- "Intercept"
  }
  if (predictions[1] == "Intercept") {
    pred_model <- intercept_model
  } else if (predictions[1] == "No intercept") {
    pred_model <- NI_model
  }
  f <- pred_model$f
  
  if (!is.null(NI_model)) {
    b_NI <- NI_model$par
    l2_t1.max_NI <- uniroot(function(x) -b_NI[2]*x/b_NI[1]-f(max(df$A2.Reads)), interval = c(0, f(max(df$A2.Reads))))$root
    dd_lines1_NI <- data.frame(t1 = seq(from=0, to=f(max(df$A1.Reads)), length.out=63)) %>%
      arrange(t1) %>%
      filter(t1 >= 0) %>%
      mutate(t2 = -b_NI[1]*t1/b_NI[2], Model = "no_intercept")
    dd_lines1 <- dd_lines1_NI # to be overwritten if both intercept_model and NI_model are supplied
  }
  if (!is.null(intercept_model)) {
    b <- intercept_model$par
    l2_t1.max_I <- uniroot(function(x) (-b[1]-b[3]*x)/b[2]-f(max(df$A2.Reads)), interval = c(0, f(max(df$A2.Reads))))$root
    # General line definitions that handle both:
    # Two full lines (the one will always intersect the x- and the other the y-axis).
    # **Two lines that fuse into the identity line if they intersect in the 1st quadrant.
    l1_root <- uniroot(function(x) (-b[1]-b[2]*x)/b[3], interval = c(-1, 1)*f(max(df$A2.Reads)))$root
    l2_root <- uniroot(function(x) (-b[1]-b[3]*x)/b[2], interval = c(-1, 1)*f(max(df$A2.Reads)))$root
    l1_l2_meet <- uniroot(function(x) (-b[1]-b[3]*x)/b[2] - (-b[1]-b[2]*x)/b[3], interval = c(-5, 5)*f(max(df$A2.Reads)))$root
    dd_lines1 <- data.frame(t1 = c(l1_root+10^(-6), l1_l2_meet+10^(-6), seq(from=0, to=f(max(df$A1.Reads)), length.out=63))) %>%
      arrange(t1) %>%
      filter(t1 >= 0) %>%
      mutate(t2 = (-b[1]-b[2]*t1)/b[3], Model = "intercept")
  }
  if ((!is.null(intercept_model)) & (!is.null(NI_model))) {
    l2_t1.max <- max(l2_t1.max_I, l2_t1.max_NI)
    dd_lines1 <- dd_lines1 |> bind_rows(dd_lines1_NI)
    dd_lines2 <- data.frame(t1 = c(l2_root+10^(-6), l1_l2_meet+10^(-6), seq(from=0, to=l2_t1.max, length.out=63))) %>%
      arrange(t1) %>%
      filter(t1 >= 0) %>%
      mutate(intercept = (-b[1]-b[3]*t1)/b[2],
             no_intercept = -b_NI[2]*t1/b_NI[1]) %>%
      pivot_longer(cols = c(intercept, no_intercept),
                   names_to = "Model",
                   values_to = "t2")
  } else if (!is.null(NI_model)) {
    dd_lines2 <- data.frame(t1 = seq(from=0, to=l2_t1.max_NI, length.out=63)) %>%
      arrange(t1) %>%
      filter(t1 >= 0) %>%
      mutate(t2 = -b_NI[2]*t1/b_NI[1], Model = "no_intercept")
  } else if (!is.null(intercept_model)) {
    dd_lines2 <- data.frame(t1 = c(l2_root+10^(-6), l1_l2_meet+10^(-6), seq(from=0, to=l2_t1.max_I, length.out=63))) %>%
      arrange(t1) %>%
      filter(t1 >= 0) %>%
      mutate(t2 = (-b[1]-b[3]*t1)/b[2], Model = "intercept")
  }
  # **Uncomment to get two lines that fuse into the identity line if they intersect in the 1st quadrant.
  # if (l1_l2_meet > 0) {
  #   dd_lines1$t1[dd_lines1$t1 < l1_l2_meet] <- 0
  #   dd_lines1$t2[dd_lines1$t2 < l1_l2_meet] <- 0
  #   dd_lines2$t1[dd_lines2$t1 < l1_l2_meet] <- 0
  #   dd_lines2$t2[dd_lines2$t2 < l1_l2_meet] <- 0
  # }
  dd_lines1 <- dd_lines1 %>% filter(t2 >= 0, t2 <= f(max(df$A2.Reads)))
  dd_lines2 <- dd_lines2 %>% filter(t2 >= 0, t2 <= f(max(df$A2.Reads)))
  
  
  if (q > 0) {
    dd_bands <- band_line(df, model = pred_model, q_threshold = q, f = f) %>%
      mutate(Band = paste0("q = ", round(q, digits = q_digits)))
  }
  
  dd_missed <- predict_prob_threshold(df, pred_model, q) %>%
    filter(Genotype_pred_AA != Genotype_true_AA) %>%
    mutate(Missed = if_else(Genotype_pred_AA == no_call, no_call, "WC")) |> 
    arrange(Missed)
  
  modelplot <- ggplot(df, aes(x = f(A1.Reads), y = f(A2.Reads))) + complete_theme
  if (q > 0) {
    # Should it be a rule of thumb not to use scales (legend space) for variables
    # with only one category, but leave such explanation to the caption?
    # geom_polygon(aes(x=x, y=y), fill = "grey75", data = dd_bands) +
    modelplot <- modelplot +
      # geom_polygon(aes(x=x, y=y, fill = Band), alpha = 0.75, data = dd_bands) +
      geom_polygon(aes(x=x, y=y), fill = "grey50", alpha = 0.75, data = dd_bands) #+
      # scale_fill_manual(name = "Threshold",
      #                   values = "grey50",
      #                   labels = unname(c(TeX(paste0("q = ", round(q, digits = q_digits))))))
  }
  if ((!is.null(intercept_model)) & (!is.null(NI_model))) {
    modelplot <- modelplot +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines1) +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines2) +
      scale_linetype_manual(values = c("solid", "dotted"),
                            labels = unname(c(TeX("$\\beta_0 \\neq 0$"),
                                              TeX("$\\beta_0 = 0$"))))
  } else if (!is.null(intercept_model)) {
    modelplot <- modelplot +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines1) +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines2) +
      scale_linetype_manual(values = "solid",
                            # labels = unname(TeX("$\\beta_0 \\neq 0$")))
                            labels = unname(TeX(paste0("$\\beta_0 \\neq 0, q = $", round(q, digits = q_digits)))))
  } else if (!is.null(NI_model)) {
    modelplot <- modelplot +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines1) +
      geom_line(aes(x=t1, y=t2, linetype=Model), data = dd_lines2) +
      scale_linetype_manual(values = "solid",
                            # labels = unname(TeX("$\\beta_0 = 0$")))
                            labels = unname(TeX(paste0("$\\beta_0 = 0, q = $", round(q, digits = q_digits)))))
  }
  
  modelplot <- modelplot +
    geom_point(aes(colour = Genotype_true_AA), size = 1, alpha = 0.15) +
    geom_point(aes(colour = Genotype_true_AA), data = dd_missed, size = 1) +
    scale_colour_manual(name="True genotype",
                        values=c("#F0E442", "#CC79A7", "#56B4E9"),
                        breaks=c("A1A1", "A1A2", "A2A2"),
                        labels=unname(c(TeX("$a_1 a_1$"), TeX("$a_1 a_2$"), TeX("$a_2 a_2$")))) +
    ylim(0, f(max(df$A2.Reads))) +
    # coord_cartesian(xlim=c(0,0.15),ylim=c(0,0.15))+
    xlab(TeX("$t_1$")) +
    ylab(TeX("$t_2$")) +
    # geoms below will use another colour scale
    new_scale_color() +
    geom_point(aes(shape = Missed, colour = Missed), data = dd_missed) +
    scale_shape_manual(values=c(3, 4)) +
    scale_colour_manual(name="Missed",
                        values=c("black", "red"),
                        breaks=c(no_call, "WC"),
                        labels=unname(c(no_call, "WC")))
  
  if (!legend) {modelplot <- modelplot + theme(legend.position="none")}
  if (is.character(title)) {modelplot <- modelplot + labs(title = title)}
  
  return(modelplot)
}


#---------------------------------#
# Bootstrap create sample size ####
#---------------------------------#
dd_boot_create_samplesize <- function(df, repititions, fun="sqrt", INT=TRUE,
                                      seq_start=1, seq_end=F,
                                      b_int=c(0,1,-2), method="Nelder-Mead", control_list=list(maxit=500),
                                      stratify="Marker", seeds=F) {
  if (fun == "sqrt") {
    f <- function(x) sqrt(x)
  } else if (fun == "log") {
    f <- function(x) log(x+1)
  } else if (fun == "I") {
    f <- function(x) x
  }
  
  individuals <- unique(df$ID)
  dilutions <- unique(df$Dilution)
  runs <- unique(df$Run)
  markers <- unique(df$Target.ID)
  # dfs <- select(Dilution, A1.Reads, A2.Reads, Genotype_true_AA) |> split(~Dilution)
  # If df has more than one dilution, the point is to look at the sampling as if
  # a lab has made a dilutions series where they want to fit using all dilutions
  # Hence, in that case we will sample an equal amount of SNPs/markers from each
  # dilution, so each sample size increment corresponds to adding a new individual.
  # This means that if we have e.g. two dilutions, then the increments of sample
  # sizes will be twice the increment as when we have just one dilution.
  # We will stop the sample size when it is as big as our original data, which
  # is the standard sample size when doing bootstrap, i.e. we let the sample
  # size go up to length(individuals)*length(runs)
  if (!is.numeric(seq_end)) {seq_end <- length(individuals)*length(runs)}
  ss <- seq(from=seq_start, to=seq_end, by=1)
  
  if (is.character(stratify)) {
    future_lapply(1:repititions, function(x){
      lapply(ss, function(s){
        # Note: when we sample per marker, we believe that markers can have
        # different drop-out rates, so not all bootstrap samples of size s will
        # become of size s, but only the size reflected in the data.
        df_test <- lapply(markers, function(m) {
          ids <- sample(individuals, size = s, replace = TRUE)
          df_individuals <- filter(df, ID %in% ids)
          # We should see Run and Dilution as independent, since if we can't do
          # that, it implies that we should fit a new model each time we want to
          # use it, which would make it unusable. One could of course split the
          # data by Run and bootstrap for the individual runs to see if the
          # result for each Run is similar to the other runs, and if the results
          # for all runs match those where each run is seen as interchangeable
          # with the others.
          # If we run the command xtabs(~ Run + Dilution, dd_GT), we see that
          # not all runs in the original data contain all individuals, and this
          # makes it quite complicated to make a the intended bootstrap using
          # all runs, so after having spend too much time thinking about it, I
          # conclude that the easiest solution is to filter away Run E and F.
          # (of course the problem can be solved by making a very specialized
          # algorithm specific for this particular dataset, but this is not
          # worth the effort)
          if (length(runs) > 1) {
            df_individuals <- lapply(ids, function(p) {
              lapply(dilutions, function(d) {
                df_dp <- filter(df_individuals, ID == p, Dilution == d)
                run_p <- sample(unique(df_dp$Run), 1)
                filter(df_dp, Run == run_p)
              }) |> bind_rows()
            }) |> bind_rows()
          }
          
          filter(df_individuals, Target.ID == m)
        }) |> bind_rows()
        
        m <- fit_symmetry_model2(df_test, f=f, intercept = INT, method = method, control_list = control_list, b_int = b_int)
        df_test <- predict_prob_threshold(df_test, fit = m, GT_nucleotide = FALSE)
        m$n_WC <- sum(df_test$Genotype_true_AA != df_test$Genotype_pred_AA)
        m$f <- fun
        m$Samplesize <- s
        m$Dilution <- paste0(dilutions, "pg", collapse = " & ")
        m$Stratification <- stratify
        
        return(m)
      })
    }, future.seed = seeds, future.scheduling = 1)
  } else if (is.numeric(stratify)) {
    df_het <- filter(df, Genotype_true == "A1A2")
    df_hom <- filter(df, Genotype_true != "A1A2")
    future_lapply(1:repititions, function(x){
      # Stratifying on heterozygous proportion excludes the possibility to also
      # sample per individual as different individuals have different het_props.
      # Thus the point about stratifying on het_prop is to see if it suffices to
      # distinguish between observations just as het and hom and NOT considering
      # which marker they come from. This also implies that if we have more
      # dilutions and sample a set of markers in one of the dilutions, then we
      # should NOT demand that it is the same markers from the same individuals that
      # we sample in the other dilution, but simply just some set of the same
      # size and same genotype (on the het/hom-level). This does, however, make
      # it difficult to handle drop-outs, but maybe it is most informative NOT
      # to consider drop-outs in the bootstrap for sample size, because
      # different labs and populations may have different drop-out rates, so the
      # rate for this dataset may be of little relevance for others, in which
      # case it may be much more illuminating to know how many individuals with all
      # SNPs observed, that should be used...
      lapply(ss, function(s){
        if (length(stratify) == 2) {stratify <- rbeta(1, stratify[1], stratify[2])}
        df_test <- lapply(dilutions, function(d) {
          df_het_d <- filter(df_het, Dilution == d)
          df_hom_d <- filter(df_hom, Dilution == d)
          idx_het <- sample(nrow(df_het_d), size = round(s*165*stratify), replace = TRUE)
          idx_hom <- sample(nrow(df_hom_d), size = round(s*165*(1-stratify)), replace = TRUE)
          bind_rows(df_het_d[idx_het,], df_hom_d[idx_hom,])
        }) |> bind_rows()
        
        m <- fit_symmetry_model2(df_test, f=f, intercept = INT, method = method, control_list = control_list, b_int = b_int)
        df_test <- predict_prob_threshold(df_test, fit = m, GT_nucleotide = FALSE)
        m$n_WC <- sum(df_test$Genotype_true_AA != df_test$Genotype_pred_AA)
        m$f <- fun
        m$Samplesize <- s
        m$Dilution <- paste0(dilutions, "pg", collapse = " & ")
        m$Stratification <- "Het-prop"
        m$Het_prop <- stratify
        
        return(m)
      })
    }, future.seed = seeds, future.scheduling = 1)
  } else {
    future_lapply(1:repititions, function(x){
      lapply(ss, function(s){
        df_test <- lapply(dilutions, function(d) {
          df_d <- filter(df, Dilution == d)
          idx <- sample(nrow(df_d), size = s*165, replace = TRUE)
          df_d[idx,]
        }) |> bind_rows()
        
        m <- fit_symmetry_model2(df_test, f=f, intercept = INT, method = method, control_list = control_list, b_int = b_int)
        df_test <- predict_prob_threshold(df_test, fit = m, GT_nucleotide = FALSE)
        m$n_WC <- sum(df_test$Genotype_true_AA != df_test$Genotype_pred_AA)
        m$f <- fun
        m$Samplesize <- s
        m$Dilution <- paste0(dilutions, "pg", collapse = " & ")
        m$Stratification <- "None"
        
        return(m)
      })
    }, future.seed = seeds, future.scheduling = 1)
  }
}


#----------------------------#
# Cross-validation create ####
#----------------------------#
# In the bootstrap for accessing the stability of model parameters etc., it is a
# matter of choice, that we want to be able to classify the genotypes solely
# based on the two allele input-signals. Of course we can criticise the
# simplicity of this approach and demand that e.g. the total coverage (for all
# four bases) should be taken into consideration too (not just the two bases
# that are assumed the only possible). If this assumption should be taken
# serious, it should be OK to fit the model on random subsets of SNPs, i.e.
# subsets where we forget about the actual marker, but just care about it being
# some kind of biallelic SNP.
# For validation, it is better to take the actual marker into consideration,
# since if certain markers are more difficult for the kit to measure, and thus
# have a higher rate of wrong calls, then this approach will give a more
# realistic idea about the model's performance.
cross_val_create <- function(df, n_ids_train, by_individual=F,
                             f=f_vst, intercept=T, b_int=c(0,1,-2),
                             method="Nelder-Mead", hessian=F, control_list=list(maxit = 500)) {
  ids <- unique(df$ID)
  # The simplest way to cross-validate would be to split up the data by individual,
  # and then use the full profiles for these individuals:
  if (by_individual) {
    ids_for_train <- sample(ids, size = n_ids_train)
    dd_train <- df |> filter(ID %in% ids_for_train)
    dd_test <- df |> filter(!(ID %in% ids_for_train))
  } else {
    # Another way to do the cross-validation is to create new profiles from the
    # existing ones. This can be done by again splitting by individual, but then do so
    # for each marker. This ensures that if the data contains repeated
    # measurements (dilutions and/or runs), then these are split together.
    markers <- unique(df$Target.ID)
    dd_traintest <- lapply(markers, function(m) {
      dd_marker <- df |> filter(Target.ID == m)
      ids_for_train <- sample(ids, size = n_ids_train)
      # Some markers failed to be measured for all individuals, so this procedure
      # doesn't necessarily sample a number of n_ids_train
      # markers for each marker, which reflects the variability in the data and
      # makes the validation more realistic.
      dd_train <- dd_marker |> filter(ID %in% ids_for_train)
      dd_test <- dd_marker |> filter(!(ID %in% ids_for_train))
      
      return(list("Train" = dd_train, "Test" = dd_test))
    })
    dd_train <- lapply(dd_traintest, function(m) m$Train) |> bind_rows()
    dd_test <- lapply(dd_traintest, function(m) m$Test) |> bind_rows()
  }
  
  m <- fit_symmetry_model2(dd_train, f=f, intercept=intercept, b_int=b_int, method=method, hessian=hessian, control_list=control_list)
  dd_test <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  N_train <- nrow(dd_train)
  N_test <- nrow(dd_test)
  
  dd_test <- dd_test |> mutate(Genotype_EQC = threshold_rule(dd_test))
  HSG_NN <- sum(dd_test$Genotype == no_call)
  HSG_WC <- sum(dd_test$Genotype != no_call & dd_test$Genotype != dd_test$Genotype_true)
  EQC_NN <- sum(dd_test$Genotype_EQC == no_call)
  EQC_WC <- sum(dd_test$Genotype_EQC != no_call & dd_test$Genotype_EQC != dd_test$Genotype_true)
  
  return(list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]], "N_train"=N_train, "N_test"=N_test, "HSG_NN" = HSG_NN, "HSG_WC" = HSG_WC, "EQC_NN" = EQC_NN, "EQC_WC" = EQC_WC))
}


#-----------------------------------#
# Cross-validation create custom ####
#-----------------------------------#
cross_val_create_custom <- function(df, dils_train, dils_test,
                                    f=f_vst, intercept=T, b_int=c(0,1,-2), train_prop=0.75,
                                    method="Nelder-Mead", hessian=F, control_list=list(maxit = 500)) {
  # We will loop trough each marker and dilution to sample n_ids_train SNPs within each
  # marker for the training data. This is done by sampling n_ids_train individuals,
  # but for each marker we will use the same set of individuals for all dilutions
  # in order to resemble a dilution series.
  # In this way, we use the original data as a template to put together n_ids_train new SNP profiles.
  # If the train data uses dilutions from different series, the individuals will be different, so
  # we can't use the same n_ids_train for all dilutions, but have to use two sets of n_ids_train
  # individuals for each dilution series.
  # If df has several runs for the same individual (and same dilution), then we will
  # treat the runs as independent, i.e. simply consider a new run as new individuals
  # by e.g. interaction(Run, ID), which have to be done before using this function,
  # i.e. the dataframe used in the 'df'-argument of this function has to contain this interaction().
  # In cases where the test and train data uses different dilutions, one could argue
  # that all data for the dilution used for testing could be used for out-of-sample
  # estimation, since that data is in principle unseen by the fitted models, BUT
  # that will have the side effect that the test data will always be the same for
  # all of the training splits, which will lead to horizontal lines in our plot and
  # inflates boxplots, which makes them more or less incomparable to the other plots.
  # Therefore, we also want to make a split on the test data.
  # This will actually make the sampling easier
  df_train1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_train)
  df_test1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_test)
  
  df_train2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_train)
  df_test2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_test)
  
  #-----------------------#
  # First dilution series #
  #-----------------------#
  # Finding the random split for each marker
  # We want to loop through all available markers, so we need to look within both
  # the train and test dilutions (some may be missing from one, but present in
  # the other if we are at low dilutions)
  markers <- unique(c(df_train1$Target.ID, df_test1$Target.ID))
  # We only want to grab the individuals in the training data, since we decide
  # that this is where we want the proportion 'train_prop' of the data to be.
  ids_train <- unique(df_train1$ID)
  n_ids_train <- floor(train_prop*length(ids_train))
  
  split_ids <- lapply(markers, function(m) {
    sample(ids_train, size = n_ids_train)
  })
  names(split_ids) <- markers
  # Doing the splitting and combining
  
  df_train1 <- lapply(names(split_ids), function(m) {
    df_train1 |> filter(ID %in% split_ids[[m]], Target.ID == m)
  }) |> bind_rows()
  df_test1 <- lapply(names(split_ids), function(m) {
    df_test1 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
  }) |> bind_rows()
  
  #------------------------#
  # Second dilution series #
  #------------------------#
  # Repeat same steps as for the first dilution series:
  markers <- unique(c(df_train2$Target.ID, df_test2$Target.ID))
  ids_train <- unique(df_train2$ID)
  n_ids_train <- floor(train_prop*length(ids_train))
  
  split_ids <- lapply(markers, function(m) {
    sample(ids_train, size = n_ids_train)
  })
  names(split_ids) <- markers
  
  df_train2 <- lapply(names(split_ids), function(m) {
    df_train2 |> filter(ID %in% split_ids[[m]], Target.ID == m)
  }) |> bind_rows()
  df_test2 <- lapply(names(split_ids), function(m) {
    df_test2 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
  }) |> bind_rows()
  
  # Combining dilution both serie
  df_train <- bind_rows(df_train1, df_train2)
  df_test <- bind_rows(df_test1, df_test2)
  
  # Doing the fitting and testing
  m <- fit_symmetry_model2(df_train, f=f, intercept=intercept, b_int=b_int, method=method, hessian=hessian, control_list=control_list)
  df_test <- predict_prob_threshold(df_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- df_test |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- df_test |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  N_train <- nrow(df_train)
  N_test <- nrow(df_test)
  
  df_test <- df_test |> mutate(Genotype_EQC = threshold_rule(df_test))
  HSG_NN <- sum(df_test$Genotype == no_call)
  HSG_WC <- sum(df_test$Genotype != no_call & df_test$Genotype != df_test$Genotype_true)
  EQC_NN <- sum(df_test$Genotype_EQC == no_call)
  EQC_WC <- sum(df_test$Genotype_EQC != no_call & df_test$Genotype_EQC != df_test$Genotype_true)
  
  return(list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]], "N_train"=N_train, "N_test"=N_test, "HSG_NN" = HSG_NN, "HSG_WC" = HSG_WC, "EQC_NN" = EQC_NN, "EQC_WC" = EQC_WC))
}


#----------------------------------------------------------------------------------------#
# Version outputting p_max for all observations (consumes more space per saved rds-file) #
#----------------------------------------------------------------------------------------#
cross_val_create_custom_all <- function(df, dils_train, dils_test,
                                        f=f_vst, intercept=T, b_int=c(0,1,-2), train_prop=0.75,
                                        method="Nelder-Mead", hessian=F, control_list=list(maxit = 500)) {
  # We will loop trough each marker and dilution to sample n_ids_train SNPs within each
  # marker for the training data. This is done by sampling n_ids_train individuals,
  # but for each marker we will use the same set of individuals for all dilutions
  # in order to resemble a dilution series.
  # In this way, we use the original data as a template to put together n_ids_train new SNP profiles.
  # If the train data uses dilutions from different series, the individuals will be different, so
  # we can't use the same n_ids_train for all dilutions, but have to use two sets of n_ids_train
  # individuals for each dilution series.
  # If df has several runs for the same individual (and same dilution), then we will
  # treat the runs as independent, i.e. simply consider a new run as new individuals
  # by e.g. interaction(Run, ID), which have to be done before using this function,
  # i.e. the dataframe used in the 'df'-argument of this function has to contain this interaction().
  # In cases where the test and train data uses different dilutions, one could argue
  # that all data for the dilution used for testing could be used for out-of-sample
  # estimation, since that data is in principle unseen by the fitted models, BUT
  # that will have the side effect that the test data will always be the same for
  # all of the training splits, which will lead to horizontal lines in our plot and
  # inflates boxplots, which makes them more or less incomparable to the other plots.
  # Therefore, we also want to make a split on the test data.
  # This will actually make the sampling easier
  df_train1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_train)
  df_test1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_test)
  
  df_train2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_train)
  df_test2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_test)
  
  #-----------------------#
  # First dilution series #
  #-----------------------#
  # Finding the random split for each marker
  # We want to loop through all available markers, so we need to look within both
  # the train and test dilutions (some may be missing from one, but present in
  # the other if we are at low dilutions)
  markers <- unique(c(df_train1$Target.ID, df_test1$Target.ID))
  # We only want to grab the individuals in the training data, since we decide
  # that this is where we want the proportion 'train_prop' of the data to be.
  ids_train <- unique(df_train1$ID)
  n_ids_train <- floor(train_prop*length(ids_train))
  
  split_ids <- lapply(markers, function(m) {
    sample(ids_train, size = n_ids_train)
  })
  names(split_ids) <- markers
  # Doing the splitting and combining
  
  df_train1 <- lapply(names(split_ids), function(m) {
    df_train1 |> filter(ID %in% split_ids[[m]], Target.ID == m)
  }) |> bind_rows()
  df_test1 <- lapply(names(split_ids), function(m) {
    df_test1 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
  }) |> bind_rows()
  
  #------------------------#
  # Second dilution series #
  #------------------------#
  # Repeat same steps as for the first dilution series:
  markers <- unique(c(df_train2$Target.ID, df_test2$Target.ID))
  ids_train <- unique(df_train2$ID)
  n_ids_train <- floor(train_prop*length(ids_train))
  
  split_ids <- lapply(markers, function(m) {
    sample(ids_train, size = n_ids_train)
  })
  names(split_ids) <- markers
  
  df_train2 <- lapply(names(split_ids), function(m) {
    df_train2 |> filter(ID %in% split_ids[[m]], Target.ID == m)
  }) |> bind_rows()
  df_test2 <- lapply(names(split_ids), function(m) {
    df_test2 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
  }) |> bind_rows()
  
  # Combining dilution both serie
  df_train <- bind_rows(df_train1, df_train2)
  df_test <- bind_rows(df_test1, df_test2)
  
  # Doing the fitting and testing
  m <- fit_symmetry_model2(df_train, f=f, intercept=intercept, b_int=b_int, method=method, hessian=hessian, control_list=control_list)
  df_test <- predict_prob_threshold(df_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- df_test |> filter(WC) |> select(p_max) |> arrange(p_max)
  p_NN <- df_test |> filter(!WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) == 0) {
    p_WC <- NULL
  }
  N_train <- nrow(df_train)
  N_test <- nrow(df_test)
  
  df_test <- df_test |> mutate(Genotype_EQC = threshold_rule(df_test))
  HSG_NN <- sum(df_test$Genotype == no_call)
  HSG_WC <- sum(df_test$Genotype != no_call & df_test$Genotype != df_test$Genotype_true)
  EQC_NN <- sum(df_test$Genotype_EQC == no_call)
  EQC_WC <- sum(df_test$Genotype_EQC != no_call & df_test$Genotype_EQC != df_test$Genotype_true)
  
  return(list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]], "N_train"=N_train, "N_test"=N_test, "HSG_NN" = HSG_NN, "HSG_WC" = HSG_WC, "EQC_NN" = EQC_NN, "EQC_WC" = EQC_WC))
}


#-------------------------------------------------------------------------------#
# Version where all data is used for test, when fit and test dils are different #
#-------------------------------------------------------------------------------#
cross_val_create_custom2 <- function(df, n_ids_train, dils_train, dils_test,
                                     f=f_vst, intercept=T, b_int=c(0,1,-2),
                                     method="Nelder-Mead", hessian=F, control_list=list(maxit = 500)) {
  # We will loop trough each marker and sample n_ids_train SNPs within each
  # marker for the training data. This is done by sampling n_ids_train individuals,
  # but for each marker we will use the same set of individuals for all dilutions
  # in order to resemble a dilution series.
  # In this way, we use the original data as a template to put together n_ids_train new SNP profiles.
  # If the train data uses dilutions from different series, the individuals will be different, so
  # we can't use the same n_ids_train for all dilutions, but have to use two sets of n_ids_train
  # individuals for each dilution series.
  # If df has several runs for the same individual (and same dilution), then we will
  # treat the runs as independent, i.e. simply consider a new run as new individuals
  # by e.g. interaction(Run, ID), which have to be done before using this function,
  # i.e. the dataframe used in the 'df'-argument of this function has to contain this interaction().
  df_train1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_train)
  df_test1 <- filter(df, Dilution %in% c(1000, 500, 250, 125, 62.5, 31.25), Dilution %in% dils_test)
  
  df_train2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_train)
  df_test2 <- filter(df, Dilution %in% c(50, 25, 12.5, 6.25), Dilution %in% dils_test)
  
  common_dils <- dils_train[dils_train %in% dils_test]
  if (length(common_dils) > 0) {
    # If there is an overlap, it can be found in both df_trainx and df_testx
    df_c1 <- filter(df_train1, Dilution %in% common_dils)
    df_c2 <- filter(df_train2, Dilution %in% common_dils)
    
    #-----------------------#
    # First dilution series #
    #-----------------------#
    # Finding the random split for each marker
    markers <- unique(df_train1$Target.ID)
    ids_train <- unique(df_train1$ID)
    # The overlap should be removed from df_trainx and df_testx, because it has
    # to be randomly split between these sets:
    df_train1 <- df_train1 |> filter(!(Dilution %in% common_dils))
    df_test1 <- df_test1 |> filter(!(Dilution %in% common_dils))
    split_ids <- lapply(markers, function(m) {
      sample(ids_train, size = n_ids_train)
    })
    names(split_ids) <- markers
    # Doing the splitting and combining
    # (remembering that the same individuals should also be used to split the
    # non-overlapping training data, while all of the non-overlapping test data
    # should NOT be split, since it is already out-of-sample data)
    df_train1_c <- lapply(names(split_ids), function(m) {
      df_c1 |> filter(ID %in% split_ids[[m]], Target.ID == m)
    }) |> bind_rows()
    df_train1 <- lapply(names(split_ids), function(m) {
      df_train1 |> filter(ID %in% split_ids[[m]], Target.ID == m)
    }) |> bind_rows() |> bind_rows(df_train1_c)
    
    df_test1 <- df_test1 |> bind_rows(lapply(names(split_ids), function(m) {
      df_c1 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
    }) |> bind_rows())
    
    #------------------------#
    # Second dilution series #
    #------------------------#
    # Repeat same steps as for the first dilution series:
    markers <- unique(df_train2$Target.ID)
    ids_train <- unique(df_train2$ID)
    
    df_train2 <- df_train2 |> filter(!(Dilution %in% common_dils))
    df_test2 <- df_test2 |> filter(!(Dilution %in% common_dils))
    split_ids <- lapply(markers, function(m) {
      sample(ids_train, size = n_ids_train)
    })
    names(split_ids) <- markers
    
    df_train2_c <- lapply(names(split_ids), function(m) {
      df_c2 |> filter(ID %in% split_ids[[m]], Target.ID == m)
    }) |> bind_rows()
    df_train2 <- lapply(names(split_ids), function(m) {
      df_train2 |> filter(ID %in% split_ids[[m]], Target.ID == m)
    }) |> bind_rows() |> bind_rows(df_train2_c)
    
    df_test2 <- df_test2 |> bind_rows(lapply(names(split_ids), function(m) {
      df_c2 |> filter(!(ID %in% split_ids[[m]]), Target.ID == m)
    }) |> bind_rows())
    
    
    df_train <- bind_rows(df_train1, df_train2)
    df_test <- bind_rows(df_test1, df_test2)
    
  } else {
    # If the dilutions to be used for train and test don't have any overlap,
    # then we will use all possible data with the desired dilutions for testing.
    df_test <- bind_rows(df_test1, df_test2)
    # And the random split for the training data is easily done, since the
    # complement of the split shouldn't be transferred to the test data.
    if (nrow(df_train1) > 0) {
      markers <- unique(df_train1$Target.ID)
      ids_train <- unique(df_train1$ID)
      lapply(markers, function(m) {
        ids_for_train <- sample(ids_train, size = n_ids_train)
        return(df_train1 |> filter(ID %in% ids_for_train, Target.ID == m))
      }) |> bind_rows()
    }
    if (nrow(df_train2) > 0) {
      markers <- unique(df_train2$Target.ID)
      ids_train <- unique(df_train2$ID)
      lapply(markers, function(m) {
        ids_for_train <- sample(ids_train, size = n_ids_train)
        return(df_train2 |> filter(ID %in% ids_for_train, Target.ID == m))
      }) |> bind_rows()
    }
    
    df_train <- bind_rows(df_train1, df_train2)
  }
  
  m <- fit_symmetry_model2(df_train, f=f, intercept=intercept, b_int=b_int, method=method, hessian=hessian, control_list=control_list)
  df_test <- predict_prob_threshold(df_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- df_test |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- df_test |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  N_train <- nrow(df_train)
  N_test <- nrow(df_test)
  
  df_test <- df_test |> mutate(Genotype_EQC = threshold_rule(df_test))
  HSG_NN <- sum(df_test$Genotype == no_call)
  HSG_WC <- sum(df_test$Genotype != no_call & df_test$Genotype != df_test$Genotype_true)
  EQC_NN <- sum(df_test$Genotype_EQC == no_call)
  EQC_WC <- sum(df_test$Genotype_EQC != no_call & df_test$Genotype_EQC != df_test$Genotype_true)
  
  return(list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]], "N_train"=N_train, "N_test"=N_test, "HSG_NN" = HSG_NN, "HSG_WC" = HSG_WC, "EQC_NN" = EQC_NN, "EQC_WC" = EQC_WC))
}


#--------------------------------------#
# Version collecting all models in one #
#--------------------------------------#
cross_val_create_custom3 <- function(df, n_ids_train, dils_train, dils_test,
                                    method="Nelder-Mead", hessian=F, control_list=list(maxit = 500)) {
  # If df has several runs for the same individual (and same dilution), then we will
  # treat the runs as independent, i.e. simply consider a new run as new individuals
  # by e.g. interaction(Run, ID)
  if (sum(dils_train %in% dils_test) > 0) {
    common_dils <- dils_train[dils_train %in% dils_test]
    dd_traintest <- lapply(common_dils, function(c_dil) {
      dd_overlap <- df |> filter(Dilution == c_dil)
      markers <- unique(dd_overlap$Target.ID)
      ids_train <- unique(dd_overlap$ID)
      lapply(markers, function(m) {
        ids_for_train <- sample(ids_train, size = n_ids_train)
        d_train <- dd_overlap |> filter(ID %in% ids_for_train, Target.ID == m)
        d_test <- dd_overlap |> filter(!(ID %in% ids_for_train), Target.ID == m)
        return(list("Train" = d_train, "Test" = d_test))
      })
    })
    
    # Independent parts of data
    dd_test <- df |> filter(Dilution %in% setdiff(dils_test, common_dils))
    
    dd_train <- lapply(setdiff(dils_train, common_dils), function(dil) {
      d_train <- df |> filter(Dilution == dil)
      markers <- unique(d_train$Target.ID)
      ids_train <- unique(d_train$ID)
      lapply(markers, function(m) {
        ids_for_train <- sample(ids_train, size = n_ids_train)
        return(d_train |> filter(ID %in% ids_for_train, Target.ID == m))
      }) |> bind_rows()
    }) |> bind_rows()
    
    dd_train <- dd_traintest |> lapply(function(d) {
      lapply(d, function(m) {
        m$Train
      }) |> bind_rows()
    }) |> bind_rows() |> 
      bind_rows(dd_train)
    
    dd_test <- dd_traintest |> lapply(function(d) {
      lapply(d, function(m) {
        m$Test
      }) |> bind_rows()
    }) |> bind_rows() |> 
      bind_rows(dd_test)
  } else {
    # If the dilutions to be used for train and test don't have any overlap,
    # then we will use all possible data with the desired dilutions for testing.
    # I.e. we don't make further filtering for dd_test, but only for dd_train.
    dd_train <- lapply(dils_train, function(dil) {
      dd_train <- df |> filter(Dilution == dil)
      markers <- unique(dd_train$Target.ID)
      ids_train <- unique(dd_train$ID)
      lapply(markers, function(m) {
        ids_for_train <- sample(ids_train, size = n_ids_train)
        return(dd_train |> filter(ID %in% ids_for_train, Target.ID == m))
      }) |> bind_rows()
    }) |> bind_rows()
    
    dd_test <- df |> filter(Dilution %in% dils_test)
  }
  
  
  #--------#
  # MODELS #
  #--------#
  
  # Sqrt NI
  m <- fit_symmetry_model2(dd_train, f=function(x){sqrt(x)}, intercept=F, b_int=c(1,-2), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  N_train <- nrow(dd_train)
  N_test <- nrow(dd_test)
  
  dd_test2 <- dd_test |> mutate(Genotype_EQC = threshold_rule(dd_test))
  HSG_NN <- sum(dd_test2$Genotype == no_call)
  HSG_WC <- sum(dd_test2$Genotype != no_call & dd_test2$Genotype != dd_test2$Genotype_true)
  EQC_NN <- sum(dd_test2$Genotype_EQC == no_call)
  EQC_WC <- sum(dd_test2$Genotype_EQC != no_call & dd_test2$Genotype_EQC != dd_test2$Genotype_true)
  
  l_sqrt_NI <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  # Sqrt
  m <- fit_symmetry_model2(dd_train, f=function(x){sqrt(x)}, intercept=T, b_int=c(0,1,-2), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  
  l_sqrt <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  # Log NI
  m <- fit_symmetry_model2(dd_train, f=function(x){log(x+1)}, intercept=F, b_int=c(1,-2), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  
  l_log_NI <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  # Log
  m <- fit_symmetry_model2(dd_train, f=function(x){log(x+1)}, intercept=T, b_int=c(0,1,-2), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  
  l_log <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  # Id NI
  m <- fit_symmetry_model2(dd_train, f=function(x){x}, intercept=F, b_int=c(0,-1), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  
  l_id_NI <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  # Id
  m <- fit_symmetry_model2(dd_train, f=function(x){x}, intercept=T, b_int=c(0,0,-1), method=method, hessian=hessian, control_list=control_list)
  dd_test2 <- predict_prob_threshold(dd_test, fit = m) |> 
    mutate(WC = Genotype_pred != Genotype_true)
  
  p_WC <- dd_test2 |> filter(WC) |> select(p_max) |> arrange(p_max)
  if (nrow(p_WC) > 0) {
    p_NN <- dd_test2 |> filter(!WC, p_max <= max(p_WC$p_max)) |> select(p_max) |> arrange(p_max)
  } else {
    p_WC <- NULL
    p_NN <- NULL
  }
  
  l_id <- list("fit"=m, "p_WC"=p_WC[[1]], "p_NN"=p_NN[[1]])
  
  return(list("Sqrt_NI" = l_sqrt_NI, "Sqrt" = l_sqrt, "Log_NI" = l_log_NI, "Log" = l_log, "Id_NI" = l_id_NI, "Id" = l_id, "N_train"=N_train, "N_test"=N_test, "HSG_NN" = HSG_NN, "HSG_WC" = HSG_WC, "EQC_NN" = EQC_NN, "EQC_WC" = EQC_WC))
}


#--------------------------------#
# Plot function for bootstrap ####
#--------------------------------#
plot_boot <- function(df, par, dils, complete_sep=T, Y_lim=F, CI=F, Q=F, full_fit=F,
                      complete_theme=theme_minimal(), main_title = F, add_legend = F, x_lab=F, y_lab=F, x_tick=F, y_tick=T,
                      plot_title_size = 7, axis_title_size = 7, axis_ticks_size = 7, legend_title_size = 7, legend_text_size = 7) {
  df <- df |> filter(Parameter == par, Dilution == dils)
  if (is.data.frame(full_fit)) {df <- df |> left_join(full_fit, by = c("Dilution", "Parameter"))}
  df <- df |> group_by(Dilution, Samplesize, Parameter)
  if (is.numeric(CI)) {
    df <- df |> mutate("CI_l" = mean(Estimate)-qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()),
                       "CI_u" = mean(Estimate)+qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()))
  }
  if (is.numeric(Q)) {
    df <- df |> mutate("q_lower" = quantile(Estimate, probs = Q[1]),
                       "q_upper" = quantile(Estimate, probs = Q[2]))
  }
  if (complete_sep %in% c("Yes", "No")) {
    df <- filter(df, `Complete separation` == complete_sep)
  }
  
  p <- ggplot(df, aes(x=Samplesize, y=Estimate)) + 
    complete_theme +
    theme(plot.title = element_text(size = plot_title_size),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_ticks_size),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_text_size))
  if (complete_sep == "Yes") {
    p <- p + geom_point(colour = "#00BFC4")
  } else if (complete_sep == "No") {
    p <- p + geom_point(colour = "#F8766D")
  } else if (complete_sep) {
    if (add_legend) {
      p <- p + geom_point(aes(colour=`Complete separation`))
    } else {
      p <- p + geom_point(aes(colour=`Complete separation`), show.legend = F)
    }
  } else {
    p <- p + geom_point()
  }
  
  
  if (is.data.frame(full_fit) & is.numeric(Q)) {
    df_lines <- df |> pivot_longer(cols = c("q_lower", "q_upper", "Sample_estimate"),
                                   names_to = "Lines",
                                   values_to = "y_value") |> 
      mutate(Lines = case_when(Lines == "q_lower" ~ "Q10",
                              Lines == "q_upper" ~ "Q90",
                              .default = "Fit to all data"))
    if (!add_legend) {
      p <- p + geom_line(aes(x=Samplesize, y=y_value, lty=Lines), colour="black", data = df_lines, show.legend = F)
    } else {
      p <- p + geom_line(aes(x=Samplesize, y=y_value, lty=Lines), colour="black", data = df_lines)
    }
    p <- p + scale_linetype_manual(values = c("Q90" = "dotted", "Fit to all data" = "dashed", "Q10" = "dotted"),
                                   breaks = c("Q90", "Fit to all data", "Q10"))
  } else {
    if (is.data.frame(full_fit)) {
      p <- p + 
        geom_hline(aes(yintercept = Sample_estimate), lty = "dashed", show.legend = F)
    }
    if (is.numeric(Q)) {
      p <- p + 
        geom_line(aes(y = q_lower), lty = "dotted", colour="black", show.legend = F) +
        geom_line(aes(y = q_upper), lty = "dotted", colour="black", show.legend = F)
    }
  }
  
  if (is.numeric(CI)) {
    p <- p + 
      geom_line(aes(y = CI_l), lty = "dotted", colour="black", show.legend = F) +
      geom_line(aes(y = CI_u), lty = "dotted", colour="black", show.legend = F)
  }
  if (is.numeric(Y_lim)) {
    p <- p + coord_cartesian(ylim = Y_lim)
  }
  
  # Title, axis labels and ticks
  if (is.character(main_title)) {
    p <- p + ggtitle(main_title)
  }
  if (is.character(x_lab) | is.expression(x_lab)) {
    p <- p + xlab(x_lab)
  } else {
    p <- p + theme(axis.title.x=element_blank())
  }
  if (is.character(y_lab) | is.expression(y_lab)) {
    p <- p + ylab(y_lab)
  } else {
    p <- p + theme(axis.title.y=element_blank())
  }
  if (is.list(x_tick)) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]], labels = x_tick[["text"]])
  } else if (is.numeric(x_tick)) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]])
  } else if (!x_tick) {
    p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  }
  if (!y_tick) {p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())}
  
  return(p)
}


#----------------------------------------------#
# Plot function for bootstrap with het_prop ####
#----------------------------------------------#
plot_boot_hetprop <- function(df, par, dil, complete_sep=F, add_legend = F, Y_lim=F, CI=F, Q=F, full_fit=F, x_lab=F, y_lab=F, x_tick=F, y_tick=T, main_title = F) {
  df <- df |> filter(Parameter == par, Dilution == dil)
  if (is.data.frame(full_fit)) {
    full_fit <- full_fit |> select(Dilution, Het_prop_pop, Parameter, Sample_estimate)
    df <- df |> left_join(full_fit, by = c("Dilution", "Parameter"))
  }
  df <- df |> group_by(Dilution, Het_prop, Parameter)
  if (is.numeric(CI)) {
    df <- df |> mutate("CI_l" = mean(Estimate)-qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()),
                       "CI_u" = mean(Estimate)+qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()))
  }
  if (is.numeric(Q)) {
    df <- df |> mutate("q_lower" = quantile(Estimate, probs = Q[1]),
                       "q_upper" = quantile(Estimate, probs = Q[2]))
  }
  
  if (complete_sep %in% c("Yes", "No")) {
    df <- filter(df, `Complete separation` == complete_sep)
    p <- ggplot(df, aes(x=Het_prop, y=Estimate)) + geom_point(colour = "#F8766D")
  } else {
    p <- ggplot(df, aes(x=Het_prop, y=Estimate, colour=`Complete separation`)) + geom_point()
    if (!add_legend) {
      p <- p + scale_colour_discrete(guide = "none")
    }
  }
  
  if (is.data.frame(full_fit)) {
    p <- p + 
      geom_point(aes(x = Het_prop_pop, y = Sample_estimate), shape = 4, colour = "red", show.legend = F)
  }
  if (is.numeric(CI)) {
    p <- p + 
      geom_line(aes(y = CI_l), lty = "dotted", colour="black", show.legend = F) +
      geom_line(aes(y = CI_u), lty = "dotted", colour="black", show.legend = F)
  }
  if (is.numeric(Q)) {
    p <- p + 
      geom_line(aes(y = q_lower), lty = "dotted", colour="black", show.legend = F) +
      geom_line(aes(y = q_upper), lty = "dotted", colour="black", show.legend = F)
  }
  if (is.numeric(Y_lim)) {
    p <- p + coord_cartesian(ylim = Y_lim)
  }
  
  # Title, axis labels and ticks
  if (is.character(main_title)) {
    p <- p + ggtitle(main_title)
  }
  if (is.character(x_lab) | is.expression(x_lab)) {
    p <- p + xlab(x_lab)
  } else {
    p <- p + theme(axis.title.x=element_blank())
  }
  if (is.character(y_lab) | is.expression(y_lab)) {
    p <- p + ylab(y_lab)
  } else {
    p <- p + theme(axis.title.y=element_blank())
  }
  if (length(x_tick) == 2) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]], labels = x_tick[["text"]])
  } else if (is.numeric(x_tick)) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]])
  } else if (!x_tick) {
    p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  }
  if (!y_tick) {p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())}
  
  return(p)
}


#-----------------------------------------------------#
# Plot function for comparing optimization methods ####
#-----------------------------------------------------#
plot_methods <- function(df, par, dil, Y_lim=F, CI=F, Q=F, full_fits=F, x_lab=F, y_lab=F, x_tick=F, y_tick=T, main_title = F, add_legend = F) {
  df <- df |> filter(Parameter == par, Dilution == dil)
  # instead of jittering, we manually move the points of each Method:
  m <- unique(df$Method)
  df <- df |> mutate(Samplesize = if_else(Method == m[1],
                                          Samplesize+1/4,
                                          Samplesize-1/4))
  if (is.data.frame(full_fits)) {df <- df |> left_join(full_fits, by = c("Dilution", "Parameter", "Method"))}
  df <- df |> group_by(Dilution, Samplesize, Parameter, Method)
  if (is.numeric(CI)) {
    df <- df |> mutate("CI_l" = mean(Estimate)-qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()),
                       "CI_u" = mean(Estimate)+qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()))
  }
  if (is.numeric(Q)) {
    df <- df |> mutate("q_lower" = quantile(Estimate, probs = Q[1]),
                       "q_upper" = quantile(Estimate, probs = Q[2]))
  }
  
  p <- ggplot(df, aes(x=Samplesize, y=Estimate, colour=Method)) + 
    geom_point(shape=20)
  if (!add_legend) {
    p <- p + scale_colour_discrete(guide = "none")
  }
  
  if (is.data.frame(full_fits)) {
    p <- p + 
      geom_hline(aes(yintercept = Sample_estimate), lty = "dashed", show.legend = F)
  }
  if (is.numeric(CI)) {
    p <- p + 
      geom_line(aes(y = CI_l, lty = Method), colour="grey", show.legend = F) +
      geom_line(aes(y = CI_u, lty = Method), colour="grey", show.legend = F) +
      scale_linetype_manual(values = c("dotted", "dashed"))
  }
  if (is.numeric(Q)) {
    p <- p + 
      geom_line(aes(y = q_lower, lty = Method), colour="grey", show.legend = F) +
      geom_line(aes(y = q_upper, lty = Method), colour="grey", show.legend = F) +
      scale_linetype_manual(values = c("dotted", "dashed"))
  }
  if (is.numeric(Y_lim)) {
    p <- p + coord_cartesian(ylim = Y_lim)
  }
  
  # Title, axis labels and ticks
  if (is.character(main_title)) {
    p <- p + ggtitle(main_title)
  }
  if (is.character(x_lab) | is.expression(x_lab)) {
    p <- p + xlab(x_lab)
  } else {
    p <- p + theme(axis.title.x=element_blank())
  }
  if (is.character(y_lab) | is.expression(y_lab)) {
    p <- p + ylab(y_lab)
  } else {
    p <- p + theme(axis.title.y=element_blank())
  }
  if (length(x_tick) == 2) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]], labels = x_tick[["text"]])
  } else if (is.numeric(x_tick)) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]])
  } else if (!x_tick) {
    p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  }
  if (!y_tick) {p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())}
  
  return(p)
}


#-----------------------------------------------------#
# Plot function for parameters in cross-validation ####
#-----------------------------------------------------#
plot_val <- function(df, par, dil, add_legend = F, Y_lim=F, CI=F, Q=F, full_fit=F, x_lab=F, y_lab=F, x_tick=F, y_tick=T, main_title = F) {
  df <- df |> filter(Parameter == par)
  if (is.character(dil)) {df <- df |> filter(Dilution == dil)}
  if (is.data.frame(full_fit)) {
    if (is.character(dil)) {
      df <- df |> left_join(full_fit, by = c("Dilution", "Parameter"), relationship = "many-to-one")
    } else {
      df <- df |> left_join(full_fit, by = "Parameter", relationship = "many-to-one")
    }
  }
  df <- df |> group_by(Parameter)
  if (is.numeric(CI)) {
    df <- df |> mutate("CI_l" = mean(Estimate)-qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()),
                       "CI_u" = mean(Estimate)+qt((1+CI)/2, n()-1)*sd(Estimate)/sqrt(n()))
  }
  if (is.numeric(Q)) {
    df <- df |> mutate("q_lower" = quantile(Estimate, probs = Q[1]),
                       "q_upper" = quantile(Estimate, probs = Q[2]))
  }
  
  p <- ggplot(df, aes(x=`Validation sample`, y=Estimate)) + geom_point()
  
  if (is.data.frame(full_fit)) {
    p <- p + 
      geom_hline(aes(yintercept = Sample_estimate), lty = "dashed", colour="red", show.legend = F)
  }
  if (is.numeric(CI)) {
    p <- p + 
      geom_hline(aes(yintercept = CI_l), lty = "dotted", colour="red", show.legend = F) +
      geom_hline(aes(yintercept = CI_u), lty = "dotted", colour="red", show.legend = F)
  }
  if (is.numeric(Q)) {
    p <- p + 
      geom_hline(aes(yintercept = q_lower), lty = "dotted", colour="red", show.legend = F) +
      geom_hline(aes(yintercept = q_upper), lty = "dotted", colour="red", show.legend = F)
  }
  if (is.numeric(Y_lim)) {
    p <- p + coord_cartesian(ylim = Y_lim)
  }
  
  # Title, axis labels and ticks
  if (is.character(main_title)) {
    p <- p + ggtitle(main_title)
  }
  if (is.character(x_lab) | is.expression(x_lab)) {
    p <- p + xlab(x_lab)
  } else {
    p <- p + theme(axis.title.x=element_blank())
  }
  if (is.character(y_lab) | is.expression(y_lab)) {
    p <- p + ylab(y_lab)
  } else {
    p <- p + theme(axis.title.y=element_blank())
  }
  if (length(x_tick) == 2) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]], labels = x_tick[["text"]])
  } else if (is.numeric(x_tick)) {
    p <- p + scale_x_continuous(breaks = x_tick[["tick"]])
  } else if (!x_tick) {
    p <- p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())
  }
  if (!y_tick) {p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())}
  
  return(p)
}