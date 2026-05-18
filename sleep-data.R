################################################################################
# SLEEP DATA APPLICATION (4D): CONFORMAL PREDICTION (PARALLELISED)
# HYBRID METHODOLOGY ADAPTED - GRID VOLUME (DELTA STEP)
################################################################################

rm(list = ls())

# ── 1. REQUIRED PACKAGES ──────────────────────────────────────────────────────
pkgs <- c("DirichletReg", "tidyverse", "foreach", "doParallel", "progress", "readxl", "dplyr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(DirichletReg)
  library(tidyverse)
  library(readxl)
  library(foreach)
  library(doParallel)
  library(progress)
})

# ── Parallel cluster setup ───────────────────────────────────────────────────
# Leave one physical core for main process.
n_cores <- parallel::detectCores(logical = FALSE) - 1
if(n_cores < 1) n_cores <- 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat(sprintf("Cluster started with %d cores.\n", n_cores))
set.seed(40)

################################################################################
# SECTION 2 — DATA LOADING AND PRE-PROCESSING
################################################################################
# Load Excel file and preprocess covariates and 4D compositional responses.
# Handles zero/negative values, normalization to simplex, and factor encoding.

path_file <- "C:/Users/Lucas/Desktop/CP-dirichlet/aplicacao/sleep/dados-boot.xlsx"
if(!file.exists(path_file)){
  if(interactive()) try(path_file <- file.choose())
}
if(file.exists(path_file)){
  df <- read_excel(path_file)
} else stop("Data file not found. Check the path.")

# Covariate processing
df$CaseControl <- as.factor(df$CaseControl)
df$TST         <- as.numeric(df$TST)

# Compose 4D response matrix
Yraw <- as.matrix(df[, c("N1_PERC", "N2_PERC", "N3_PERC", "REM_PERC")])
eps_num <- 1e-6
Yraw[!is.finite(Yraw)] <- NA

vars_check <- c("CaseControl", "TST")
ok <- complete.cases(Yraw, df[, vars_check])
df    <- df[ok, , drop = FALSE]
Yraw <- Yraw[ok, , drop = FALSE]

# Replace zeros/negatives and normalize to simplex
Yraw[Yraw <= 0] <- eps_num
Yraw <- Yraw / rowSums(Yraw)

# Attach compositional columns and DirichletReg response object
df$y1 <- Yraw[,1]; df$y2 <- Yraw[,2]; df$y3 <- Yraw[,3]; df$y4 <- Yraw[,4]
df$Y  <- DR_data(df[, c("y1","y2","y3","y4")])

cat(sprintf("Original dataset: N = %d observations (3D simplex / tetrahedron).\n", nrow(df)))

################################################################################
# SECTION 3 — DATA AUGMENTATION (SYNTHETIC OBSERVATIONS)
################################################################################
# Fit a preliminary Dirichlet model to generate additional synthetic points
# to reach N_TARGET observations via Dirichlet sampling.

fit_gen <- DirichReg(Y ~ TST | TST, data = df, model = "alternative")
mu_hat   <- fitted(fit_gen)
coefs_phi <- fit_gen$coefficients[7:8]
phi_hat_vec <- as.numeric(exp(cbind(1, df$TST) %*% coefs_phi))

N_TARGET <- 1000
n_orig   <- nrow(df)
n_new    <- N_TARGET - n_orig

if(n_new > 0){
  cat(sprintf("Generating %d synthetic observations...\n", n_new))
  idx_rep <- sample(seq_len(n_orig), size = n_new, replace = TRUE)
  df_new  <- df[idx_rep, , drop = FALSE]
  mu_new  <- mu_hat[idx_rep, ]
  phi_new <- phi_hat_vec[idx_rep]
  
  Y_new_raw <- t(sapply(1:n_new, function(i){
    DirichletReg::rdirichlet(1, alpha = phi_new[i] * mu_new[i, ])
  }))
  
  df_new$y1 <- Y_new_raw[,1]; df_new$y2 <- Y_new_raw[,2]
  df_new$y3 <- Y_new_raw[,3]; df_new$y4 <- Y_new_raw[,4]
  df_new$Y  <- DR_data(df_new[, c("y1","y2","y3","y4")])
  
  full_data_raw <- rbind(df, df_new)
} else full_data_raw <- df

full_data <- full_data_raw %>% dplyr::select(CaseControl, TST, y1, y2, y3, y4)
full_data$Y <- DR_data(full_data[, c("y1","y2","y3","y4")])

cat(sprintf("Final dataset ready: N = %d\n", nrow(full_data)))

################################################################################
# SECTION 4 — FUNCTIONS, SCORES, AND GRID (ADAPTED TO 4D)
################################################################################

# ── Quantile score calculation
calc_score_quantile <- function(model, data){
  if(nrow(data) == 0) return(numeric(0))
  pred <- predict(model, newdata = data, alpha = TRUE)
  alphas <- if(is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  n <- nrow(Y); out <- numeric(n)
  for(i in 1:n){
    a <- alphas[i,]; s0 <- sum(a)
    res_j <- numeric(length(a))
    for(j in seq_along(a)){
      p <- pbeta(Y[i,j], shape1 = a[j], shape2 = s0 - a[j])
      res_j[j] <- qnorm(pmin(1-1e-8, pmax(1e-8, p)))
    }
    out[i] <- max(abs(res_j))
  }
  out
}

# ── Negative log-likelihood score
calc_score_loglik <- function(model, data, eps=1e-12){
  if(nrow(data) == 0) return(numeric(0))
  pred <- predict(model, newdata = data, alpha = TRUE)
  A <- if(is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  Y <- pmax(Y, eps); Y <- Y / rowSums(Y)
  logfy <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A-1)*log(Y))
  -logfy
}

# ── Conformal quantile threshold
get_qhat <- function(scores, alpha){
  sc <- sort(na.omit(scores)); n <- length(sc)
  if(n <= 0) return(Inf)
  k <- ceiling((n+1)*(1-alpha))
  sc[min(k,n)]
}

# ── Tau-related helper functions for HDR
.tstar_from_alpha_qhat <- function(alpha, qhat){ (sum(lgamma(alpha)) - lgamma(sum(alpha))) - qhat }
.solve_rho_for_i <- function(w, W, tstar, i){
  wi0 <- w[i]; Sww <- sum(ifelse(w > 0, w*log(w), 0))
  F_obj <- function(rho){
    if(rho <= 0) return(-Inf)
    Wmi <- W - wi0; denom <- wi0*rho + Wmi*(1+rho)
    if(denom <= 0) return(-Inf)
    wi0*log(rho) + Wmi*log1p(rho) - W*log(denom) + Sww - tstar
  }
  tryCatch(uniroot(F_obj, c(1e-12, 1e12))$root, error=function(e) NA_real_)
}
tau_min_simplex <- function(alpha, qhat){
  alpha <- as.numeric(alpha); qhat <- as.numeric(qhat)
  if(any(alpha<1)) return(rep(0,length(alpha)))
  w <- alpha-1; W <- sum(w); tstar <- .tstar_from_alpha_qhat(alpha, qhat)
  tau <- numeric(length(alpha))
  for(i in seq_along(alpha)){
    rho <- .solve_rho_for_i(w,W,tstar,i)
    if(is.na(rho)) tau[i] <- 0
    else { denom <- w[i]/(1+rho) + (W-w[i])/rho; tau[i] <- (1/denom)*w[i]/(1+rho) }
  }
  tau <- pmax(pmin(tau,1),0)
  if(sum(tau)>1) tau <- tau/(sum(tau)+1e-12)
  tau
}

# ── Recursive 4D grid builder with delta step
build_recursive_grid_4d <- function(tau, step, eps=1e-10){
  if(sum(tau) > 1 + eps) return(NULL)
  y1_min <- tau[1]; y1_max <- 1 - sum(tau[2:4])
  if(y1_max < y1_min - eps) return(NULL)
  y1_seq <- seq(y1_min, y1_max, by=step)
  
  grid_list <- vector("list", length(y1_seq))
  for(i in seq_along(y1_seq)){
    val_y1 <- y1_seq[i]
    y2_min <- tau[2]; y2_max <- 1 - sum(tau[3:4]) - val_y1
    if(y2_max >= y2_min - eps){
      y2_seq <- seq(y2_min, y2_max, by=step)
      sub_list <- vector("list", length(y2_seq))
      for(j in seq_along(y2_seq)){
        val_y2 <- y2_seq[j]
        y3_min <- tau[3]; y3_max <- 1 - tau[4] - val_y1 - val_y2
        if(y3_max >= y3_min - eps){
          y3_seq <- seq(y3_min, y3_max, by=step)
          val_y4 <- 1 - val_y1 - val_y2 - y3_seq
          valid <- val_y4 >= tau[4] - eps
          if(any(valid)) sub_list[[j]] <- cbind(y1=val_y1, y2=val_y2, y3=y3_seq[valid], y4=val_y4[valid])
        }
      }
      grid_list[[i]] <- do.call(rbind, sub_list)
    }
  }
  final_grid <- do.call(rbind, grid_list)
  if(is.null(final_grid)) return(NULL)
  as.data.frame(final_grid)
}

# ── Evaluate log-likelihood scores at each grid point for a test observation
score_grid_at_xstar <- function(model, xstar_row, grid){
  covs <- xstar_row[rep(1, nrow(grid)), setdiff(names(xstar_row), c("y1","y2","y3","y4","Y")), drop=FALSE]
  newdata <- cbind(covs, grid)
  newdata$Y <- DR_data(newdata[, c("y1","y2","y3","y4")])
  calc_score_loglik(model, newdata)
}

# ── Monte Carlo generation for 4D volume estimation
make_mc_newdata_4d <- function(te_row, Ysim){
  B <- nrow(Ysim)
  df_y <- data.frame(y1=Ysim[,1], y2=Ysim[,2], y3=Ysim[,3], y4=Ysim[,4])
  cov_names <- setdiff(names(te_row), c("y1","y2","y3","y4","Y"))
  covs <- te_row[, cov_names, drop=FALSE]
  covs_rep <- covs[rep(1,B), , drop=FALSE]
  cbind(covs_rep, df_y)
}

# ── Monte Carlo volume estimation
volume_mc_from_mcdata <- function(model, qhat, mc_newdata, type=c("quantile","loglik")){
  type <- match.arg(type)
  scores_sim <- if(type=="quantile") calc_score_quantile(model, mc_newdata) else calc_score_loglik(model, mc_newdata)
  mean(scores_sim <= qhat)
}

# ── Analytical 4D HDR volume
volume_hdr_split_formula_4d <- function(tau){
  s <- sum(tau); if(!is.finite(s) || s>=1) return(0); (1-s)^3
}

################################################################################
# SECTION 5 — METHOD WRAPPERS
################################################################################

run_quantile <- function(model, cal, te_row, alpha, mc_newdata){
  scores <- calc_score_quantile(model, cal)
  qhat <- get_qhat(scores, alpha)
  score_te <- calc_score_quantile(model, te_row)
  vol_val <- volume_mc_from_mcdata(model, qhat, mc_newdata, type="quantile")
  
  data.frame(method="Quantile", coverage=as.logical(score_te <= qhat), volume_relative=as.numeric(vol_val))
}

run_hdr_split <- function(model, cal, te_row, alpha, mc_newdata){
  scores <- calc_score_loglik(model, cal)
  qhat <- get_qhat(scores, alpha)
  pred <- predict(model, newdata=te_row, alpha=TRUE)
  a <- as.numeric(if(is.list(pred)) pred$alpha else pred)
  tau <- tau_min_simplex(a, qhat)
  y <- as.numeric(te_row[, c("y1","y2","y3","y4")])
  covered <- all(y >= tau)
  vol_val <- volume_hdr_split_formula_4d(tau)
  
  data.frame(method="HDR_Split", coverage=as.logical(covered), volume_relative=as.numeric(vol_val))
}

run_hdr_grid <- function(model, cal, te_row, alpha, step){
  scores <- calc_score_loglik(model, cal)
  qhat <- get_qhat(scores, alpha)
  pred <- predict(model, newdata=te_row, alpha=TRUE)
  a <- as.numeric(if(is.list(pred)) pred$alpha else pred)
  tau <- tau_min_simplex(a, qhat)
  
  grid_df <- build_recursive_grid_4d(tau, step=step)
  vol_val <- 0
  if(!is.null(grid_df) && nrow(grid_df) > 0){
    grid_scores <- score_grid_at_xstar(model, te_row, grid_df)
    ratio <- sum(grid_scores <= qhat) / nrow(grid_df)
    s_tau <- sum(tau); scale_factor <- if(s_tau >= 1) 0 else (1-s_tau)^3
    vol_val <- ratio * scale_factor
  }
  score_te <- calc_score_loglik(model, te_row)
  data.frame(method="HDR_Grid", coverage=as.logical(score_te <= qhat), volume_relative=as.numeric(vol_val))
}

################################################################################
# SECTION 6 — CROSS-VALIDATION LOOP (PARALLEL)
################################################################################

ALPHA_CONFORMAL <- 0.10
K_FOLDS <- 10
B_MC_AREA <- 10000
GRID_STEP <- 0.005 # Fine step for 4D grid

cat(paste0("\n--- Starting Cross-Validation (K=", K_FOLDS, ") ---\n"))
cat(paste0("--- Grid Step: ", GRID_STEP, " ---\n"))

N_TOTAL <- nrow(full_data)
fold_ids <- sample(rep(1:K_FOLDS, length.out=N_TOTAL))
worker_pkgs <- c("DirichletReg", "tidyverse", "dplyr")

results_df <- foreach(k=1:K_FOLDS, .combine=rbind, .packages=worker_pkgs) %dopar% {
  set.seed(13+k)
  test_idx <- which(fold_ids == k)
  rest_idx <- setdiff(seq_len(N_TOTAL), test_idx)
  n_tr_target <- floor(0.70*N_TOTAL)
  train_idx <- sample(rest_idx, n_tr_target)
  calib_idx <- setdiff(rest_idx, train_idx)
  
  dat_tr <- full_data[train_idx, , drop=FALSE]
  dat_ca <- full_data[calib_idx, , drop=FALSE]
  dat_te <- full_data[test_idx, , drop=FALSE]
  
  dat_tr$Y <- DR_data(dat_tr[, c("y1","y2","y3","y4")])
  dat_ca$Y <- DR_data(dat_ca[, c("y1","y2","y3","y4")])
  
  model <- tryCatch({ DirichReg(Y ~ TST | TST, data=dat_tr, model="alternative") }, error=function(e) NULL)
  if(is.null(model)) return(NULL)
  
  Ysim <- DirichletReg::rdirichlet(B_MC_AREA, c(1,1,1,1))
  
  fold_res <- list()
  for(i in 1:nrow(dat_te)){
    pt <- dat_te[i,, drop=FALSE]
    mc_newdata <- make_mc_newdata_4d(pt, Ysim)
    r1 <- run_quantile(model, dat_ca, pt, ALPHA_CONFORMAL, mc_newdata)
    r2 <- run_hdr_split(model, dat_ca, pt, ALPHA_CONFORMAL, mc_newdata)
    r3 <- run_hdr_grid(model, dat_ca, pt, ALPHA_CONFORMAL, step=GRID_STEP)
    fold_res[[i]] <- rbind(r1,r2,r3)
  }
  do.call(rbind, fold_res)
}

stopCluster(cl)

################################################################################
# SECTION 7 — FINAL RESULTS SUMMARY
################################################################################

summary_tab <- results_df %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(Coverage = mean(coverage, na.rm=TRUE),
                   Relative_Volume = mean(volume_relative, na.rm=TRUE))

cat("\n=======================================================\n")
cat("FINAL RESULTS (SLEEP DATA 4D) - DELTA GRID STEP\n")
print(summary_tab)
