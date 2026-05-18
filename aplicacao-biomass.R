################################################################################
# R6DATA APPLICATION: CONFORMAL PREDICTION (PARALLELISED + VISUALISATION)
# HYBRID METHODOLOGY
#
# Overview:
#   This script applies three split-conformal prediction methods to the R6data
#   plant biomass allocation dataset. Responses are three-component compositional
#   vectors (LMF, SMF, RMF) modelled via Dirichlet regression.
#
# Methods compared:
#   1. Quantile Split — conformal score based on max Beta-CDF probit transform.
#   2. HDR Split      — highest density region via analytical tau boundary.
#   3. HDR Grid       — same HDR boundary but area estimated on a regular grid.
#
# Workflow:
#   1.  Load and pre-process R6data.csv.
#   2.  Define score functions, tau solver, and grid builder.
#   3.  Wrap each method in a self-contained function.
#   4.  Run K-fold cross-validation in parallel (foreach / doParallel).
#   5.  Summarise coverage and relative area per method.
#   6.  Visualise prediction regions for four randomly selected test points
#       using the Ternary package.
################################################################################

rm(list = ls())

# ── Required packages ─────────────────────────────────────────────────────────
pkgs <- c(
  "DirichletReg",  # Dirichlet regression model and rdirichlet()
  "tidyverse",     # Data manipulation and pipe operator
  "foreach",       # Parallel loop syntax
  "doParallel",    # Parallel backend registration
  "progress",      # (Optional) progress bar utilities
  "dplyr",         # Explicit namespace use inside foreach workers
  "ggtern",        # Ternary plots (ggplot2-based; loaded for compatibility)
  "gridExtra"      # Multi-panel plot layout helpers
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(DirichletReg)
  library(tidyverse)
  library(foreach)
  library(doParallel)
  library(progress)
  library(ggtern)
  library(gridExtra)
})

# ── Parallel cluster setup ────────────────────────────────────────────────────
# Leave one physical core for the main process.
n_cores <- parallel::detectCores(logical = FALSE) - 1
if (n_cores < 1) n_cores <- 1

cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat(sprintf("Cluster started with %d cores.\n", n_cores))

set.seed(2013)

################################################################################
# SECTION 1 — DATA LOADING AND PRE-PROCESSING (R6DATA)
#
# R6data.csv contains repeated measurements of biomass fraction allocation in
# plants. The three response columns (LMF, SMF, RMF) must form a valid
# compositional vector (non-negative, sum to 1) after cleaning.
#
# Pre-processing steps:
#   - Convert categorical variables to factors with consistent levels.
#   - Scale the Day covariate (z-score) to improve numerical stability.
#   - Remove rows with any missing covariate or response value.
#   - Push zero/negative fractions to eps_num to keep the log-likelihood finite.
#   - Re-normalise rows so the three components always sum to exactly 1.
################################################################################

R6data <- read_csv("R6data.csv")
df     <- R6data

df$Trt     <- as.factor(df$Trt)
df$Species <- as.factor(df$Species)
df$Day     <- as.numeric(scale(df$Day))   # z-score scaling

# Extract the raw compositional matrix and set non-finite values to NA.
Yraw    <- as.matrix(df[, c("LMF", "SMF", "RMF")])
eps_num <- 1e-6
Yraw[!is.finite(Yraw)] <- NA

# Keep only rows where both the response and all model covariates are complete.
vars_check <- c("Day", "lnTDM", "Trt", "Species")
ok   <- complete.cases(Yraw, df[, vars_check])
df   <- df[ok, , drop = FALSE]
Yraw <- Yraw[ok, , drop = FALSE]

# Replace zeros/negatives and re-normalise to the simplex.
Yraw[Yraw <= 0] <- eps_num
Yraw            <- Yraw / rowSums(Yraw)

# Attach cleaned response columns and the DR_data response object.
df$y1 <- Yraw[, 1]
df$y2 <- Yraw[, 2]
df$y3 <- Yraw[, 3]
df$Y  <- DR_data(df[, c("y1", "y2", "y3")])

# Build the modelling data.frame with only the variables needed downstream.
full_data    <- df %>% dplyr::select(Species, Day, Trt, lnTDM, y1, y2, y3)
full_data$Y  <- DR_data(full_data[, c("y1", "y2", "y3")])

# Store factor levels so that train/calibration/test subsets always use the
# same encoding (prevents prediction errors from dropped levels in subsets).
LVL_SPECIES <- levels(full_data$Species)
LVL_TRT     <- levels(full_data$Trt)

cat(sprintf("Dataset ready: N = %d observations.\n", nrow(full_data)))

################################################################################
# SECTION 2 — SCORE FUNCTIONS, TAU SOLVER, AND GRID BUILDER
################################################################################

# ── calc_score_quantile ───────────────────────────────────────────────────────
# Computes the Quantile conformal non-conformity score for each row of 'data'.
#
# For observation i and component j:
#   F_j(y_ij) = pbeta(y_ij, alpha_ij, sum_alpha_i - alpha_ij)   [marginal Beta CDF]
#   z_ij      = qnorm(clamp(F_j(y_ij), 1e-8, 1-1e-8))           [probit transform]
#   score_i   = max_j |z_ij|
#
# Parameters:
#   model  fitted DirichReg object.
#   data   data.frame with covariate and y1/y2/y3 columns.
#
# Returns: numeric vector of length nrow(data).
calc_score_quantile <- function(model, data) {
  if (nrow(data) == 0) return(numeric(0))
  
  pred   <- predict(model, newdata = data, alpha = TRUE)
  alphas <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  
  # Extract response columns dynamically (pattern "y" followed by digits).
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  
  n   <- nrow(Y)
  out <- numeric(n)
  
  for (i in 1:n) {
    a  <- alphas[i, ]
    s0 <- sum(a)
    
    res_j <- numeric(length(a))
    
    for (j in seq_along(a)) {
      # Marginal Beta CDF value for component j.
      p       <- pbeta(Y[i, j], shape1 = a[j], shape2 = s0 - a[j])
      # Probit (normal quantile) transform; clamped to avoid ±Inf.
      res_j[j] <- qnorm(pmin(1 - 1e-8, pmax(1e-8, p)))
    }
    
    out[i] <- max(abs(res_j))
  }
  
  out
}

# ── calc_score_loglik ─────────────────────────────────────────────────────────
# Computes the negative log-likelihood (NLL) conformal score for each row.
#
#   score_i = -log f(y_i | alpha_i)
#           = log B(alpha_i) - sum_j (alpha_ij - 1) log(y_ij)
#
# Parameters:
#   model  fitted DirichReg object.
#   data   data.frame with covariate and y1/y2/y3 columns.
#   eps    small constant for numerical safety in log(y) (default 1e-12).
#
# Returns: numeric vector of length nrow(data).
calc_score_loglik <- function(model, data, eps = 1e-12) {
  if (nrow(data) == 0) return(numeric(0))
  
  pred <- predict(model, newdata = data, alpha = TRUE)
  A    <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  Y <- pmax(Y, eps)
  Y <- Y / rowSums(Y)   # re-normalise after clamping
  
  # Vectorised log-density of the Dirichlet.
  logfy <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A - 1) * log(Y))
  
  -logfy   # return negative log-likelihood as the score
}

# ── get_qhat ──────────────────────────────────────────────────────────────────
# Conformal quantile from a vector of calibration scores.
#
# Formula: qhat = score[ min(ceil((n+1)(1-alpha)), n) ]
#
# Parameters:
#   scores  numeric vector (NAs are removed before sorting).
#   alpha   nominal miscoverage level.
#
# Returns: scalar threshold; Inf if no valid scores are available.
get_qhat <- function(scores, alpha) {
  sc <- sort(na.omit(scores))
  n  <- length(sc)
  
  if (n <= 0) return(Inf)
  
  k <- ceiling((n + 1) * (1 - alpha))
  sc[min(k, n)]
}

# ── .tstar_from_alpha_qhat (internal) ────────────────────────────────────────
# Computes the log-density threshold:
#   tstar = log B(alpha) - qhat
#         = [sum lgamma(alpha_j) - lgamma(sum alpha)] - qhat
#
# Parameters:
#   alpha  numeric vector of Dirichlet parameters.
#   qhat   conformal NLL threshold.
#
# Returns: scalar tstar.
.tstar_from_alpha_qhat <- function(alpha, qhat) {
  (sum(lgamma(alpha)) - lgamma(sum(alpha))) - qhat
}

# ── .solve_rho_for_i (internal) ───────────────────────────────────────────────
# Finds rho > 0 such that the HDR boundary equation for component i is zero,
# using R's uniroot() over the interval [1e-12, 1e12].
#
# The objective function F(rho):
#   F(rho) = w_i*log(rho) + (W-w_i)*log(1+rho) - W*log(w_i*rho + (W-w_i)*(1+rho))
#             + Sww - tstar
# where w_j = alpha_j - 1, W = sum(w), Sww = sum_{w_j>0} w_j*log(w_j).
#
# Parameters:
#   w      numeric vector of shape excesses (alpha - 1).
#   W      sum of w.
#   tstar  log-density threshold.
#   i      integer index of the component being solved for.
#
# Returns: scalar rho root, or NA_real_ if no root is found.
.solve_rho_for_i <- function(w, W, tstar, i) {
  wi0 <- w[i]
  Sww <- sum(ifelse(w > 0, w * log(w), 0))
  
  F_obj <- function(rho) {
    if (rho <= 0) return(-Inf)
    Wmi   <- W - wi0
    denom <- wi0 * rho + Wmi * (1 + rho)
    if (denom <= 0) return(-Inf)
    wi0 * log(rho) + Wmi * log1p(rho) - W * log(denom) + Sww - tstar
  }
  
  tryCatch(
    uniroot(F_obj, c(1e-12, 1e12))$root,
    error = function(e) NA_real_
  )
}

# ── tau_min_simplex ───────────────────────────────────────────────────────────
# Computes the minimum-vertex boundary point tau of the HDR in the 2-simplex.
#
# tau_j is the smallest value of y_j such that the Dirichlet log-density still
# equals the threshold -qhat on the simplex face y_j = tau_j.
#
# Algorithm:
#   1. Compute tstar = log B(alpha) - qhat.
#   2. For each component j, call .solve_rho_for_i to obtain rho_j.
#   3. Convert rho_j to tau_j:
#        tau_j = w_j / [denom * (1 + rho_j)]
#        where denom = w_j/(1+rho_j) + (W-w_j)/rho_j.
#   4. Clamp to [0, 1] and project onto simplex if sum(tau) > 1.
#
# Parameters:
#   alpha  numeric vector of Dirichlet parameters (must satisfy all alpha_j >= 1).
#   qhat   conformal NLL threshold.
#
# Returns: length-K numeric vector of boundary coordinates (all in [0, 1]).
tau_min_simplex <- function(alpha, qhat) {
  alpha <- as.numeric(alpha)
  qhat  <- as.numeric(qhat)
  
  # Method requires all alpha_j >= 1 (mode exists on the interior).
  if (any(alpha < 1)) return(rep(0, length(alpha)))
  
  w     <- alpha - 1
  W     <- sum(w)
  tstar <- .tstar_from_alpha_qhat(alpha, qhat)
  
  tau <- numeric(length(alpha))
  
  for (i in seq_along(alpha)) {
    rho <- .solve_rho_for_i(w, W, tstar, i)
    
    if (is.na(rho)) {
      tau[i] <- 0
    } else {
      denom  <- w[i] / (1 + rho) + (W - w[i]) / rho
      tau[i] <- (1 / denom) * w[i] / (1 + rho)
    }
  }
  
  # Clamp and project onto the simplex.
  tau <- pmax(pmin(tau, 1), 0)
  if (sum(tau) > 1) tau <- tau / (sum(tau) + 1e-12)
  
  tau
}

# ── build_triangle_grid_fixed_step ────────────────────────────────────────────
# Constructs a uniform lattice over the HDR sub-simplex { y : y_j >= tau_j }.
#
# The grid is generated by iterating y1 and y2 with the given step size and
# computing y3 = 1 - y1 - y2, keeping only points where y3 >= tau[3].
#
# Parameters:
#   tau   length-3 numeric vector of boundary coordinates.
#   step  step size between grid lines (e.g. 1 / (GRID_N1 - 1)).
#   eps   small tolerance for boundary comparisons (default 1e-12).
#
# Returns:
#   data.frame with columns y1, y2, y3 (one row per valid grid point), or
#   NULL if no valid points exist or the sub-simplex is degenerate.
build_triangle_grid_fixed_step <- function(tau, step, eps = 1e-12) {
  tau              <- as.numeric(tau)
  tau[!is.finite(tau)] <- 0
  tau              <- pmax(pmin(tau, 1), 0)
  
  S <- sum(tau)
  
  # Degenerate case: tau fills the simplex; return the single normalised point.
  if (S > 1 - eps) {
    y <- if (S > 0) tau / S else c(1/3, 1/3, 1/3)
    return(data.frame(y1 = y[1], y2 = y[2], y3 = y[3]))
  }
  
  y1_min <- tau[1]
  y1_max <- 1 - tau[2] - tau[3]
  
  if (!is.finite(y1_max) || y1_max < y1_min + eps) return(NULL)
  
  # Regular sequence along y1; 'by = step' ensures uniform density.
  y1_seq   <- seq(y1_min, y1_max, by = step)
  out_list <- vector("list", length(y1_seq))
  
  for (k in seq_along(y1_seq)) {
    y1k <- y1_seq[k]
    
    y2_min <- tau[2]
    y2_max <- 1 - tau[3] - y1k
    
    if (!is.finite(y2_max) || y2_max < y2_min + eps) next
    
    y2_seq <- seq(y2_min, y2_max, by = step)
    y3_seq <- 1 - y1k - y2_seq
    
    # Keep only points where y3 >= tau[3].
    keep <- (y3_seq >= tau[3] - eps)
    
    if (any(keep)) {
      out_list[[k]] <- data.frame(
        y1 = rep(y1k, sum(keep)),
        y2 = y2_seq[keep],
        y3 = y3_seq[keep]
      )
    }
  }
  
  grid <- do.call(rbind, out_list)
  if (is.null(grid) || nrow(grid) == 0) return(NULL)
  
  grid
}

# ── score_grid_at_xstar ───────────────────────────────────────────────────────
# Evaluates the NLL score at every grid point, holding covariates fixed at the
# values of a single test observation (xstar_row).
#
# Parameters:
#   model     fitted DirichReg object.
#   xstar_row single-row data.frame of the test observation.
#   grid      data.frame with columns y1, y2, y3 (from build_triangle_grid_fixed_step).
#
# Returns: numeric vector of NLL scores, one per row of grid.
score_grid_at_xstar <- function(model, xstar_row, grid) {
  # Replicate the covariate values (excluding response columns) for every grid point.
  covs    <- xstar_row[rep(1, nrow(grid)), setdiff(colnames(xstar_row), c("y1","y2","y3","Y")), drop = FALSE]
  newdata <- cbind(covs, grid)
  calc_score_loglik(model, newdata)
}

# ── make_mc_newdata_3d ────────────────────────────────────────────────────────
# Builds the Monte Carlo evaluation dataset for area estimation by replicating
# the covariate block of a test observation across B simulated simplex points.
#
# Parameters:
#   te_row   single-row data.frame of the test observation.
#   Ysim     [B x 3] matrix of Dirichlet(1,1,1) samples (uniform on simplex).
#
# Returns: data.frame with B rows containing covariates + y1, y2, y3.
make_mc_newdata_3d <- function(te_row, Ysim) {
  B    <- nrow(Ysim)
  df_y <- data.frame(y1 = Ysim[, 1], y2 = Ysim[, 2], y3 = Ysim[, 3])
  
  cov_names <- setdiff(names(te_row), c("y1","y2","y3","Y"))
  covs      <- te_row[, cov_names, drop = FALSE]
  covs_rep  <- covs[rep(1, B), , drop = FALSE]
  
  cbind(covs_rep, df_y)
}

# ── area_mc_from_mcdata ───────────────────────────────────────────────────────
# Estimates the relative area of a conformal prediction region by Monte Carlo.
#
# Scores each row of mc_newdata (pre-computed uniform simplex samples with
# fixed covariates) and returns the fraction that do not exceed qhat.
#
# Parameters:
#   model      fitted DirichReg object.
#   qhat       conformal threshold.
#   mc_newdata data.frame produced by make_mc_newdata_3d.
#   type       "quantile" or "loglik" — which score function to use.
#
# Returns: scalar in [0, 1] approximating the region's area relative to the
#          simplex.
area_mc_from_mcdata <- function(model, qhat, mc_newdata, type = c("quantile", "loglik")) {
  type <- match.arg(type)
  
  scores_sim <- if (type == "quantile") {
    calc_score_quantile(model, mc_newdata)
  } else {
    calc_score_loglik(model, mc_newdata)
  }
  
  mean(scores_sim <= qhat)
}

# ── area_hdr_split_formula ────────────────────────────────────────────────────
# Analytical relative area of the HDR region in the 2-simplex given tau.
#
# The region { y in simplex : y_j >= tau_j for all j } has area proportional
# to (1 - sum(tau))^2 when K = 3.
#
# Parameters:
#   tau  length-K boundary vector.
#
# Returns: scalar in [0, 1]; 0 if sum(tau) >= 1.
area_hdr_split_formula <- function(tau) {
  s <- sum(tau)
  if (!is.finite(s) || s >= 1) return(0)
  (1 - s)^2
}

################################################################################
# SECTION 3 — METHOD WRAPPERS
#
# Each wrapper receives a fitted model, calibration data, a single test row,
# the miscoverage level alpha, and any method-specific extras. It returns a
# one-row data.frame with columns:
#   method         — character method name.
#   coverage       — logical; TRUE if the test point is inside the region.
#   area_relative  — numeric estimate of the region's area relative to the simplex.
################################################################################

# ── run_quantile ──────────────────────────────────────────────────────────────
# Applies the Quantile split-conformal method to a single test point.
#
# Parameters:
#   model      fitted DirichReg object.
#   cal        calibration data.frame.
#   te_row     single-row test data.frame.
#   alpha      nominal miscoverage level.
#   mc_newdata Monte Carlo evaluation dataset (from make_mc_newdata_3d).
#
# Returns: one-row data.frame.
run_quantile <- function(model, cal, te_row, alpha, mc_newdata) {
  scores   <- calc_score_quantile(model, cal)
  qhat     <- get_qhat(scores, alpha)
  score_te <- calc_score_quantile(model, te_row)
  
  # Area estimated as the fraction of uniform simplex samples with score <= qhat.
  area_val <- area_mc_from_mcdata(model, qhat, mc_newdata, type = "quantile")
  
  data.frame(
    method        = "Quantile",
    coverage      = as.logical(score_te <= qhat),
    area_relative = as.numeric(area_val)
  )
}

# ── run_hdr_split ─────────────────────────────────────────────────────────────
# Applies the HDR Split conformal method to a single test point.
#
# Coverage rule: all components of y_te must be >= the corresponding tau_j.
# Area: analytical formula (1 - sum(tau))^2.
#
# Parameters:
#   model      fitted DirichReg object.
#   cal        calibration data.frame.
#   te_row     single-row test data.frame.
#   alpha      nominal miscoverage level.
#   mc_newdata (unused here; kept for uniform interface).
#
# Returns: one-row data.frame.
run_hdr_split <- function(model, cal, te_row, alpha, mc_newdata) {
  scores <- calc_score_loglik(model, cal)
  qhat   <- get_qhat(scores, alpha)
  
  pred <- predict(model, newdata = te_row, alpha = TRUE)
  a    <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  tau  <- tau_min_simplex(a, qhat)
  
  y       <- as.numeric(te_row[, c("y1","y2","y3")])
  covered <- all(y >= tau)
  
  area_val <- area_hdr_split_formula(tau)
  
  data.frame(
    method        = "HDR_Split",
    coverage      = as.logical(covered),
    area_relative = as.numeric(area_val)
  )
}

# ── run_hdr_grid ──────────────────────────────────────────────────────────────
# Applies the HDR Grid conformal method to a single test point.
#
# Area is estimated by evaluating the NLL score at every grid point inside the
# tau-bounded sub-simplex and computing:
#   area_rel = (1 - sum(tau))^2 * (fraction of grid points with score <= qhat)
#
# Parameters:
#   model  fitted DirichReg object.
#   cal    calibration data.frame.
#   te_row single-row test data.frame.
#   alpha  nominal miscoverage level.
#   step   grid step size.
#
# Returns: one-row data.frame.
run_hdr_grid <- function(model, cal, te_row, alpha, step) {
  scores <- calc_score_loglik(model, cal)
  qhat   <- get_qhat(scores, alpha)
  
  pred <- predict(model, newdata = te_row, alpha = TRUE)
  a    <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  tau  <- tau_min_simplex(a, qhat)
  
  grid_df <- build_triangle_grid_fixed_step(tau, step = step)
  
  if (is.null(grid_df) || nrow(grid_df) == 0) {
    area_rel <- 0
  } else {
    grid_scores  <- score_grid_at_xstar(model, te_row, grid_df)
    accept       <- (grid_scores <= qhat)
    p_hat        <- mean(accept, na.rm = TRUE)
    
    s_tau        <- sum(tau)
    factor_poly  <- if (!is.finite(s_tau) || s_tau >= 1) 0 else (1 - s_tau)^2
    area_rel     <- factor_poly * p_hat
  }
  
  # Coverage uses the NLL score of the actual test observation.
  score_te <- calc_score_loglik(model, te_row)
  
  data.frame(
    method        = "HDR_Grid",
    coverage      = as.logical(score_te <= qhat),
    area_relative = as.numeric(area_rel)
  )
}

################################################################################
# SECTION 4 — K-FOLD CROSS-VALIDATION (PARALLELISED)
#
# The dataset is split into K_FOLDS folds. For each fold k:
#   - Test set   : observations with fold_ids == k.
#   - Training   : 70% of the remaining observations (randomly sampled).
#   - Calibration: the other 30% of remaining observations.
#
# A DirichReg model is fitted on the training set and all three methods are
# evaluated on each test point. Results are combined via rbind across folds.
#
# Key considerations for parallel workers:
#   - Factor levels are explicitly re-applied so that subsets never drop levels.
#   - B_MC_AREA uniform simplex samples are drawn once per fold and reused
#     for all test points in that fold (avoiding B_MC_AREA * N_TEST extra draws).
#   - The model formula includes quadratic Day, lnTDM, and a full interaction
#     of Species, Trt, and (linear) Day in the mean sub-model, with a
#     Species*Day + Trt precision sub-model.
################################################################################

ALPHA_CONFORMAL <- 0.10          # Nominal miscoverage level (target coverage = 90%)
K_FOLDS         <- 20            # Number of cross-validation folds
B_MC_AREA       <- 10000         # Uniform simplex samples for Quantile area estimation
SEED_BASE       <- 2013          # Base seed; each fold uses SEED_BASE + k
GRID_N1         <- 200           # Grid lines along y1 for HDR Grid
GRID_STEP       <- 1 / (GRID_N1 - 1)  # Derived step size

cat(paste0("\n--- Starting Cross-Validation (K=", K_FOLDS, ") ---\n"))

N_TOTAL  <- nrow(full_data)
fold_ids <- sample(rep(1:K_FOLDS, length.out = N_TOTAL))

# Packages to load in each parallel worker.
worker_pkgs <- c("DirichletReg", "tidyverse", "dplyr")

results_df <- foreach(
  k           = 1:K_FOLDS,
  .combine    = rbind,
  .packages   = worker_pkgs
) %dopar% {
  
  set.seed(SEED_BASE + k)
  
  # ── Data splitting ────────────────────────────────────────────────────────
  test_idx     <- which(fold_ids == k)
  rest_idx     <- setdiff(seq_len(N_TOTAL), test_idx)
  n_tr_target  <- floor(0.70 * length(rest_idx))
  train_idx    <- sample(rest_idx, n_tr_target)
  calib_idx    <- setdiff(rest_idx, train_idx)
  
  dat_tr <- full_data[train_idx, , drop = FALSE]
  dat_ca <- full_data[calib_idx, , drop = FALSE]
  dat_te <- full_data[test_idx,  , drop = FALSE]
  
  # Re-apply factor levels to avoid dropped-level prediction errors.
  dat_tr$Species <- factor(dat_tr$Species, levels = LVL_SPECIES)
  dat_ca$Species <- factor(dat_ca$Species, levels = LVL_SPECIES)
  dat_te$Species <- factor(dat_te$Species, levels = LVL_SPECIES)
  dat_tr$Trt     <- factor(dat_tr$Trt, levels = LVL_TRT)
  dat_ca$Trt     <- factor(dat_ca$Trt, levels = LVL_TRT)
  dat_te$Trt     <- factor(dat_te$Trt, levels = LVL_TRT)
  
  # Attach DR_data response object required by DirichReg.
  dat_tr$Y <- DirichletReg::DR_data(dat_tr[, c("y1","y2","y3")])
  dat_ca$Y <- DirichletReg::DR_data(dat_ca[, c("y1","y2","y3")])
  
  # ── Model fitting ─────────────────────────────────────────────────────────
  # Mean sub-model: quadratic Day, lnTDM, full Species*Trt*Day interaction.
  # Precision sub-model: Species*Day + Trt.
  model <- tryCatch({
    DirichReg(
      Y ~ I(Day^2) + lnTDM + Species * Trt * Day | Species * Day + Trt,
      data  = dat_tr,
      model = "alternative"
    )
  }, error = function(e) NULL)
  
  # Skip this fold if the model failed to converge.
  if (is.null(model)) return(NULL)
  
  # Pre-generate uniform simplex samples once and reuse for all test points.
  Ysim <- DirichletReg::rdirichlet(B_MC_AREA, c(1, 1, 1))
  
  # ── Evaluate all three methods on each test point ─────────────────────────
  fold_res <- list()
  
  for (i in 1:nrow(dat_te)) {
    pt <- dat_te[i, , drop = FALSE]
    
    # Skip observations with missing factor levels.
    if (any(is.na(pt$Species)) || any(is.na(pt$Trt))) next
    
    mc_newdata <- make_mc_newdata_3d(pt, Ysim)
    
    r1 <- run_quantile(model, dat_ca, pt, ALPHA_CONFORMAL, mc_newdata)
    r2 <- run_hdr_split(model, dat_ca, pt, ALPHA_CONFORMAL, mc_newdata)
    r3 <- run_hdr_grid(model, dat_ca, pt, ALPHA_CONFORMAL, step = GRID_STEP)
    
    fold_res[[i]] <- rbind(r1, r2, r3)
  }
  
  do.call(rbind, fold_res)
}

# Shut down the parallel cluster after the cross-validation loop.
stopCluster(cl)

################################################################################
# SECTION 5 — FINAL RESULTS SUMMARY
#
# Aggregate coverage and relative area across all folds and test points,
# grouped by method.
################################################################################

summary_tab <- results_df %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    Coverage      = mean(coverage,      na.rm = TRUE),
    Relative_Area = mean(area_relative, na.rm = TRUE)
  )

cat("\n=======================================================\n")
cat("FINAL RESULTS (R6data) — PARALLEL K-FOLD CV\n")
print(summary_tab)

################################################################################
# SECTION 6 — VISUALISATION: PREDICTION REGIONS FOR FOUR TEST POINTS
#
# Four observations are selected at random from full_data. A single
# train/calibration split (70/30 of the non-test rows) is fitted and used to
# generate ternary-plot visualisations of the three prediction regions.
#
# Two sets of 2×2 panels are produced:
#   Panel set 1 — Quantile method: filled polygon (convex hull of accepted MC points).
#   Panel set 2 — HDR methods:     gradient-coloured grid (accepted HDR Grid points)
#                                  with the HDR Split triangle boundary overlaid.
#
# Axis mapping (triangular plot, point = "up"):
#   Top vertex  (a)  — SMF (y2)
#   Left vertex (b)  — LMF (y1)
#   Right vertex (c) — RMF (y3)
# All TernaryXxx() calls therefore supply coordinates in the order (y2, y1, y3).
################################################################################

cat("\n--- Generating prediction region plots for 4 points (Ternary package) ---\n")

# Install and load the Ternary package if needed.
if (!"Ternary" %in% installed.packages()[, "Package"]) install.packages("Ternary")
library(Ternary)

# ── Select four test points ───────────────────────────────────────────────────
set.seed(456)
valid_indices <- sample(1:nrow(full_data), 4)

# ── Single model fit for visualisation ───────────────────────────────────────
idx_viz_test  <- valid_indices
idx_viz_rest  <- setdiff(1:nrow(full_data), idx_viz_test)
n_tr_viz      <- floor(0.70 * length(idx_viz_rest))
idx_viz_train <- sample(idx_viz_rest, n_tr_viz)
idx_viz_calib <- setdiff(idx_viz_rest, idx_viz_train)

dat_tr_viz <- full_data[idx_viz_train, ]
dat_tr_viz$Y <- DR_data(dat_tr_viz[, c("y1","y2","y3")])

dat_ca_viz <- full_data[idx_viz_calib, ]
dat_ca_viz$Y <- DR_data(dat_ca_viz[, c("y1","y2","y3")])

dat_te_viz <- full_data[idx_viz_test, ]

model_viz <- DirichReg(
  Y ~ I(Day^2) + lnTDM + Species * Trt * Day | Species * Day + Trt,
  data  = dat_tr_viz,
  model = "alternative"
)

# ── Pre-compute conformal thresholds for both score types ─────────────────────
sc_qt_cal  <- calc_score_quantile(model_viz, dat_ca_viz)
qhat_qt    <- get_qhat(sc_qt_cal,  ALPHA_CONFORMAL)

sc_hdr_cal <- calc_score_loglik(model_viz, dat_ca_viz)
qhat_hdr   <- get_qhat(sc_hdr_cal, ALPHA_CONFORMAL)

# ── generate_plot_data ────────────────────────────────────────────────────────
# Prepares all data structures needed to draw one test point's prediction
# regions in ternary coordinates.
#
# Parameters:
#   pt_idx  row index in full_data of the test observation.
#
# Returns: named list with elements
#   $qt       — data.frame of convex-hull vertices for the Quantile region.
#   $split    — data.frame of the four HDR Split triangle vertices (closed).
#   $grid     — data.frame of accepted HDR Grid points with a 'score' column.
#   $true_pt  — single-row data.frame with columns y1, y2, y3.
generate_plot_data <- function(pt_idx) {
  pt <- full_data[pt_idx, , drop = FALSE]
  
  # ── Quantile region ────────────────────────────────────────────────────────
  # Simulate 50 000 uniform simplex points, score each one, and keep those
  # whose score does not exceed qhat_qt.
  Ysim    <- DirichletReg::rdirichlet(50000, c(1, 1, 1))
  mc_data <- make_mc_newdata_3d(pt, Ysim)
  
  scores_qt   <- calc_score_quantile(model_viz, mc_data)
  accepted_qt <- Ysim[scores_qt <= qhat_qt, ]
  
  df_qt_poly <- data.frame(
    y1 = accepted_qt[, 1],
    y2 = accepted_qt[, 2],
    y3 = accepted_qt[, 3]
  )
  
  # Compute convex hull in (y1, y2) space to draw a clean polygon.
  # (The accepted simplex cloud is convex for the Quantile score.)
  if (nrow(df_qt_poly) > 0) {
    hull_idx   <- chull(df_qt_poly[, c("y1","y2")])
    df_qt_poly <- df_qt_poly[c(hull_idx, hull_idx[1]), ]  # close the polygon
  }
  
  # ── HDR Split region ───────────────────────────────────────────────────────
  # The HDR Split boundary is a triangle with vertices obtained by sending
  # each component to its maximum while holding the others at tau.
  pred <- predict(model_viz, newdata = pt, alpha = TRUE)
  a    <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  tau  <- tau_min_simplex(a, qhat_hdr)
  
  # Three simplex vertices of the HDR triangle (one per corner of the region).
  v1 <- c(1 - tau[2] - tau[3], tau[2],             tau[3])
  v2 <- c(tau[1],               1 - tau[1] - tau[3], tau[3])
  v3 <- c(tau[1],               tau[2],              1 - tau[1] - tau[2])
  
  df_split <- data.frame(rbind(v1, v2, v3, v1))  # close the polygon
  colnames(df_split) <- c("y1","y2","y3")
  
  # ── HDR Grid region ────────────────────────────────────────────────────────
  # Build a regular grid inside the tau-bounded sub-simplex with a finer step
  # (1/120) than the CV step to get a smooth visual appearance.
  grid_df <- build_triangle_grid_fixed_step(tau, step = 1/120)
  
  if (!is.null(grid_df) && nrow(grid_df) > 0) {
    grid_scores  <- score_grid_at_xstar(model_viz, pt, grid_df)
    grid_df$score <- grid_scores
    grid_df       <- grid_df[grid_df$score <= qhat_hdr, ]  # keep accepted points
  }
  
  list(
    qt       = df_qt_poly,
    split    = df_split,
    grid     = grid_df,
    true_pt  = pt[, c("y1","y2","y3")]
  )
}

# ── score_to_color ────────────────────────────────────────────────────────────
# Maps a numeric score vector to a colour gradient for HDR Grid visualisation.
#
# Scores are linearly normalised to [0, 1] and mapped through a
# yellow → blue → green colour ramp (mimics the target figure style).
#
# Parameters:
#   scores  numeric vector.
#
# Returns: character vector of hex colour codes.
score_to_color <- function(scores) {
  if (length(scores) == 0) return(character(0))
  
  s_min <- min(scores)
  s_max <- max(scores)
  
  if (s_min == s_max) return(rep("blue", length(scores)))
  
  norm_score <- (scores - s_min) / (s_max - s_min)
  color_fun  <- colorRamp(c("yellow", "blue", "green"))
  rgb_vals   <- color_fun(norm_score)
  
  rgb(rgb_vals, maxColorValue = 300)
}

# ── Panel set 1: Quantile method (2×2 grid) ───────────────────────────────────
# Each panel shows the convex hull of accepted Monte Carlo points (shaded red)
# and the true test observation (cyan dot).
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))

for (idx in idx_viz_test) {
  pdata <- generate_plot_data(idx)
  
  # Ternary plot scaffold.
  # point = "up": top vertex = a-axis (SMF), left = b-axis (LMF), right = c-axis (RMF).
  TernaryPlot(
    alab        = "SMF (Stem) \u2192",
    blab        = "LMF (Leaf) \u2192",
    clab        = "\u2190 RMF (Root)",
    point       = "up",
    lab.cex     = 0.8,
    grid.lines  = 5,
    grid.col    = "grey80",
    axis.col    = "grey60",
    axis.labels = seq(0, 100, by = 20),
    main        = paste0("Test #", idx)
  )
  
  # Data order for TernaryXxx(): c(a, b, c) = c(y2, y1, y3).
  if (nrow(pdata$qt) > 0) {
    TernaryPolygon(
      pdata$qt[, c("y2","y1","y3")],
      col    = rgb(0.9, 0.6, 0.6, 0.7),  # semi-transparent red
      border = "red"
    )
  }
  
  # True test point (cyan dot).
  TernaryPoints(pdata$true_pt[, c("y2","y1","y3")], pch = 20, col = "cyan4", cex = 1.2)
}

# ── Panel set 2: HDR methods (2×2 grid) ──────────────────────────────────────
# Each panel shows:
#   - Gradient-coloured accepted HDR Grid points (yellow→blue→green by score).
#   - Dark-blue triangle border for the HDR Split region.
#   - True test point (red dot).
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))

for (idx in idx_viz_test) {
  pdata <- generate_plot_data(idx)
  
  TernaryPlot(
    alab        = "SMF (Stem) \u2192",
    blab        = "LMF (Leaf) \u2192",
    clab        = "\u2190 RMF (Root)",
    point       = "up",
    lab.cex     = 0.8,
    grid.lines  = 5,
    grid.col    = "grey80",
    axis.col    = "grey60",
    axis.labels = seq(0, 100, by = 20),
    main        = paste0("Test #", idx)
  )
  
  # HDR Grid: gradient-coloured accepted interior points.
  # Data order: c(y2, y1, y3).
  if (!is.null(pdata$grid) && nrow(pdata$grid) > 0) {
    cols <- score_to_color(pdata$grid$score)
    TernaryPoints(pdata$grid[, c("y2","y1","y3")], col = cols, pch = 17, cex = 0.6)
  }
  
  # HDR Split: triangle boundary (no fill).
  TernaryPolygon(pdata$split[, c("y2","y1","y3")], border = "darkblue", lwd = 2, col = NA)
  
  # True test point (red dot).
  TernaryPoints(pdata$true_pt[, c("y2","y1","y3")], pch = 20, col = "red", cex = 1.2)
}

# Restore single-panel layout.
par(mfrow = c(1, 1))
