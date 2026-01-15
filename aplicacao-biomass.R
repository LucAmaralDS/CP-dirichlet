# =================================================================
# ANEXO: Análise de Regressão Dirichlet e Predição Conformal
# (tau fechado + uniroot + NA->0) + build_triangle_grid estável
# =================================================================

rm(list = ls())

# ----------------------------------------------------------
# 1. PACOTES E CONFIGURAÇÕES INICIAIS
# ----------------------------------------------------------
pkgs <- c("DirichletReg", "progress", "tidyverse",
          "ggtern", "gridExtra", "grid", "viridis", "dplyr")

for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

set.seed(2013)

# Parâmetros para Predição Conformal
ALPHA_CONFORMAL <- 0.10
K_FOLDS         <- 20
EPS_LOG         <- 1e-9
GRID_N1         <- 100
GRID_N2         <- 100

# Parâmetros para Diagnóstico (Bootstrap/Envelope)
DO_BOOT     <- TRUE
DO_ENVELOPE <- TRUE   # Alterar para FALSE se necessário
B_BOOT      <- 1000   # Réplicas para suavização do resíduo
n_sim       <- 19     # Número de envelopes
B_ENV       <- 100    # Bootstrap interno do envelope

# ----------------------------------------------------------
# 2. CARREGAMENTO E PRÉ-PROCESSAMENTO DOS DADOS
# ----------------------------------------------------------
df <- read.csv("R6data.csv")
df <- R6data
df$Trt     <- as.factor(df$Trt)
df$Species <- as.factor(df$Species)

# Padronização de Day
df$Day <- as.numeric(scale(df$Day))

# Definição da composição (Simplex)
Yraw <- as.matrix(df[, c("LMF", "SMF", "RMF")])

eps_num <- 1e-6
Yraw[!is.finite(Yraw)] <- NA

vars_check <- c("Day", "lnTDM", "Trt", "Species")
ok <- complete.cases(Yraw, df[, vars_check])

df   <- df[ok, , drop = FALSE]
Yraw <- Yraw[ok, , drop = FALSE]

# fechamento + evita zeros
Yraw[Yraw <= 0] <- eps_num
Yraw <- Yraw / rowSums(Yraw)

df$y1 <- Yraw[, 1]
df$y2 <- Yraw[, 2]
df$y3 <- Yraw[, 3]
df$Y  <- DR_data(df[, c("y1", "y2", "y3")])

full_data <- df %>% dplyr::select(Species, Day, Trt, lnTDM, y1, y2, y3)

cat(sprintf("Dados processados. N total: %d observações.\n", nrow(full_data)))

# níveis globais (para não criar NA por fold)
LVL_SPECIES <- levels(full_data$Species)
LVL_TRT     <- levels(full_data$Trt)

# ----------------------------------------------------------
# 3. AJUSTE DO MODELO INFERENCIAL
# ----------------------------------------------------------
model_alt <- DirichReg(
  Y ~ I(Day^2) + lnTDM + Species * Trt * Day | Species * Day + Trt,
  data = df,
  model = "alternative"
)

cat("\n--- Resumo do Modelo Ajustado ---\n")
print(summary(model_alt))

# ----------------------------------------------------------
# 4. FUNÇÕES AUXILIARES: DIAGNÓSTICO E PREDIÇÃO
# ----------------------------------------------------------

# 4.1. Resíduos PIT Multivariado (Rosenblatt para Dirichlet)
pit_resid_dirichlet <- function(Y, alpha, eps = 1e-12) {
  Y <- as.matrix(Y); alpha <- as.matrix(alpha)
  n <- nrow(Y); k <- ncol(Y)
  U <- matrix(NA, nrow = n, ncol = k - 1)
  
  alpha_0 <- rowSums(alpha)
  U[, 1] <- pbeta(Y[, 1], shape1 = alpha[, 1], shape2 = alpha_0 - alpha[, 1])
  
  if (k > 2) {
    current_sum_Y <- Y[, 1]
    current_sum_alpha <- alpha_0 - alpha[, 1]
    for (j in 2:(k - 1)) {
      Y_cond <- Y[, j] / (1 - current_sum_Y)
      Y_cond <- pmin(pmax(Y_cond, 0), 1)
      shape2_new <- current_sum_alpha - alpha[, j]
      U[, j] <- pbeta(Y_cond, shape1 = alpha[, j], shape2 = shape2_new)
      current_sum_Y <- current_sum_Y + Y[, j]
      current_sum_alpha <- shape2_new
    }
  }
  pmin(pmax(U, eps), 1 - eps)
}

# 4.2. Colapso escalar (li)
collapse_li <- function(R) {
  Z <- qnorm(R)
  list(li = rowMeans(Z))
}

# 4.3. Quantil conformal
get_qhat <- function(scores, alpha) {
  sc <- sort(na.omit(scores))
  n <- length(sc)
  if (n <= 0) return(Inf)
  k <- ceiling((n + 1) * (1 - alpha))
  sc[min(k, n)]
}

# 4.3a Score quantílico (marginais Beta)
calc_score_quantile <- function(model, data) {
  if (nrow(data) == 0) return(numeric(0))
  
  pred <- predict(model, newdata = data, alpha = TRUE)
  alphas <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  
  Y <- as.matrix(data[, c("y1", "y2", "y3")])
  n <- nrow(Y); out <- numeric(n)
  
  for (i in 1:n) {
    a <- as.numeric(alphas[i, ])
    s0 <- sum(a)
    p <- pbeta(Y[i, ], shape1 = a, shape2 = s0 - a)
    p <- pmin(1 - 1e-8, pmax(1e-8, p))
    out[i] <- max(abs(qnorm(p)))
  }
  out
}

# 4.3b Score -loglik
calc_score_loglik <- function(model, data, eps = 1e-12) {
  if (nrow(data) == 0) return(numeric(0))
  pred <- predict(model, newdata = data, alpha = TRUE)
  A <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  
  Y <- as.matrix(data[, c("y1", "y2", "y3")])
  Y <- pmax(Y, eps)
  Y <- Y / rowSums(Y)
  Y <- pmin(pmax(Y, eps), 1 - eps)
  Y <- Y / rowSums(Y)
  
  logfy <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A - 1) * log(Y))
  -logfy
}

# ----------------------------------------------------------
# 4.4. tau : fechado + uniroot + NA->0 (igual simulação)
# ----------------------------------------------------------
.tstar_from_alpha_qhat <- function(alpha, qhat) {
  alpha <- as.numeric(alpha)
  qhat  <- as.numeric(qhat)[1]
  (sum(lgamma(alpha)) - lgamma(sum(alpha))) - qhat
}

.solve_rho_for_i <- function(w, W, tstar, i, tol = 1e-10) {
  w     <- as.numeric(w)
  W     <- as.numeric(W)[1]
  tstar <- as.numeric(tstar)[1]
  i     <- as.integer(i)[1]
  
  if (!isTRUE(is.finite(W)) || !isTRUE(is.finite(tstar)) || length(w) == 0L) return(NA_real_)
  if (i < 1L || i > length(w)) return(NA_real_)
  
  wi0 <- as.numeric(w[i])[1]
  if (!isTRUE(is.finite(wi0))) return(NA_real_)
  
  Sww <- as.numeric(sum(ifelse(w > 0, w * base::log(w), 0)))
  
  F <- function(rho) {
    rho <- as.numeric(rho)[1]
    if (!isTRUE(is.finite(rho)) || rho <= 0) return(-Inf)
    
    wi  <- as.numeric(w[i])[1]
    Wmi <- as.numeric(W - wi)[1]
    
    denom <- wi * rho + Wmi * (1 + rho)
    if (!isTRUE(is.finite(denom) && denom > 0)) return(-Inf)
    
    val <- wi * base::log(rho) +
      Wmi * base::log1p(rho) -
      W   * base::log(denom) +
      Sww - tstar
    
    if (!isTRUE(is.finite(val))) return(-Inf)
    val
  }
  
  lo <- 1e-12; hi <- 1e12
  flo <- F(lo); fhi <- F(hi)
  if (!isTRUE(is.finite(flo)) || !isTRUE(is.finite(fhi)) || flo * fhi > 0) return(NA_real_)
  
  ur <- try(stats::uniroot(F, c(lo, hi), tol = tol), silent = TRUE)
  if (inherits(ur, "try-error")) return(NA_real_)
  as.numeric(ur$root)
}

tau_from_alpha_qhat_closed <- function(alpha, qhat) {
  alpha <- as.numeric(alpha)
  qhat  <- as.numeric(qhat)[1]
  
  # mesma regra da simulação: se alpha < 1, não tenta fechado
  if (any(alpha < 1 - 1e-12) || any(!is.finite(alpha)))
    return(rep(NA_real_, length(alpha)))
  
  w <- alpha - 1
  W <- sum(w)
  tstar <- .tstar_from_alpha_qhat(alpha, qhat)
  
  tau <- numeric(length(alpha))
  for (i in seq_along(alpha)) {
    rho <- .solve_rho_for_i(w, W, tstar, i)
    if (!is.finite(rho)) {
      tau[i] <- NA_real_
    } else {
      denom <- w[i]/(1+rho) + (sum(w)-w[i])/rho
      if (!is.finite(denom) || denom <= 0) {
        tau[i] <- NA_real_
      } else {
        lambda <- 1 / denom
        tau[i] <- (lambda * w[i])/(1 + rho)
      }
    }
  }
  tau
}

# DROP-IN no lugar do tau antigo
tau_min_triangle <- function(alpha, qhat) {
  alpha <- as.numeric(alpha)
  qhat  <- as.numeric(qhat)[1]
  
  tau <- tau_from_alpha_qhat_closed(alpha, qhat)
  
  # FIEL AO ORIGINAL: NA -> 0
  tau[!is.finite(tau)] <- 0
  tau <- pmax(pmin(tau, 1), 0)
  
  s <- sum(tau)
  if (s > 1) tau <- tau / (s + 1e-12)
  tau
}

# ----------------------------------------------------------
# 4.5. build_triangle_grid (versão estável da simulação)
# ----------------------------------------------------------
build_triangle_grid <- function(tau, n1, n2, eps = 1e-12) {
  stopifnot(length(tau) == 3)
  
  tau <- as.numeric(tau)
  tau[!is.finite(tau)] <- 0
  tau <- pmax(pmin(tau, 1), 0)
  
  S <- sum(tau)
  if (S > 1 - eps) {
    y <- if (S > 0) tau / S else c(1/3, 1/3, 1/3)
    return(data.frame(y1 = y[1], y2 = y[2], y3 = y[3]))
  }
  
  y1_min <- tau[1]
  y1_max <- 1 - tau[2] - tau[3]
  
  if (!is.finite(y1_min) || !is.finite(y1_max) || y1_max < y1_min + eps) {
    y <- c(tau[1], tau[2], 1 - tau[1] - tau[2])
    y <- pmax(pmin(y, 1), 0)
    y <- y / sum(y)
    return(data.frame(y1 = y[1], y2 = y[2], y3 = y[3]))
  }
  
  y1_seq <- seq(y1_min, y1_max, length.out = n1)
  
  out_list <- vector("list", length(y1_seq))
  for (k in seq_along(y1_seq)) {
    y1k <- y1_seq[k]
    y2_min <- tau[2]
    y2_max <- 1 - tau[3] - y1k
    
    if (!is.finite(y2_max) || y2_max < y2_min + eps) next
    
    y2_seq <- seq(y2_min, y2_max, length.out = n2)
    y3_seq <- 1 - y1k - y2_seq
    
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
  
  if (is.null(grid) || nrow(grid) == 0) {
    y <- c(tau[1], tau[2], 1 - tau[1] - tau[2])
    y <- pmax(pmin(y, 1), 0)
    y <- y / sum(y)
    grid <- data.frame(y1 = y[1], y2 = y[2], y3 = y[3])
  }
  
  grid
}

# (sem mudança na sua avaliação do HDR_Grid)
is_inside_valid_region <- function(valid_df, y_test, eps = 1e-12) {
  if (is.null(valid_df) || nrow(valid_df) == 0) return(FALSE)
  y_test <- as.numeric(y_test)
  if (length(y_test) != 3 || any(!is.finite(y_test))) return(FALSE)
  if (abs(sum(y_test) - 1) > 1e-8) return(FALSE)
  
  y1_vals <- sort(unique(valid_df$y1))
  y1t <- y_test[1]; y2t <- y_test[2]
  if (y1t < min(y1_vals) - eps || y1t > max(y1_vals) + eps) return(FALSE)
  
  y2_min <- numeric(length(y1_vals)); y2_max <- numeric(length(y1_vals))
  for (k in seq_along(y1_vals)) {
    sl <- valid_df[abs(valid_df$y1 - y1_vals[k]) <= 1e-15, , drop = FALSE]
    y2_min[k] <- min(sl$y2); y2_max[k] <- max(sl$y2)
  }
  
  k_eq <- which(abs(y1_vals - y1t) <= 1e-12)
  if (length(k_eq) > 0) {
    k <- k_eq[1]
    return(y2t >= y2_min[k] - eps && y2t <= y2_max[k] + eps)
  }
  k_hi <- which(y1_vals > y1t)[1]
  if (is.na(k_hi) || k_hi <= 1) return(FALSE)
  
  k_lo <- k_hi - 1
  w <- (y1t - y1_vals[k_lo]) / (y1_vals[k_hi] - y1_vals[k_lo])
  min_i <- (1 - w) * y2_min[k_lo] + w * y2_min[k_hi]
  max_i <- (1 - w) * y2_max[k_lo] + w * y2_max[k_hi]
  
  (y2t >= min_i - eps && y2t <= max_i + eps)
}

# ----------------------------------------------------------
# 5. DIAGNÓSTICO DO MODELO (RESÍDUOS E ENVELOPE)
# ----------------------------------------------------------
if (DO_BOOT) {
  cat("\n--- Iniciando Diagnóstico de Resíduos ---\n")
  
  # CONSISTENTE com o ajuste: usa y1:y3 (não LMF:RMF)
  Y_obs <- as.matrix(df[, c("y1", "y2", "y3")])
  alpha_hat <- as.matrix(model_alt$fitted.values$alpha)
  n <- nrow(Y_obs)
  
  # A. resíduo observado (suavizado)
  R0  <- pit_resid_dirichlet(Y_obs, alpha_hat)
  li0 <- collapse_li(R0)$li
  
  cat("Gerando réplicas bootstrap para suavização...\n")
  pb <- progress::progress_bar$new(total = B_BOOT, format = " [:bar] :percent :eta")
  li_boot <- matrix(NA_real_, n, B_BOOT)
  
  for (b in 1:B_BOOT) {
    pb$tick()
    Yb <- DirichletReg::rdirichlet(n, alpha_hat)
    Yb <- pmax(Yb, 1e-12); Yb <- Yb / rowSums(Yb)
    
    df_b <- df; df_b$Y <- DR_data(Yb)
    
    fit_b <- tryCatch({
      DirichReg(Y ~ I(Day^2) + lnTDM + Species * Trt * Day |
                  Species * Day + Trt,
                data = df_b, model = "alternative",
                control = list(iterlim = 5000, tol1 = 1e-7))
    }, error = function(e) NULL)
    
    if (!is.null(fit_b)) {
      alpha_b <- as.matrix(fit_b$fitted.values$alpha)
      Rb <- pit_resid_dirichlet(Yb, alpha_b)
      li_boot[, b] <- collapse_li(Rb)$li
    }
  }
  
  # Quantil-normal
  ai <- integer(n); si <- numeric(n)
  li_boot_clean <- li_boot[, colSums(is.na(li_boot)) == 0, drop = FALSE]
  B_EFF <- ncol(li_boot_clean)
  
  for (i in 1:n) {
    ai[i] <- sum(li_boot_clean[i, ] < li0[i])
    si[i] <- qnorm(runif(1, ai[i] / (B_EFF + 1), (ai[i] + 1) / (B_EFF + 1)))
  }
  
  # Envelope
  s_ord <- sort(si)
  x_env <- qnorm((seq_len(n) - 3/8) / (n + 1/4))
  yrange <- range(s_ord, na.rm = TRUE)
  
  if (DO_ENVELOPE) {
    cat("\nGerando Envelope de Confiança (Simulação sob H0)...\n")
    pb2 <- progress::progress_bar$new(total = n_sim, format = " [:bar] :percent :eta")
    res_sim <- matrix(NA_real_, n, n_sim)
    
    for (e in 1:n_sim) {
      pb2$tick()
      Ye <- DirichletReg::rdirichlet(n, alpha_hat)
      Ye <- pmax(Ye, 1e-12); Ye <- Ye / rowSums(Ye)
      
      df_e <- df; df_e$Y <- DR_data(Ye)
      
      fit_e <- tryCatch({
        DirichReg(Y ~ I(Day^2) + lnTDM + Species * Trt * Day |
                    Species * Day + Trt,
                  data = df_e, model = "alternative",
                  control = list(iterlim = 5000, tol1 = 1e-7))
      }, error = function(x) NULL)
      
      if (is.null(fit_e)) next
      
      alpha_e <- as.matrix(fit_e$fitted.values$alpha)
      Re  <- pit_resid_dirichlet(Ye, alpha_e)
      lie <- collapse_li(Re)$li
      
      # Bootstrap aninhado
      lie_boot_in <- matrix(NA_real_, n, B_ENV)
      for (b in 1:B_ENV) {
        Yeb <- DirichletReg::rdirichlet(n, alpha_e)
        Yeb <- pmax(Yeb, 1e-12); Yeb <- Yeb / rowSums(Yeb)
        df_eb <- df_e; df_eb$Y <- DR_data(Yeb)
        
        fit_eb <- tryCatch({
          DirichReg(Y ~ I(Day^2) + lnTDM + Species * Trt * Day |
                      Species * Day + Trt,
                    data = df_eb, model = "alternative",
                    control = list(iterlim = 2000, tol1 = 1e-6))
        }, error = function(x) NULL)
        
        if (!is.null(fit_eb)) {
          lie_boot_in[, b] <- collapse_li(
            pit_resid_dirichlet(Yeb, as.matrix(fit_eb$fitted.values$alpha))
          )$li
        }
      }
      
      valid_in <- colSums(is.na(lie_boot_in)) == 0
      lie_clean <- lie_boot_in[, valid_in, drop = FALSE]
      
      if (ncol(lie_clean) > 10) {
        sie <- numeric(n); B_in <- ncol(lie_clean)
        for (i in 1:n) {
          aie <- sum(lie_clean[i, ] < lie[i])
          sie[i] <- qnorm(runif(1, aie / (B_in + 1), (aie + 1) / (B_in + 1)))
        }
        res_sim[, e] <- sie
      }
    }
    
    res_sim <- res_sim[, colSums(is.na(res_sim)) == 0, drop = FALSE]
    if (ncol(res_sim) > 2) {
      s_all_ord <- apply(cbind(si, res_sim), 2, sort)
      min_env <- apply(s_all_ord[, -1, drop = FALSE], 1, min)
      max_env <- apply(s_all_ord[, -1, drop = FALSE], 1, max)
      yrange <- range(s_all_ord, na.rm = TRUE)
      
      par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
      plot(x_env, s_ord, pch = 16, cex = 0.6, ylim = yrange,
           xlab = "Quantis Teóricos N(0,1)", ylab = "Resíduos do Modelo",
           main = "Diagnóstico: QQ-Plot com Envelope")
      polygon(c(x_env, rev(x_env)), c(min_env, rev(max_env)),
              col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
      lines(x_env, min_env, lty = 2, col = "gray50")
      lines(x_env, max_env, lty = 2, col = "gray50")
      abline(0, 1, col = "red", lwd = 1.5)
      points(x_env, s_ord, pch = 16, cex = 0.6)
    }
  } else {
    par(mfrow = c(1, 1))
    plot(x_env, s_ord, pch = 16, cex = 0.6,
         xlab = "Quantis Teóricos N(0,1)",
         ylab = "Resíduos do Modelo", main = "QQ-Plot de Resíduos")
    abline(0, 1, col = "red", lwd = 1.5, lty = 2)
  }
}

# ------------------------------------------------------------------------------
# 6. VALIDAÇÃO PREDITIVA (CONFORMAL PREDICTION - CV)
# ------------------------------------------------------------------------------
cat(paste0("\n--- Iniciando Validação Cruzada (K=", K_FOLDS, ") ---\n"))

N_TOTAL <- nrow(full_data)
fold_ids <- sample(rep(1:K_FOLDS, length.out = N_TOTAL))

# Wrappers para os métodos (MANTÉM A FORMA DE AVALIAÇÃO)
run_quantile <- function(model, cal, te, alpha) {
  qhat <- get_qhat(calc_score_quantile(model, cal), alpha)
  pred <- predict(model, newdata = te, alpha = TRUE)
  a <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  s0 <- rowSums(a)
  
  L <- qbeta(pnorm(-qhat), a, s0 - a)
  U <- qbeta(pnorm(qhat), a, s0 - a)
  y <- as.numeric(te[, c("y1", "y2", "y3")])
  
  covered <- (y[1] >= L[1] && y[1] <= U[1] &&
                y[2] >= L[2] && y[2] <= U[2] &&
                y[3] >= L[3] && y[3] <= U[3])
  
  data.frame(method = "Quantile", cov_global = covered,
             cov_y1 = (y[1] >= L[1] && y[1] <= U[1]),
             cov_y2 = (y[2] >= L[2] && y[2] <= U[2]),
             cov_y3 = (y[3] >= L[3] && y[3] <= U[3]),
             w_y1 = U[1] - L[1], w_y2 = U[2] - L[2], w_y3 = U[3] - L[3])
}

run_hdr_split <- function(model, cal, te, alpha) {
  qhat <- get_qhat(calc_score_loglik(model, cal), alpha)
  pred <- predict(model, newdata = te, alpha = TRUE)
  tau <- tau_min_triangle(if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred, qhat)
  
  L <- tau
  U <- pmax(0, 1 - c(sum(tau[-1]), sum(tau[-2]), sum(tau[-3])))
  y <- as.numeric(te[, c("y1", "y2", "y3")])
  
  covered <- (y[1] >= L[1] && y[1] <= U[1] &&
                y[2] >= L[2] && y[2] <= U[2] &&
                y[3] >= L[3] && y[3] <= U[3])
  
  data.frame(method = "HDR_Split", cov_global = covered,
             cov_y1 = (y[1] >= L[1] && y[1] <= U[1]),
             cov_y2 = (y[2] >= L[2] && y[2] <= U[2]),
             cov_y3 = (y[3] >= L[3] && y[3] <= U[3]),
             w_y1 = U[1] - L[1], w_y2 = U[2] - L[2], w_y3 = U[3] - L[3])
}

run_hdr_grid <- function(model, cal, te, alpha) {
  qhat <- get_qhat(calc_score_loglik(model, cal), alpha)
  pred <- predict(model, newdata = te, alpha = TRUE)
  tau <- tau_min_triangle(if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred, qhat)
  
  grid_df <- build_triangle_grid(tau, n1 = GRID_N1, n2 = GRID_N2)
  w <- c(NA, NA, NA); c_glob <- FALSE; c_marg <- rep(FALSE, 3)
  
  if (nrow(grid_df) > 0) {
    covs <- te[rep(1, nrow(grid_df)), setdiff(names(te), c("y1","y2","y3","Y")), drop = FALSE]
    newdata <- cbind(covs, grid_df)
    s_grid <- calc_score_loglik(model, newdata)
    valid <- grid_df[s_grid <= qhat, , drop = FALSE]
    
    if (nrow(valid) > 0) {
      y <- as.numeric(te[, c("y1", "y2", "y3")])
      c_glob <- is_inside_valid_region(valid, y)
      rg <- apply(valid, 2, range)
      w <- rg[2, ] - rg[1, ]
      c_marg <- (y >= rg[1, ] & y <= rg[2, ])
    }
  }
  
  data.frame(method = "HDR_Grid", cov_global = c_glob,
             cov_y1 = c_marg[1], cov_y2 = c_marg[2], cov_y3 = c_marg[3],
             w_y1 = w[1], w_y2 = w[2], w_y3 = w[3])
}

# Loop de Cross-Validation
results_list <- list()
pb <- progress::progress_bar$new(total = K_FOLDS,
                                 format = " Fold :current/:total [:bar] :percent :eta")

for (k in 1:K_FOLDS) {
  
  test_idx <- which(fold_ids == k)
  rest_idx <- setdiff(seq_len(N_TOTAL), test_idx)
  
  n_tr <- floor(0.7 * N_TOTAL)
  set.seed(13 + k)
  train_idx <- sample(rest_idx, n_tr)
  calib_idx <- setdiff(rest_idx, train_idx)
  
  dat_tr <- full_data[train_idx, , drop = FALSE]
  dat_ca <- full_data[calib_idx, , drop = FALSE]
  dat_te <- full_data[test_idx,  , drop = FALSE]
  
  # fixa níveis globais (evita NA por fold)
  dat_tr$Species <- factor(dat_tr$Species, levels = LVL_SPECIES)
  dat_ca$Species <- factor(dat_ca$Species, levels = LVL_SPECIES)
  dat_te$Species <- factor(dat_te$Species, levels = LVL_SPECIES)
  
  dat_tr$Trt <- factor(dat_tr$Trt, levels = LVL_TRT)
  dat_ca$Trt <- factor(dat_ca$Trt, levels = LVL_TRT)
  dat_te$Trt <- factor(dat_te$Trt, levels = LVL_TRT)
  
  dat_tr$Y <- DR_data(dat_tr[, c("y1", "y2", "y3")])
  dat_ca$Y <- DR_data(dat_ca[, c("y1", "y2", "y3")])
  
  model <- tryCatch({
    DirichReg(Y ~ I(Day^2) + lnTDM + Species * Trt * Day |
                Species * Day + Trt,
              data = dat_tr, model = "alternative")
  }, error = function(e) NULL)
  
  if (is.null(model)) { pb$tick(); next }
  
  fold_res <- vector("list", nrow(dat_te))
  for (i in 1:nrow(dat_te)) {
    pt <- dat_te[i, , drop = FALSE]
    # mantém seu critério original (não altera avaliação)
    if (any(is.na(pt$Species)) || any(is.na(pt$Trt))) next
    
    r1 <- run_quantile(model, dat_ca, pt, ALPHA_CONFORMAL)
    r2 <- run_hdr_split(model, dat_ca, pt, ALPHA_CONFORMAL)
    r3 <- run_hdr_grid(model, dat_ca, pt, ALPHA_CONFORMAL)
    fold_res[[i]] <- dplyr::bind_rows(r1, r2, r3)
  }
  
  results_list[[k]] <- dplyr::bind_rows(fold_res)
  pb$tick()
}

# ----------------------------------------------------------
# 7. RESULTADOS FINAIS
# ----------------------------------------------------------
final_df <- dplyr::bind_rows(results_list)

summary_tab <- final_df %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    Global_Cov = mean(cov_global, na.rm = TRUE),
    Cov_Y1 = mean(cov_y1, na.rm = TRUE),
    Cov_Y2 = mean(cov_y2, na.rm = TRUE),
    Cov_Y3 = mean(cov_y3, na.rm = TRUE),
    Width_Y1 = mean(w_y1, na.rm = TRUE),
    Width_Y2 = mean(w_y2, na.rm = TRUE),
    Width_Y3 = mean(w_y3, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=======================================================\n")
cat(sprintf("RESULTADOS FINAIS (R6data) - %d-Fold CV\n", K_FOLDS))
cat(sprintf("Nível de Cobertura Alvo: %.1f%%\n", (1 - ALPHA_CONFORMAL) * 100))
cat("=======================================================\n")
print(summary_tab)
