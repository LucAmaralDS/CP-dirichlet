################################################################################
# (MODELO TST | TST)
################################################################################

rm(list = ls())

# 1. PACOTES
pkgs <- c("DirichletReg", "tidyverse", "readxl", "progress", "dplyr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(DirichletReg)
  library(tidyverse)
  library(readxl)
  library(progress)
})

set.seed(13)

# ------------------------------------------------------------------------------
# 2. CARREGAMENTO E PRÉ-PROCESSAMENTO
# ------------------------------------------------------------------------------
# Ajuste o caminho se necessário
dados_boot <- read_excel("C:/Users/Lucas/Desktop/CP-dirichlet/aplicacao/sleep/dados-boot.xlsx")

if (!exists("dados_boot")) stop("Objeto dados_boot não encontrado.")

df <- dados_boot

# Tratamento de Covariáveis
df$CaseControl <- as.factor(df$CaseControl)
df$TST  <- as.numeric(df$TST) 

# Definição da composição
Yraw <- as.matrix(df[, c("N1_PERC", "N2_PERC", "N3_PERC", "REM_PERC")])
eps_num <- 1e-6
Yraw[!is.finite(Yraw)] <- NA
vars_check <- c("CaseControl", "TST") 
ok <- complete.cases(Yraw, df[, vars_check])

df   <- df[ok, , drop = FALSE]
Yraw <- Yraw[ok, , drop = FALSE]
Yraw[Yraw <= 0] <- eps_num
Yraw <- Yraw / rowSums(Yraw)

df$y1 <- Yraw[, 1]; df$y2 <- Yraw[, 2]; df$y3 <- Yraw[, 3]; df$y4 <- Yraw[, 4]
df$Y <- DR_data(df[, c("y1", "y2", "y3", "y4")])

cat(sprintf("Dataset Original: N = %d observações.\n", nrow(df)))

# -----------------------------------------
# 3. DATA AUGMENTATION 
# -----------------------------------------

# Ajuste do Modelo Gerador
fit_gen <- DirichReg(
  Y ~ TST|TST ,
  data  = df,
  model = "alternative"
)
fit_gen

# --- EXTRAÇÃO DOS PARÂMETROS ---
# Mu (Médias) para cada observação
mu_hat <- fitted(fit_gen) 

# Phi (Precisão) para cada observação
# 
coefs_phi <- fit_gen$coefficients[7:8]
phi_hat_vec <- as.numeric(exp(cbind(1, df$TST) %*% coefs_phi))

# Definição do Tamanho
N_TARGET <- 1000
n_orig   <- nrow(df)
n_new    <- N_TARGET - n_orig

if (n_new > 0) {
  cat(sprintf("Gerando %d novas observações sintéticas...\n", n_new))
  
  # 1. Bootstrap das Covariáveis (índices)
  idx_rep <- sample(seq_len(n_orig), size = n_new, replace = TRUE)
  
  # 2. Cria dataframe novo replicando as covariáveis
  df_new  <- df[idx_rep, , drop = FALSE]
  
  # 3. Recupera os parâmetros (Mu e Phi) associados a esses índices
  mu_new  <- mu_hat[idx_rep, ]
  phi_new <- phi_hat_vec[idx_rep] # Vetor de precisão correspondente
  
  # 4. Gera Y via Dirichlet
  # (phi) *(mu)
  Y_new_raw <- t(sapply(1:n_new, function(i) {
    DirichletReg::rdirichlet(1, alpha = phi_new[i] * mu_new[i, ])
  }))
  
  # 5. Organiza df_new
  df_new$y1 <- Y_new_raw[, 1]; df_new$y2 <- Y_new_raw[, 2]
  df_new$y3 <- Y_new_raw[, 3]; df_new$y4 <- Y_new_raw[, 4]
  df_new$Y <- DR_data(df_new[, c("y1", "y2", "y3", "y4")])
  
  # 6. Junta Tudo
  full_data_raw <- rbind(df, df_new)
} else {
  full_data_raw <- df
}

full_data <- full_data_raw %>% dplyr::select(CaseControl, TST, y1, y2, y3, y4)
full_data$Y <- DR_data(full_data[, c("y1", "y2", "y3", "y4")])
LVL_CASE <- levels(full_data$CaseControl)

cat(sprintf("Dataset Final Pronto: N = %d\n", nrow(full_data)))

summary(full_data)
# ------------------------------------------------------------------------------
# 4. FUNÇÕES E CONFIGURAÇÕES 
# ------------------------------------------------------------------------------
calc_score_quantile <- function(model, data) {
  if (nrow(data) == 0) return(numeric(0))
  pred <- predict(model, newdata = data, alpha = TRUE)
  alphas <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  n <- nrow(Y); out <- numeric(n)
  for (i in 1:n) {
    a <- alphas[i, ]; s0 <- sum(a)
    res_j <- numeric(length(a))
    for(j in seq_along(a)) {
      p <- pbeta(Y[i, j], shape1 = a[j], shape2 = s0 - a[j])
      res_j[j] <- qnorm(pmin(1-1e-8, pmax(1e-8, p)))
    }
    out[i] <- max(abs(res_j))
  }
  out
}

calc_score_loglik <- function(model, data, eps = 1e-12) {
  if (nrow(data) == 0) return(numeric(0))
  pred <- predict(model, newdata = data, alpha = TRUE)
  A <- if (is.list(pred) && !is.null(pred$alpha)) pred$alpha else pred
  Y <- as.matrix(data[, grep("^y[0-9]+$", names(data))])
  Y <- pmax(Y, eps); Y <- Y / rowSums(Y)
  logfy <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A - 1) * log(Y))
  -logfy
}

get_qhat <- function(scores, alpha) {
  sc <- sort(na.omit(scores)); n <- length(sc)
  if (n <= 0) return(Inf)
  k <- ceiling((n + 1) * (1 - alpha))
  sc[min(k, n)]
}

# Funções Tau e Grid (Mantidas iguais)
.tstar_from_alpha_qhat <- function(alpha, qhat) { (sum(lgamma(alpha)) - lgamma(sum(alpha))) - qhat }
.solve_rho_for_i <- function(w, W, tstar, i) {
  wi0 <- w[i]; Sww <- sum(ifelse(w > 0, w * log(w), 0))
  F_obj <- function(rho) {
    if(rho <= 0) return(-Inf)
    Wmi <- W - wi0; denom <- wi0 * rho + Wmi * (1 + rho)
    if(denom <= 0) return(-Inf)
    wi0 * log(rho) + Wmi * log1p(rho) - W * log(denom) + Sww - tstar
  }
  tryCatch(uniroot(F_obj, c(1e-12, 1e12))$root, error = function(e) NA_real_)
}
tau_min_simplex <- function(alpha, qhat) {
  alpha <- as.numeric(alpha); qhat <- as.numeric(qhat)
  if(any(alpha < 1)) return(rep(0, length(alpha))) 
  w <- alpha - 1; W <- sum(w); tstar <- .tstar_from_alpha_qhat(alpha, qhat)
  tau <- numeric(length(alpha))
  for (i in seq_along(alpha)) {
    rho <- .solve_rho_for_i(w, W, tstar, i)
    if (is.na(rho)) tau[i] <- 0
    else {
      denom <- w[i]/(1+rho) + (W-w[i])/rho
      tau[i] <- (1/denom) * w[i] / (1 + rho)
    }
  }
  tau <- pmax(pmin(tau, 1), 0)
  if(sum(tau) > 1) tau <- tau / (sum(tau) + 1e-12)
  tau
}

build_polytope_grid_4d <- function(tau, res=25) {
  tau <- as.numeric(tau); tau[!is.finite(tau)] <- 0
  if(sum(tau) > 0.99) { pt <- tau/sum(tau); return(data.frame(y1=pt[1], y2=pt[2], y3=pt[3], y4=pt[4])) }
  max_y1 <- 1 - sum(tau[2:4]); if(max_y1 < tau[1]) max_y1 <- tau[1] + 1e-5
  seq_y1 <- seq(tau[1], max_y1, length.out=res)
  grid_list <- list()
  for(val_y1 in seq_y1) {
    max_y2 <- 1 - val_y1 - sum(tau[3:4]); if(max_y2 < tau[2]) next
    seq_y2 <- seq(tau[2], max_y2, length.out=res)
    for(val_y2 in seq_y2) {
      max_y3 <- 1 - val_y1 - val_y2 - tau[4]; if(max_y3 < tau[3]) next
      seq_y3 <- seq(tau[3], max_y3, length.out=res)
      val_y4 <- 1 - val_y1 - val_y2 - seq_y3
      valid <- val_y4 >= tau[4] - 1e-9
      if(any(valid)) grid_list[[length(grid_list)+1]] <- data.frame(y1=val_y1, y2=val_y2, y3=seq_y3[valid], y4=val_y4[valid])
    }
  }
  do.call(rbind, grid_list)
}

# Wrappers
run_quantile <- function(model, cal, te_row, alpha) {
  qhat <- get_qhat(calc_score_quantile(model, cal), alpha)
  pred <- predict(model, newdata = te_row, alpha = TRUE)
  a <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  L <- qbeta(pnorm(-qhat), a, sum(a) - a)
  U <- qbeta(pnorm(qhat), a, sum(a) - a)
  y <- as.numeric(te_row[, c("y1","y2","y3","y4")])
  cov_marg <- (y >= L & y <= U)
  data.frame(method="Quantile", cov_global=all(cov_marg),
             cov_y1=cov_marg[1], cov_y2=cov_marg[2], cov_y3=cov_marg[3], cov_y4=cov_marg[4],
             w_y1=U[1]-L[1], w_y2=U[2]-L[2], w_y3=U[3]-L[3], w_y4=U[4]-L[4])
}

run_hdr_split <- function(model, cal, te_row, alpha) {
  qhat <- get_qhat(calc_score_loglik(model, cal), alpha)
  pred <- predict(model, newdata = te_row, alpha = TRUE)
  a <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  tau <- tau_min_simplex(a, qhat)
  L <- tau; U <- numeric(4)
  for(k in 1:4) U[k] <- pmax(0, 1 - sum(tau[-k]))
  y <- as.numeric(te_row[, c("y1","y2","y3","y4")])
  cov_marg <- (y >= L & y <= U)
  data.frame(method="HDR_Split", cov_global=all(cov_marg),
             cov_y1=cov_marg[1], cov_y2=cov_marg[2], cov_y3=cov_marg[3], cov_y4=cov_marg[4],
             w_y1=U[1]-L[1], w_y2=U[2]-L[2], w_y3=U[3]-L[3], w_y4=U[4]-L[4])
}

run_hdr_grid <- function(model, cal, te_row, alpha) {
  qhat <- get_qhat(calc_score_loglik(model, cal), alpha)
  pred <- predict(model, newdata = te_row, alpha = TRUE)
  a <- as.numeric(if (is.list(pred)) pred$alpha else pred)
  tau <- tau_min_simplex(a, qhat)
  grid_df <- build_polytope_grid_4d(tau, res=25)
  c_glob <- FALSE; c_marg <- rep(FALSE, 4); w <- rep(NA, 4)
  if(!is.null(grid_df) && nrow(grid_df) > 0) {
    covs <- te_row[rep(1, nrow(grid_df)), setdiff(names(te_row), c("y1","y2","y3","y4","Y")), drop=FALSE]
    new_dat <- cbind(covs, grid_df)
    new_dat$Y <- DR_data(new_dat[, c("y1","y2","y3","y4")])
    s_grid <- calc_score_loglik(model, new_dat)
    valid <- grid_df[s_grid <= qhat, , drop=FALSE]
    if(nrow(valid) > 0) {
      y_real <- as.numeric(te_row[, c("y1","y2","y3","y4")])
      rg <- apply(valid, 2, range)
      w <- rg[2,] - rg[1,]
      c_marg <- (y_real >= rg[1,] & y_real <= rg[2,])
      c_glob <- all(c_marg)
    }
  }
  data.frame(method="HDR_Grid", cov_global=c_glob,
             cov_y1=c_marg[1], cov_y2=c_marg[2], cov_y3=c_marg[3], cov_y4=c_marg[4],
             w_y1=w[1], w_y2=w[2], w_y3=w[3], w_y4=w[4])
}

# ------------------------------------------------------------------------------
# 5. LOOP CV (10-FOLD)
# ------------------------------------------------------------------------------
ALPHA_CONFORMAL <- 0.10
K_FOLDS  <- 10
cat(paste0("\n--- Iniciando Validação Cruzada (K=", K_FOLDS, ") ---\n"))
N_TOTAL <- nrow(full_data)
fold_ids <- sample(rep(1:K_FOLDS, length.out = N_TOTAL))
results_list <- list()
pb <- progress::progress_bar$new(total = K_FOLDS, format = " [:bar] :percent :eta")

for (k in 1:K_FOLDS) {
  test_idx <- which(fold_ids == k)
  rest_idx <- setdiff(seq_len(N_TOTAL), test_idx)
  n_train_target <- floor(0.70 * N_TOTAL)
  
  set.seed(100 + k)
  train_idx <- sample(rest_idx, n_train_target)
  calib_idx <- setdiff(rest_idx, train_idx)
  
  dat_tr <- full_data[train_idx, , drop=FALSE]
  dat_ca <- full_data[calib_idx, , drop=FALSE]
  dat_te <- full_data[test_idx,  , drop=FALSE]
  
  # Recria objetos DR_data
  dat_tr$Y <- DR_data(dat_tr[, c("y1","y2","y3","y4")])
  dat_ca$Y <- DR_data(dat_ca[, c("y1","y2","y3","y4")])
  
  # Ajuste do Modelo com TST_scaled
  model <- tryCatch({
    DirichReg(Y ~ TST|TST, data = dat_tr, model = "alternative")
  }, error = function(e) NULL)
  
  if (is.null(model)) { pb$tick(); next }
  
  fold_res <- list()
  for(i in 1:nrow(dat_te)) {
    pt <- dat_te[i, , drop=FALSE]
    r1 <- run_quantile(model, dat_ca, pt, ALPHA_CONFORMAL)
    r2 <- run_hdr_split(model, dat_ca, pt, ALPHA_CONFORMAL)
    r3 <- run_hdr_grid(model, dat_ca, pt, ALPHA_CONFORMAL)
    fold_res[[i]] <- bind_rows(r1, r2, r3)
  }
  results_list[[k]] <- bind_rows(fold_res)
  pb$tick()
}

# 6. RESULTADOS
final_df <- bind_rows(results_list)
summary_tab <- final_df %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    Global_Cov = mean(cov_global, na.rm = TRUE),
    Cov_N1 = mean(cov_y1, na.rm=TRUE), Width_N1 = mean(w_y1, na.rm=TRUE),
    Cov_N2 = mean(cov_y2, na.rm=TRUE), Width_N2 = mean(w_y2, na.rm=TRUE),
    Cov_N3 = mean(cov_y3, na.rm=TRUE), Width_N3 = mean(w_y3, na.rm=TRUE),
    Cov_REM = mean(cov_y4, na.rm=TRUE), Width_REM = mean(w_y4, na.rm=TRUE)
  )

cat("\nRESULTADOS FINAIS\n")
print(summary_tab)



# =======================================
# 6. DIAGNÓSTICO DO MODELO 
# =======================================
# 4.1 Funções Auxiliares de Diagnóstico
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

collapse_li <- function(R) {
  Z <- qnorm(R)
  list(li = rowMeans(Z))
}

# 4.2 Execução do Diagnóstico
DO_DIAGNOSTICS <- TRUE
B_BOOT <- 1000  # Aumente para 1000 para publicação final
n_sim  <- 19
B_ENV  <- 1000   # Bootstrap interno

if (DO_DIAGNOSTICS) {
  cat("\n--- Iniciando Diagnóstico de Ajuste (Dados Aumentados) ---\n")
  
  # Ajuste do modelo no dataset completo
  model_diag <- DirichReg(Y ~ TST | TST, data = full_data, model = "alternative")
  
  Y_obs <- as.matrix(full_data[, c("y1", "y2", "y3", "y4")])
  alpha_hat <- as.matrix(model_diag$fitted.values$alpha)
  n <- nrow(Y_obs)
  
  # Resíduo Observado
  R0 <- pit_resid_dirichlet(Y_obs, alpha_hat)
  li0 <- collapse_li(R0)$li
  
  cat("Gerando Envelope Simulado (aguarde)...\n")
  pb <- progress::progress_bar$new(total = n_sim, format = " [:bar] :percent :eta")
  
  # Matriz para guardar os quantis simulados
  res_sim <- matrix(NA_real_, n, n_sim)
  
  # Loop do Envelope
  for (e in 1:n_sim) {
    pb$tick()
    # 1. Simular Y sob H0 (usando modelo ajustado)
    Ye <- DirichletReg::rdirichlet(n, alpha_hat)
    Ye <- pmax(Ye, 1e-12); Ye <- Ye / rowSums(Ye)
    
    df_e <- full_data; df_e$Y <- DR_data(Ye)
    
    # 2. Reajustar modelo
    fit_e <- tryCatch({
      DirichReg(Y ~ TST|TST, data = df_e, model = "alternative", control=list(iterlim=1000))
    }, error = function(x) NULL)
    
    if (is.null(fit_e)) next
    
    alpha_e <- as.matrix(fit_e$fitted.values$alpha)
    Re <- pit_resid_dirichlet(Ye, alpha_e)
    lie <- collapse_li(Re)$li
    
    # 3. Bootstrap Interno 
    lie_boot_in <- matrix(NA_real_, n, B_ENV)
    for (b in 1:B_ENV) {
      Yeb <- DirichletReg::rdirichlet(n, alpha_e)
      Yeb <- pmax(Yeb, 1e-12); Yeb <- Yeb / rowSums(Yeb)
      df_eb <- df_e; df_eb$Y <- DR_data(Yeb)
      
      fit_eb <- tryCatch({
        DirichReg(Y ~  TST|TST, data = df_eb, model = "alternative", control=list(iterlim=1000))
      }, error = function(x) NULL)
      
      if (!is.null(fit_eb)) {
        lie_boot_in[, b] <- collapse_li(pit_resid_dirichlet(Yeb, as.matrix(fit_eb$fitted.values$alpha)))$li
      }
    }
    
    # Calcula posição relativa 
    valid_in <- colSums(is.na(lie_boot_in)) == 0
    lie_clean <- lie_boot_in[, valid_in, drop = FALSE]
    if (ncol(lie_clean) > 5) {
      sie <- numeric(n); B_in_eff <- ncol(lie_clean)
      for (i in 1:n) {
        aie <- sum(lie_clean[i, ] < lie[i])
        sie[i] <- qnorm(runif(1, aie / (B_in_eff + 1), (aie + 1) / (B_in_eff + 1)))
      }
      res_sim[, e] <- sie
    }
  }
  
  # Processamento final para o plot
  
  ai <- rank(li0) - 1 # Aproximação rápida
  si <- qnorm(runif(n, ai/(n+1), (ai+1)/(n+1))) 
  
  # Ordenar
  s_ord <- sort(si)
  res_sim_ord <- apply(res_sim, 2, sort, na.last = NA)
  
  if (ncol(res_sim_ord) > 0) {
    min_env <- apply(res_sim_ord, 1, min, na.rm=TRUE)
    max_env <- apply(res_sim_ord, 1, max, na.rm=TRUE)
    x_env <- qnorm((1:n - 0.375) / (n + 0.25))
    
    # Plot
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
    yrange <- range(c(s_ord, min_env, max_env), na.rm=TRUE)
    plot(x_env, s_ord, pch = 16, cex = 0.5, ylim = yrange,
         xlab = "Quantis Teóricos N(0,1)", ylab = "Resíduos do Modelo",
         main = "Diagnóstico: QQ-Plot com Envelope (Simulado)")
    polygon(c(x_env, rev(x_env)), c(min_env, rev(max_env)),
            col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
    lines(x_env, min_env, lty = 2, col = "gray50")
    lines(x_env, max_env, lty = 2, col = "gray50")
    abline(0, 1, col = "red", lwd = 1.5)
    points(x_env, s_ord, pch = 16, cex = 0.5)
  }
  cat("\nDiagnóstico finalizado.\n")
}
