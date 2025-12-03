############################################################
# Capítulo 3 – Detrending con KRR (RBF) + ADF / KPSS
# Script completo
############################################################

rm(list = ls())
set.seed(123)

#-----------------------------------------------------------
# BLOQUE A: Librerías y funciones base
#-----------------------------------------------------------

library(tseries)
library(ggplot2)

## Tendencia verdadera g(t)  (puedes cambiarla aquí)
g_fun <- function(t) {
  g <- 0.01 * t + 2 * sin(0.01 * t)
  # Otras opciones que puedes probar luego:
  # 0.0000001 * t^3 - 0.02 * t^2 + 5 * sin(0.03 * t)
  # 0.1 * t + 5 * sin(0.2 * t) + 5 * cos(0.15 * t)
  # ifelse(t < length(t)/2, 0.1 * t, 0.1 * t + 50)
  # 0.5 * sqrt(t) + 7 * sin(0.01 * t) * cos(0.02 * t)
}

## Generar serie Y_t = g(t) + u_t, con u_t ~ AR(1)
gen_series <- function(n, phi, sigma2 = 1) {
  t <- 1:n
  g <- g_fun(t)
  
  u <- numeric(n)
  u[1] <- rnorm(1)
  for (i in 2:n) {
    u[i] <- phi * u[i - 1] + rnorm(1, sd = sqrt(sigma2))
  }
  
  y <- g + u
  list(t = t, y = y, g = g, u = u)
}

## Kernel Ridge Regression (KRR) con kernel RBF
krr_fit <- function(t_train, y_train, ell, lambda, sigma2 = 1) {
  D     <- outer(t_train, t_train, "-")^2
  K_mat <- sigma2 * exp(- D / (2 * ell^2))
  K_reg <- K_mat + lambda * diag(length(t_train))
  alpha <- solve(K_reg, y_train)
  list(alpha = alpha, t_train = t_train, ell = ell, sigma2 = sigma2)
}

krr_predict <- function(model, t_new) {
  D     <- outer(t_new, model$t_train, "-")^2
  K_new <- model$sigma2 * exp(- D / (2 * model$ell^2))
  as.numeric(K_new %*% model$alpha)
}

#-----------------------------------------------------------
# BLOQUE B: Cross-validation de KRR por tamaño muestral n
#-----------------------------------------------------------

## Grids de hiperparámetros
ell_grid    <- c(5, 10, 15, 20, 25, 30, 35, 40,
                 45, 50, 55, 60, 65, 70, 75, 80,
                 100, 120, 140, 150)
lambda_grid <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

## CV para un n dado (phi_cv = 0.6 como ruido moderado)
cv_krr_for_n <- function(n, phi_cv = 0.6,
                         sigma2 = 1,
                         train_frac = 0.6,
                         valid_frac = 0.3) {
  
  dat <- gen_series(n, phi_cv, sigma2)
  t <- dat$t
  y <- dat$y
  
  n_train <- floor(train_frac * n)
  n_valid <- floor(valid_frac * n)
  
  t_train <- t[1:n_train]
  y_train <- y[1:n_train]
  
  t_valid <- t[(n_train + 1):(n_train + n_valid)]
  y_valid <- y[(n_train + 1):(n_train + n_valid)]
  
  # usa los vectores globales ell_grid y lambda_grid
  results <- expand.grid(ell = ell_grid, lambda = lambda_grid)
  results$mse <- NA_real_
  
  for (i in seq_len(nrow(results))) {
    ell_i    <- results$ell[i]
    lambda_i <- results$lambda[i]
    
    model_i <- krr_fit(t_train, y_train, ell_i, lambda_i, sigma2)
    y_pred_valid <- krr_predict(model_i, t_valid)
    
    results$mse[i] <- mean((y_valid - y_pred_valid)^2)
  }
  
  ell_star    <- results$ell[which.min(results$mse)]
  lambda_star <- results$lambda[which.min(results$mse)]
  
  list(ell_star = ell_star,
       lambda_star = lambda_star,
       cv_table = results)
}

## Grillas globales de n y phi (mismas que en el capítulo 2)
n_grid   <- c(30, 50, 100, 200, 400, 700)
phi_grid <- c(0, 0.5, 0.8, 0.9, 0.95, 0.99)

## Hiperparámetros óptimos por cada n (se calcula una sola vez)
hyper_list <- lapply(n_grid, function(nn) cv_krr_for_n(nn))

hyper_df <- data.frame(
  n           = n_grid,
  ell_star    = sapply(hyper_list, `[[`, "ell_star"),
  lambda_star = sapply(hyper_list, `[[`, "lambda_star")
)

cat("Hiperparámetros óptimos por n:\n")
print(hyper_df)

#-----------------------------------------------------------
# BLOQUE C: Simulación Monte Carlo con KRR + ADF / KPSS
#-----------------------------------------------------------

## Una simulación Monte Carlo para un (n, phi) dado
## Devuelve:
##  - power_adf: proporción de rechazos de ADF (potencia)
##  - size_kpss: proporción de rechazos de KPSS (tamaño empírico)
sim_krr_tests <- function(n, phi, B,
                          ell_star, lambda_star,
                          sigma2 = 1,
                          alpha = 0.05) {
  
  rej_adf  <- 0
  rej_kpss <- 0
  
  for (b in 1:B) {
    dat <- gen_series(n, phi, sigma2)
    t <- dat$t
    y <- dat$y
    
    # Ajuste KRR con hiperparámetros fijados
    model <- krr_fit(t, y, ell_star, lambda_star, sigma2)
    g_hat <- krr_predict(model, t)
    r     <- y - g_hat  # residuos
    
    # Test ADF (H0: raíz unitaria)
    adf_res <- adf.test(r)
    if (adf_res$p.value < alpha) rej_adf <- rej_adf + 1
    
    # Test KPSS (H0: estacionariedad en nivel)
    kpss_res <- kpss.test(r, null = "Level")
    if (kpss_res$p.value < alpha) rej_kpss <- rej_kpss + 1
  }
  
  c(
    power_adf = rej_adf / B,
    size_kpss = rej_kpss / B
  )
}

#-----------------------------------------------------------
# BLOQUE D: Loop sobre toda la grilla (n, phi)
#-----------------------------------------------------------

## Número de simulaciones por combinación (n, phi)
B <- 100  # puedes bajar a 500 para probar más rápido

res_list <- list()
idx <- 1

for (nn in n_grid) {
  ell_star    <- hyper_df$ell_star[hyper_df$n == nn]
  lambda_star <- hyper_df$lambda_star[hyper_df$n == nn]
  
  for (ph in phi_grid) {
    sim_res <- sim_krr_tests(
      n           = nn,
      phi         = ph,
      B           = B,
      ell_star    = ell_star,
      lambda_star = lambda_star
    )
    
    res_list[[idx]] <- data.frame(
      n         = nn,
      phi       = ph,
      power_adf = sim_res["power_adf"],
      size_kpss = sim_res["size_kpss"]
    )
    
    cat("n =", nn,
        "phi =", ph,
        "→ ADF:",  round(sim_res["power_adf"], 3),
        "| KPSS:", round(sim_res["size_kpss"], 3), "\n")
    
    idx <- idx + 1
  }
}

res_krr <- do.call(rbind, res_list)

cat("\nResumen de resultados (post-KRR):\n")
print(res_krr)

#-----------------------------------------------------------
# BLOQUE E: Gráficas de resultados (post-KRR)
#-----------------------------------------------------------

## Curvas de potencia del ADF sobre residuos KRR
p_adf <- ggplot(res_krr,
                aes(x = n, y = power_adf,
                    color = factor(phi),
                    group = factor(phi))) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = n_grid) +
  labs(
    title = "Curvas de potencia del test ADF\nsobre residuos KRR (Y_t = g(t) + u_t, u_t AR(1) estacionario)",
    x = "Tamaño de muestra (n)",
    y = "Potencia (proporción de rechazos)",
    color = expression(phi)
  ) +
  theme_minimal()

## Proporción de rechazos de KPSS (tamaño empírico) sobre residuos KRR
p_kpss <- ggplot(res_krr,
                 aes(x = n, y = size_kpss,
                     color = factor(phi),
                     group = factor(phi))) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(breaks = n_grid) +
  labs(
    title = "Proporción de rechazos del test KPSS\nsobre residuos KRR (Y_t = g(t) + u_t, u_t AR(1) estacionario)",
    x = "Tamaño de muestra (n)",
    y = "Proporción de rechazos",
    color = expression(phi)
  ) +
  theme_minimal()

## Mostrar en pantalla
print(p_adf)
print(p_kpss)

## Si quieres guardarlas como PNG:
# ggsave("curvas_potencia_adf_post_krr.png", p_adf,  width = 7, height = 5, dpi = 300)
# ggsave("rechazos_kpss_post_krr.png",        p_kpss, width = 7, height = 5, dpi = 300)
