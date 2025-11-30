library(tseries)

sim_adf_kpss_grid <- function(phi_vals,
                              n_vals,
                              B      = 2000,
                              alpha  = 0.05,
                              k_lags = 1,
                              kpss_null = c("Level", "Trend"),
                              seed   = 123) {
  kpss_null <- match.arg(kpss_null)
  if (!is.null(seed)) set.seed(seed)
  
  grid <- expand.grid(phi = phi_vals, n = n_vals)
  grid$rej_adf  <- NA_real_
  grid$rej_kpss <- NA_real_
  
  for (k in seq_len(nrow(grid))) {
    phi_k <- grid$phi[k]
    n_k   <- grid$n[k]
    
    rej_adf  <- 0L
    rej_kpss <- 0L
    
    for (b in seq_len(B)) {
      # AR(1) estacionario: u_t = phi * u_{t-1} + eps_t
      u <- arima.sim(model = list(ar = phi_k), n = n_k)
      
      # ADF: H0 = raíz unitaria
      adf_res <- adf.test(u, k = k_lags)
      if (adf_res$p.value < alpha) rej_adf  <- rej_adf  + 1L
      
      # KPSS: H0 = estacionariedad (en nivel o tendencia)
      kpss_res <- kpss.test(u, null = kpss_null)
      if (kpss_res$p.value < alpha) rej_kpss <- rej_kpss + 1L
    }
    
    grid$rej_adf[k]  <- rej_adf  / B
    grid$rej_kpss[k] <- rej_kpss / B
    
    cat("phi =", phi_k, "n =", n_k,
        "→ ADF:",  round(grid$rej_adf[k],  3),
        "| KPSS:", round(grid$rej_kpss[k], 3), "\n")
  }
  
  grid
}

phi_vals <- c(0, 0.5, 0.8, 0.9, 0.95, 0.99)
n_vals   <- c(30, 50, 100, 200, 400, 700)

res_grid <- sim_adf_kpss_grid(phi_vals, n_vals,
                              B = 2000,
                              alpha = 0.05,
                              k_lags = 1,
                              kpss_null = "Level")  # aquí H0 = estacionario en nivel


library(ggplot2)

ggplot(res_grid,
       aes(x = n, y = rej_adf,
           colour = factor(phi), group = factor(phi))) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Tamaño de muestra (n)",
       y = "Proporción de rechazos",
       colour = expression(phi),
       title = "Rechazos del test KPSS bajo AR(1) estacionario") +
  theme_minimal()

