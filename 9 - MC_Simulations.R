rm(list = ls())        
source("8 - Code_ready_for_MontecarloSimulations.R")

library(tseries)
N_sim <- 100
count_adf <- 0
count_kpss <- 0 
count_both <- 0
set.seed(123)
for (sim in 1:N_sim) {
phi = 1.2
  n <- 700
  t <- 1:n
  g <- 0.1*t + 5*sin(0.2 * t) + 5*cos(0.15 * t)

  u<-numeric(n)
  u[1] <- rnorm(1)
  for (i in 2:n){
    u[i] <- phi*u[i-1] + rnorm(1)
  }
  y<-g+u
  n_train <- 400
  n_valid <- 200
  t_train <- t[1:n_train]
  y_train <- y[1:n_train]
  
  t_valid <- t[(n_train+1):(n_train+n_valid)]
  y_valid <- y[(n_train+1):(n_train+n_valid)]
  results <- expand.grid(ell = ell_grid,lambda=lambda_grid)
  results$mse <- NA_real_
  for(i in seq_len(nrow(results))) {
    ell_i <- results$ell[i]
    lambda_i <- results$lambda[i]
    model_i <- krr_fit(t_train,y_train,ell_i,lambda_i,sigma2)
    y_pred_valid <- krr_predict(model_i,t_valid)
    results$mse[i] <- mean((y_valid - y_pred_valid)^2)
  }
  ell_star <- results$ell[which.min(results$mse)]
  lambda_star <- results$lambda[which.min(results$mse)]
  
  model_final <- krr_fit(t,y,ell_star,lambda_star,sigma2 =1)
  g_hat <- krr_predict(model_final,t)
  r <- y - g_hat
  
  adf_res <- adf.test(r)
  adf_stationary <-(adf_res$p.value <0.05)
  kpss_res <- kpss.test(r,null = 'Level')
  kpss_stationary <- (kpss_res$p.value >0.05)
  
  if (adf_stationary) count_adf <- count_adf +1
  if (kpss_stationary) count_kpss <- count_kpss+1
  if(adf_stationary && kpss_stationary) count_both <- count_both+1
  
  cat('Sim',sim,'ADF',adf_stationary,
      'KPSS',kpss_stationary, '\n')

}

cat("\nRESULTADOS SOBRE", N_sim, "SIMULACIONES:\n\n")
cat("ADF detectó estacionariedad en:",  count_adf,  "simulaciones\n")
cat("KPSS detectó estacionariedad en:", count_kpss, "simulaciones\n")
cat("Ambos detectaron estacionariedad:", count_both, "simulaciones\n")

cat("\nProporción ADF  =", count_adf  / N_sim, "\n")
cat("Proporción KPSS =", count_kpss / N_sim, "\n")
cat("Proporción ambos=", count_both / N_sim, "\n")
