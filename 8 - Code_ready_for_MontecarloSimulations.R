K <- function(x, z, sigma2, ell) {
  sigma2 * exp(-((x - z)^2) / (2 * ell^2))
}

krr_fit <- function(t_train, y_train, ell, lambda, sigma2){
  D <- outer(t_train, t_train, "-")^2
  K_mat <- sigma2 * exp( - D / (2 * ell^2) )
  K_reg <- K_mat + lambda * diag(length(t_train))
  alpha <- solve(K_reg, y_train)
  list(
    alpha   = alpha,
    t_train = t_train,
    ell     = ell,
    sigma2  = sigma2
  )
}


krr_predict <- function(model, t_new){
  D <- outer(t_new, model$t_train, "-")^2
  K_new <- model$sigma2 * exp( - D / (2 * model$ell^2) )
  as.numeric(K_new %*% model$alpha)
}

ell_grid    <- c(5,10,15,20,25,30,35, 40,45,50,55,60,65,70,75, 80)
lambda_grid <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)
sigma2 <- 1
