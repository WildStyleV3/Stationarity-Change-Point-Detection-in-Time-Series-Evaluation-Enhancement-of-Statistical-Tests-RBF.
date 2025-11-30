set.seed(123)
n <- 700
t <- 1:n
#g <- 0.01 * t + 2 * sin(0.01 * t)
#g <- 0.0000001 * t^3 - 0.02 * t^2 + 5*sin(0.03*t)
#g <- 0.1*t + 5*sin(0.2 * t) + 5*cos(0.15 * t)
#g <- ifelse(t < 350, 0.1*t, 0.1*t + 50)
#g <- 0.5*sqrt(t) + 7*sin(0.01*t)*cos(0.02*t)

phi <- 0.6
u <- numeric(n)
u[1] <- rnorm(1)
for (i in 2:n) {
  u[i] <- phi * u[i-1] + rnorm(1, sd = 1)
}
y <- g + u #Funcion 
plot(y)

n_train <- 400
n_valid <- 200

t_train <- t[1:n_train] #Tiempos del 1 al 400
y_train <- y[1:n_train] #Valores de la serie en esos tiempos

t_valid <- t[(n_train+1):(n_train+n_valid)] #con esto validamos el modelo
y_valid <- y[(n_train+1):(n_train+n_valid)] #con esto validamos el modelo 

sigma2 = 1

#KERNEL RBF 
K <- function(x, z, sigma2, ell) {
  sigma2 * exp(-((x - z)^2) / (2 * ell^2))
}



ell_grid    <- c(5,10,15,20,25,30,35, 40,45,50,55,60,65,70,75, 80)
lambda_grid <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)


#Álgebra Matricial 
######################################################################
krr_fit <- function(t_train, y_train, ell, lambda, sigma2){
  D     <- outer(t_train, t_train, "-")^2
  K_mat <- sigma2 * exp(- D / (2 * ell^2))
  K_reg <- K_mat + lambda * diag(length(t_train))
  alpha <- solve(K_reg, y_train)
  list(alpha = alpha, t_train = t_train, ell = ell, sigma2 = sigma2)
}
######################################################################
  
krr_predict <- function(model, t_new){
  D     <- outer(t_new, model$t_train, "-")^2
  K_new <- model$sigma2 * exp(- D / (2 * model$ell^2))
  as.numeric(K_new %*% model$alpha)
}
results <- expand.grid(ell = ell_grid, lambda = lambda_grid)
results$mse <- NA_real_

#CROSS VALIDATION
####################################################################
for (i in seq_len(nrow(results))) {                                #
  ell_i    <- results$ell[i]                                       #
  lambda_i <- results$lambda[i]                                    #
                                                                   #
  model_i <- krr_fit(t_train, y_train, ell_i, lambda_i, sigma2)    #
  y_pred_valid <- krr_predict(model_i, t_valid)#DATOS DE VALIDACION#
                                                                   #
  results$mse[i] <- mean((y_valid - y_pred_valid)^2) #MSE CALCULAR #
}                                                                  #
results[order(results$mse), ][1:5, ] #Ver las mejores combinaciones#
####################################################################
print(results)


#HIPERPARAMETROS MÁS ÓPTIMOS
##############################################################
ell_star    <- results$ell[which.min(results$mse)]
lambda_star <- results$lambda[which.min(results$mse)]
##############################################################

model_final <- krr_fit(t, y, ell_star, lambda_star, sigma2)
g_hat <- krr_predict(model_final, t)  # K(RAW DATA)
r <- y - g_hat    # RAW DATA - K

par(mfrow = c(3,1), mar = c(3,4,3,1))

# 1) RAW DATA
plot(t, y, type = "l", col = "gray40",
     main = "RAW DATA", xlab = "t", ylab = "y")
lines(t, g, col = "blue", lwd = 2)   # tendencia verdadera

# 2) KRR (tendencia estimada)
plot(t, g_hat, type = "l", col = "red", lwd = 2,
     main = "KRR (tendencia estimada)", xlab = "t", ylab = "g_hat")
lines(t, g, col = "blue", lwd = 1, lty = 2)  # para comparar

# 3) Residuos
plot(t, r, type = "l", col = "black",
     main = "Residuos: y - g_hat", xlab = "t", ylab = "r")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1,1))


library(tseries)

#Testing Stationarity

adf_result <- adf.test(r)
adf_result


kpss_result <- kpss.test(r,null='Level')
kpss_result





