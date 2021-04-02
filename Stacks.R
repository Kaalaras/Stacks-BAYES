library(MASS)
stacks_func <- function(nchain, data, priors, sigma2_init){
  #Pour les priors : d'abord ceux de sigma2_beta puis sigma2_y
  #les 4 beta j puis sigma2_y puis sigma2_beta
  chain <- matrix(NA, nrow = nchain + 1, ncol = 6)
  
  nobs <- nrow(data)
  mu <- rep(NA, nobs)
  
  #on initialise la chaine
  beta <- rep(NA, 4)
  sigma2_y <- sigma2_init[1]
  sigma2_beta <- sigma2_init[2]
  
  update_var <- sigma2_beta*sigma2_y/(sigma2_y + nobs*sigma2_beta)
  update_mean <- sum(data[,4])/sigma2_y*update_var
    for(j in 1:4){  
    beta[j] <- update_mean
  }
  init <- c(beta, sigma2_init)
  chain[1,] <- init
  
  #pour calculer les zij
  sd_var <- rep(NA, 3)
  means <- rep(NA, 3)
  for (j in 1:3){
    sd_var[j] <- sd(data[,j])
    means[j] <- mean(data[,j])
  }
  #calcul des zij
  z <- matrix(NA, nobs, 3)
  for (i in 1:nobs){
    for (j in 1:3){
      z[i,j] <- (data[i,j] - means[j])/sd_var[j]
    }
  }
  z1 <- cbind(matrix(1,nobs,1),z)
  #on initialise aussi les outliers
  mu_init <-  beta[1] + beta[2]*z[,1] + beta[3]*z[,2] + beta[4]*z[,3]
  outliers <- matrix(NA, nchain+1, nobs)
  values <- (data[,4] - mu_init)/chain[1,4]
  outliers[1,] <- (values > 2.5)*1 + (values < -2.5 )*1
  
  #on ajoute en sortie une chaine pour les MSE
  MSE <- matrix(0, nchain+1, 1)
  MSE[1] <- 1/nobs*sum((data[,4]- mu_init)^2)
  for (iter in 1:nchain){
    #On récupère les éléments de la chaîne
    sigma2_y <- chain[iter, 5]
    sigma2_beta <- chain[iter, 6]
    beta <- chain[iter, 1:4]
    #MAJ des beta
    
    #Methode 1 : multivarié (ça ne marche pas)
    # sum <- 0
    # for (i in 1:nobs){
    #   sum <- sum + data[i,4]*z1[i,]
    # }
    # mu_star <- 1/sigma2_y * sum
    # sum_zi <- matrix(0,4,4)
    # for (i in 1:nobs){
    #   vect <- matrix(z1[i,], 4, 1)
    #   sum_zi <- sum_zi + vect%*%t(vect)
    # }
    # sigma_star_inv <- 1/sigma2_beta*diag(4) + 1/sigma2_y*sum_zi
    # 
    # beta <- mvrnorm(1, mu_star, solve(sigma_star_inv))
    
    
    
    for(j in 1:4){
      update_var <- 1/(sum(z1[,j]^2)/sigma2_y + 1/sigma2_beta)
      sum_temp <- 0
      
      update_mean <- sum(z1[,j]*(data[,4] - beta[-j]*z1[,-j]))* update_var/sigma2_y

      beta[j] <- rnorm(1 , update_mean, sqrt(update_var))
    }
    
    #MAJ de sigma2_y
    ##On calcule les mui
    for (i in 1:nobs){
      mu[i] <- beta[1] + beta[2]*z[i,1] + beta[3]*z[i,2] + beta[4]*z[i,3]
    }
    update_shape <- priors[3] + nobs*0.5
    update_rate <- priors[4] + 0.5*sum((data[,4]-mu)^2)

    sigma2_y <- 1/ rgamma(1, update_shape, update_rate)
    
    #MAJ de sigma2_beta
    update_shape <- priors[1] + 4*0.5
    update_rate <- priors[2] + 0.5*sum(beta^2)
    
    sigma2_beta <- 1/ rgamma(1, update_shape, update_rate)
    
    ###Finalement
    #MAJ de la chaine!!
    chain[iter+1,] <- c(beta, sigma2_y, sigma2_beta)
    #Et des outliers !
    values <- (data[,4] - mu)/sigma2_y
    outliers[iter+1,] <- (values > 2.5)*1 + (values < -2.5 )*1
    
    #MSE
    MSE[iter+1] <- 1/nobs*sum((data[,4]-mu)^2)

  }
  #les b normalisés
  b <- matrix(NA, nchain+1, 4)
  for (j in 2:4){
    b[,j] <- chain[,j]/sd_var[j-1]}
  b[,1] <- chain[,1] - b[,2]*means[1] - b[,3]*means[2] - b[,4]*means[3]
  
  
  return(list(chain=chain, outliers=outliers, b=b, MSE=MSE))
}

### Application
p <-3
N <-21
Y <-c(43, 37, 37, 28, 18, 18, 19, 20, 15, 14, 14, 13, 11, 12, 8, 
    7, 8, 8, 9, 15, 15)
x <-structure(c(80, 80, 75, 62, 62, 62, 62, 62, 59, 58, 58, 58, 58, 
              58, 50, 50, 50, 50, 50, 56, 70, 27, 27, 25, 24, 22, 23, 24, 24, 
              23, 18, 18, 17, 18, 19, 18, 18, 19, 19, 20, 20, 20, 89, 88, 90, 
              87, 87, 87, 93, 93, 87, 80, 89, 88, 82, 93, 89, 86, 72, 79, 80, 
              82, 91), .Dim = c(21, 3))

data <- cbind(x, Y)
out <- stacks_func(10^4, data, priors = c(100, 100, 20, 320), sigma2_init = c(4,4))
chain <- out$chain[-(1:500),]

par(mfrow = c(2, 3), mar = c(4, 5, 0.5, 0.5))
ylabs <- c(expression(beta[0]), expression(beta[1]),expression(beta[2]),expression(beta[3]), 'sigma2_y', 'sigma2_beta')
for (j in 1:6)
  plot(chain[,j], type = "l", ylab = ylabs[j])

blabs <- c('b0', 'b1', 'b2', 'b3')
par(mfrow = c(2, 2), mar = c(4, 5, 0.5, 0.5))
for (j in 1:4)
  plot(out$b[,j], type = "l", ylab = blabs[j])
