library(mvtnorm)
set.seed(123)

#Going to Try Baysian Regression on my own using the algorithm
#Subjective Second

#Simulating the data
obs_num <- 100
x_mat <- cbind(1, rnorm(obs_num), rnorm(obs_num))  # Intercept + two predictors
k <- ncol(x_mat)
beta_true <- c(3, 2, 4)
sigma_true <- sqrt(1.5)
y_vec <- x_mat %*% beta_true + rnorm(obs_num,sd =  sigma_true) #Observation Vector

#Step 1 (Prior Hyperparameters)
mu_prior_vec <- rep(0, k) #mu = 0 since we have no idea
cov_prior_mat <- diag(1000, k) #How strong we hold our higher belief, larger diags = more vague
a_prior <- 0.01 #very uninformative
b_prior <- 0.01 #very uninformative

  
#Step 2 (Initial Value)
beta_prior <- solve(t(x_mat) %*% x_mat) %*% t(x_mat) %*% y_vec #OLS estimate
sigma2 <- sigma2 <- as.numeric(t(y_vec - x_mat %*% beta_prior) %*% (y_vec - x_mat %*% beta_prior)) / (obs_num - k) #OLS estimate


#Step 3
N <- 10000
B <- 1000
beta_posterior_mat <- matrix(NA, nrow = N, ncol = k)
sigma2_posterior_vec <- numeric(N)
for (i in seq_len(N)) {
  cov_mat_posterior <- solve(solve(cov_prior_mat) + (1 / sigma2) * t(x_mat) %*% x_mat)
  mu_posterior <- cov_mat_posterior %*% (solve(cov_prior_mat) %*% mu_prior_vec + (1 / sigma2) * t(x_mat) %*% y_vec)
  #Now need to draw Beta (Just use mvt norm for simplicity)
  beta_i <-  as.vector(rmvnorm(1, mean = mu_posterior, sigma = cov_mat_posterior))
  residuals <- y_vec - (x_mat %*% beta_i)
  a_n <- a_prior + obs_num / 2
  b_n <- b_prior + sum(residuals^2) / 2
  sigma2_i <- 1/ rgamma(1, a_n, b_n)
  sigma2 <- sigma2_i
  beta_posterior_mat[i,] <- beta_i
  sigma2_posterior_vec[i] <- sigma2_i 
}

#Remove Burn in
beta_posterior_mat <- beta_posterior_mat[-(1:B),]
sigma2_posterior_vec <- sigma2_posterior_vec[-(1:B)]


#Analyses + Plots
apply(beta_posterior_mat, 2, mean)
mean(sigma2_posterior_vec)
par(mfrow = c(2,2))
plot(beta_posterior_mat[,1], type = "l", main = "Trace: beta_0")
plot(beta_posterior_mat[,2], type = "l", main = "Trace: beta_1")
plot(sigma2_posterior_vec, type = "l", main = "Trace: sigma^2")
hist(sigma2_posterior_vec, breaks = 30, main = "Posterior of sigma^2")



par(mfrow = c(2,1))
