set.seed(123)

#Going to Try Baysian Regression on my own using the algorithm
#Objective One First

#Simulating the data
obs_num <- 100
x_mat <- cbind(1, rnorm(obs_num), rnorm(obs_num))  # Intercept + two predictors
k <- ncol(x_mat)
beta_true <- c(3, 2, 4)
sigma_true <- sqrt(1.5)
y_vec <- x_mat %*% beta_true + rnorm(obs_num,sd =  sigma_true) #Observation Vector


# Step 1
beta_hat_vec <- solve(t(x_mat) %*% x_mat) %*% t(x_mat) %*% y_vec

#step 2
s2 <- as.numeric((1/(obs_num - k )) * ((t(y_vec - x_mat %*% beta_hat_vec)) %*% (y_vec - x_mat %*% beta_hat_vec)))

#Step 3
l_mat <- t(chol(solve(t(x_mat) %*% x_mat)))

#Step 4 ...
N            <- 10000
beta_draws   <- matrix(NA, nrow = N, ncol = k)
sigma2_draws <- numeric(N)

for (i in seq_len(N)) {
  u               <- rchisq(1, df = obs_num - k)
  sigma2_i        <- (obs_num - k) * s2 / u
  z_vec           <- rnorm(k)
  beta_i          <- beta_hat_vec + sqrt(sigma2_i) * l_mat %*% z_vec
  beta_draws[i, ] <- beta_i
  sigma2_draws[i] <- sigma2_i
}



#Plots + Analyses
apply(beta_draws, 2, mean)
mean(sigma2_draws)

par(mfrow = c(2,2))
plot(beta_draws[,1], type = "l", main = "Trace: beta_0")
plot(beta_draws[,2], type = "l", main = "Trace: beta_1")
plot(sigma2_draws, type = "l", main = "Trace: sigma^2")
hist(sigma2_draws, breaks = 30, main = "Posterior of sigma^2")

#Confidence Intervals
ci_bayes <- apply(beta_draws, 2, quantile, probs = c(0.025, 0.975))
rownames(ci_bayes) <- c("2.5%", "97.5%")
colnames(ci_bayes) <- c("beta_0", "beta_1", "beta_2")
ci_bayes
