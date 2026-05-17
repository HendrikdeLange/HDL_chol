set.seed(123)
df <- readRDS("data/nhanes_subset.rds")
colnames(df)

#Loading Data
obs_num <- length(df$DirectChol)
#Need to create x_mat 
df$AlcoholCategory <- relevel(factor(df$AlcoholCategory), ref = "Abstainer")
x_mat <- model.matrix(~ ., data = df[, colnames(df) != "DirectChol"])
k <- ncol(x_mat)
y_vec <- df$DirectChol


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
colnames(beta_draws) <- colnames(x_mat)

# Posterior summaries
post_means <- colMeans(beta_draws)
post_ci    <- apply(beta_draws, 2, quantile, probs = c(0.025, 0.975))

round(rbind(mean = post_means, post_ci), 4)
mean(sigma2_draws)

# Trace plots to verify sampler behaviour 
par(mfrow = c(2, 2))
for (j in 1:min(4, k)) {
  plot(beta_draws[, j], type = "l",
       main = paste("Trace:", colnames(x_mat)[j]),
       ylab = colnames(x_mat)[j], xlab = "Iteration")
}

source("final/part2/subjective_bayesian_regression.R")
par(mfrow = c(2,1))
plot(sigma2_draws, type = "l",
     ylab = "Variance", main = "Objective Bayesian Regression Residual Variance")

plot(sigma2_posterior_vec, type = "l",
     ylab = "Variance", main = "Subjective Bayesian Regression Residual Variance")
