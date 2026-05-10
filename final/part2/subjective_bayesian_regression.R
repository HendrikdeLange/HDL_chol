library(mvtnorm)
set.seed(123)

#Loading Data
obs_num <- length(df$DirectChol)
#Need to create x_mat 
df$AlcoholCategory <- relevel(factor(df$AlcoholCategory), ref = "Abstainer")
x_mat <- model.matrix(~ ., data = df[, colnames(df) != "DirectChol"])
k <- ncol(x_mat)
y_vec <- df$DirectChol

#Step 1 (Prior Hyperparameters)
#We need to decide on these, just keeping uninformative as is
mu_prior_vec <- rep(0, k) #mu = 0 since we have no idea
cov_prior_mat <- diag(1, k) #How strong we hold our higher belief, larger diags = more vague
#!!!Calibration Trick
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
colnames(beta_posterior_mat) <- colnames(x_mat)

# Posterior summaries
post_means <- colMeans(beta_posterior_mat)
post_ci    <- apply(beta_posterior_mat, 2, quantile, probs = c(0.025, 0.975))

round(rbind(mean = post_means, post_ci), 4)
mean(sigma2_posterior_vec)

# Trace plots to verify sampler behaviour (Just say its riht in the presentation)
par(mfrow = c(2, 2))
for (j in 1:min(4, k)) {
  plot(beta_posterior_mat[, j], type = "l",
       main = paste("Trace:", colnames(x_mat)[j]),
       ylab = colnames(x_mat)[j], xlab = "Iteration")
}



