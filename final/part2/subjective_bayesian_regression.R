library(mvtnorm)
set.seed(123)

#Loading Data
obs_num <- length(df$DirectChol)
#Need to create x_mat 
df$AlcoholCategory <- relevel(factor(df$AlcoholCategory), ref = "Abstainer")
x_mat <- model.matrix(~ ., data = df[, colnames(df) != "DirectChol"])
colnames(x_mat)
k <- ncol(x_mat)
y_vec <- df$DirectChol
mu_prior_vec <- c(
  0.65,   # Intercept — aligns with frequentist estimate at covariate zero
  0.08,   # AlcoholCategoryHeavy
  0.07,   # AlcoholCategoryModerate
  0.002,  # Age — POSITIVE (HDL rises slightly with age)
  -0.22,   # Gendermale
  -0.018,  # BMI
  0.001,  # BPSysAve — POSITIVE (per frequentist estimate)
  -0.07,   # DiabetesYes — attenuated from -0.15
  -0.09,   # SmokeNowYes
  0.05    # PhysActiveYes — attenuated from 0.10
)

cov_prior_mat <- diag(c(
  0.09,     # Intercept       SD = 0.30
  0.04,     # Heavy           unchanged
  0.02,     # Moderate        unchanged
  0.000025, # Age             unchanged — effect per year is tiny
  0.016,    # Gendermale      unchanged
  0.0004,   # BMI             unchanged
  0.000025, # BPSysAve        unchanged
  0.02,     # DiabetesYes     tightened SD = 0.14 (was 0.20)
  0.02,     # SmokeNowYes     unchanged
  0.015     # PhysActiveYes   tightened SD = 0.12 (was 0.14)
), k)

a_prior <- 10       # unchanged
b_prior <- 1.1871   # unchanged
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



