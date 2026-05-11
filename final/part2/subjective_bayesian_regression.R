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

#Step 1 (Prior Hyperparameters)
mu_prior_vec <- c(
  # Intercept — baseline HDL for reference group (Female, Non-smoker,
  # Non-diabetic, Inactive Abstainer) literature suggests ~1.3 mmol/L
  1.30,   # (Intercept)
  
  # Treatment effects on HDL (vs Abstainer reference)
  # Heavy drinking: mixed evidence, net slightly positive but attenuated
  0.05,  # AlcoholCategoryHeavy
  # Moderate drinking: most consistent positive HDL effect in literature
  0.10,  # AlcoholCategoryModerate
  
  # Covariates
  # Age: HDL tends to rise slightly with age in women, small overall effect
  -0.002, # Age (per year, small)
  # Gender (male): men have systematically lower HDL than women (~0.25 mmol/L)
  -0.25,  # Gendermale
  # BMI: well-established inverse relationship, ~0.02 mmol/L per unit BMI
  -0.02,  # BMI (per kg/m²)
  # Systolic BP: weak negative association with HDL
  -0.002, # BPSysAve (per mmHg)
  # Diabetes: diabetics have markedly lower HDL
  -0.15,  # DiabetesYes
  # Smoking: smokers have lower HDL, well-replicated
  -0.10,  # SmokeNowYes
  # Physical activity: active individuals have higher HDL
  0.10   # PhysActiveYes
)

cov_prior_mat <- diag(c(
  # Intercept — fairly confident baseline HDL is near 1.3, allow ± 0.5 range
  0.09,   # (Intercept)       SD = 0.30 mmol/L
  
  # Treatment effects — moderate uncertainty, direction known from literature
  0.04,   # AlcoholCategoryHeavy      SD = 0.20; evidence mixed, wider
  0.02,   # AlcoholCategoryModerate   SD = 0.14; more consistent literature
  
  # Covariates
  0.000025, # Age          SD = 0.005; effect per year is tiny, tight
  0.016,    # Gendermale   SD = 0.13;  well-established, fairly tight
  0.0004,   # BMI          SD = 0.02;  per unit effect small, consistent
  0.000025, # BPSysAve     SD = 0.005; very small per-mmHg effect
  0.04,     # DiabetesYes  SD = 0.20;  meaningful but variable across studies
  0.02,     # SmokeNowYes  SD = 0.14;  reasonably consistent literature
  0.02      # PhysActiveYes SD = 0.14; moderate evidence, some variability
), k)



a_prior <- 10
b_prior <- 1.1871


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



