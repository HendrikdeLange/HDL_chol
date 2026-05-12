library(mice)
library(dplyr)
library(ggplot2)
set.seed(123)

# --- Load data ---
df <- readRDS("data/nhanes_subset.rds")

# --- Step 1: Treatment setup ---
df$AlcoholCategory <- relevel(factor(df$AlcoholCategory), ref = "Abstainer")
trt_levels         <- levels(df$AlcoholCategory)   # "Abstainer", "Moderate", "Heavy"

# --- Step 2: Create dummy potential outcome columns ---
for (treat in trt_levels) {
  df[[treat]] <- ifelse(df$AlcoholCategory == treat, df$DirectChol, NA)
}

# --- Step 3: Build imputation dataset ---
# Include covariates + all potential outcome columns
# Exclude original response and treatment label
imp_vars <- c(colnames(df)[!(colnames(df) %in% c("DirectChol", "AlcoholCategory",
                                                 trt_levels))],
              trt_levels)
imp_data <- df[, imp_vars]

# --- Step 4: Run MICE ---
mice_data <- mice(imp_data, m = 50, maxit = 10, method = "pmm", seed = 123) 




#Bullets
#Bullet 2

collapsed <- complete(mice_data, action = "long", include = FALSE) %>%
  group_by(.id) %>%                        
  summarise(
    Abstainer = mean(Abstainer, na.rm = TRUE),
    Heavy     = mean(Heavy,     na.rm = TRUE),
    Moderate  = mean(Moderate,  na.rm = TRUE),
    .groups = "drop"
  )

#This is now one datasets with the mean of impuations
#calculate the causal effects
collapsed$moderate_abstainer <- collapsed$Moderate - collapsed$Abstainer
collapsed$heavy_abstainer <- collapsed$Heavy - collapsed$Abstainer

#Histogram of Moderate - Abstainer
ggplot(collapsed, aes(x = moderate_abstainer)) +
  geom_histogram()

#Histogram of Heavy - Abstainer
ggplot(collapsed, aes(x = heavy_abstainer)) +
  geom_histogram()

#Percentiles
#MICE Moderate - Abstainer
quantile(sort(collapsed$moderate_abstainer), probs=c(0.025,0.5, 0.975))
# MICE Heavy - Abstainer
quantile(sort(collapsed$heavy_abstainer), probs=c(0.025,0.5, 0.975))


#Bullet 3
#Use entire data set again
df_vectors <- complete(mice_data, action = "long", include = FALSE)
df_vectors$moderate_abstainer <- df_vectors$Moderate - df_vectors$Abstainer
df_vectors$heavy_abstainer <- df_vectors$Heavy - df_vectors$Abstainer


#Histogram of Moderate - Abstainer 
ggplot(df_vectors, aes(x = moderate_abstainer)) +
  geom_histogram()

#Histogram of Heavy - Abstainer
ggplot(df_vectors, aes(x = heavy_abstainer)) +
  geom_histogram()

#Percentiles
# ICE Moderate - Abstainer
quantile(sort(df_vectors$moderate_abstainer), probs=c(0.025,0.5, 0.975))
# ICE Heavy - Abstainer
quantile(sort(df_vectors$heavy_abstainer), probs=c(0.025,0.5, 0.975))




#Manual Rubins Rules
#Going to use Rubins Rules to estimate the average causal effect of alcohol consumption
df_vectors <- complete(mice_data, action = "long", include = FALSE)
df_vectors$moderate_abstainer <- df_vectors$Moderate - df_vectors$Abstainer
df_vectors$heavy_abstainer <- df_vectors$Heavy - df_vectors$Abstainer
# 1. Pool the Estimates (Q_bar) and Within-Imputation Variance (W)
stats_per_imp <- df_vectors %>%
  group_by(.imp) %>%
  summarise(
    m_mod = mean(moderate_abstainer, na.rm = TRUE),
    m_hvy = mean(heavy_abstainer, na.rm = TRUE),
    # Variance of the mean (SE^2)
    v_mod = var(moderate_abstainer, na.rm = TRUE) / n(),
    v_hvy = var(heavy_abstainer, na.rm = TRUE) / n(),
    .groups = "drop"
  )

m <- nrow(stats_per_imp) # Number of imputations

# Combined Estimates (Q_bar)
est_mod <- mean(stats_per_imp$m_mod)
est_hvy <- mean(stats_per_imp$m_hvy)

# Within-Imputation Variance (W)
W_mod <- mean(stats_per_imp$v_mod)
W_hvy <- mean(stats_per_imp$v_hvy)

# 2. Between-Imputation Variance (B)
B_mod <- var(stats_per_imp$m_mod)
B_hvy <- var(stats_per_imp$m_hvy)

# 3. Total Variance (T)
# Formula: T = W + B + B/m
T_mod <- W_mod + B_mod + (B_mod / m)
T_hvy <- W_hvy + B_hvy + (B_hvy / m)

# 4. Degrees of Freedom (nu) - Barnard-Rubin adjustment
# Using the formula you provided (which is the standard version)
nu_mod <- (m - 1) * (1 + (W_mod / ((1 + 1/m) * B_mod)))^2
nu_hvy <- (m - 1) * (1 + (W_hvy / ((1 + 1/m) * B_hvy)))^2

# 5. Result Intervals
mod_ci <- est_mod + qt(c(0.025, 0.5, 0.975), df = nu_mod) * sqrt(T_mod)
hvy_ci <- est_hvy + qt(c(0.025, 0.5, 0.975), df = nu_hvy) * sqrt(T_hvy)

#Rubins Moderate - Abstainer
mod_ci

#Rubins Heavy - Abstainer
hvy_ci