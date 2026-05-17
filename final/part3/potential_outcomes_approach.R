library(mice)
library(dplyr)
library(ggplot2)
set.seed(123)

# --- Load data ---
df <- readRDS("data/nhanes_subset.rds")

min(df$DirectChol, na.rm = TRUE)
sum(df$DirectChol < 0, na.rm = TRUE)
# --- Step 1: Treatment setup ---
df$AlcoholCategory <- relevel(factor(df$AlcoholCategory), ref = "Abstainer")
trt_levels         <- levels(df$AlcoholCategory)   # "Abstainer", "Moderate", "Heavy"





#####PLOTTING PURPOSES
obs_abstainer_vec <- df$DirectChol[df$AlcoholCategory == "Abstainer"]
obs_moderate_vec <- df$DirectChol[df$AlcoholCategory == "Moderate"]
obs_heavy_vec <- df$DirectChol[df$AlcoholCategory == "Heavy"]
#############

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
mice_data <- mice(imp_data, m = 200, maxit = 20, method = "norm", seed = 100) 




#TESTING
#this is mice
df_vectors <- complete(mice_data, action = "long", include = FALSE) 
df_vectors$moderate_abs <- exp(df_vectors$Moderate) - exp(df_vectors$Abstainer)
df_vectors$heavy_abs <- exp(df_vectors$Heavy) - exp(df_vectors$Abstainer)

collapsed <- df_vectors %>%
  group_by(.imp) %>%                        
  summarise(
    moderate_abs = mean(moderate_abs, na.rm = TRUE),
    heavy_abs    = mean(heavy_abs,     na.rm = TRUE),
    .groups = "drop"
  )


#Percentiles
# ICE Moderate - Abstainer
quantile(sort(collapsed$moderate_abs), probs=c(0.025,0.5, 0.975))
# ICE Heavy - Abstainer
quantile(sort(collapsed$heavy_abs), probs=c(0.025,0.5, 0.975))


#Bullet 3
#Use entire data set again
df_vectors <- complete(mice_data, action = "long", include = FALSE)
df_vectors$moderate_abstainer <- exp(df_vectors$Moderate) - exp(df_vectors$Abstainer)
df_vectors$heavy_abstainer <- exp(df_vectors$Heavy) - exp(df_vectors$Abstainer)


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
df_vectors$moderate_abstainer <- exp(df_vectors$Moderate) - exp(df_vectors$Abstainer)
df_vectors$heavy_abstainer <- exp(df_vectors$Heavy) - exp(df_vectors$Abstainer)

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

sqrt(T_mod)
sqrt(T_hvy)

#Making the plots myself
df_vectors <- complete(mice_data, action = "long", include = FALSE)
df_vectors$Abstainer <- df_vectors$Abstainer
df_vectors$Moderate <- df_vectors$Moderate
df_vectors$Heavy <- df_vectors$Heavy
plot_df <- data.frame(
  value = c(df_vectors$Abstainer, df_vectors$Moderate, df_vectors$Heavy,
            obs_abstainer_vec, obs_moderate_vec, obs_heavy_vec),
  
  type = c(rep("Imputed",  length(df_vectors$Abstainer)),
           rep("Imputed",  length(df_vectors$Moderate)),
           rep("Imputed",  length(df_vectors$Heavy)),
           rep("Observed", length(obs_abstainer_vec)),
           rep("Observed", length(obs_moderate_vec)),
           rep("Observed", length(obs_heavy_vec))),
  
  category = c(rep("Abstainer", length(df_vectors$Abstainer)),
               rep("Moderate",  length(df_vectors$Moderate)),
               rep("Heavy",     length(df_vectors$Heavy)),
               rep("Abstainer", length(obs_abstainer_vec)),
               rep("Moderate",  length(obs_moderate_vec)),
               rep("Heavy",     length(obs_heavy_vec)))
)
plot_df$category <- factor(
  plot_df$category,
  levels = c("Abstainer", "Moderate", "Heavy")
)

plot_df$type <- factor(plot_df$type, levels = c("Observed", "Imputed"))

ggplot(plot_df, aes(x = value, fill = type, colour = type)) +
  geom_density(alpha = 0.4, linewidth = 0.9) +
  facet_wrap(~ category) +
  scale_fill_manual(values = c("Imputed" = "#EC058E", "Observed" = "#41521F")) +
  scale_colour_manual(values = c("Imputed" = "#EC058E", "Observed" = "#41521F")) +
  labs(x = "log(DirectChol (mmol/L))", y = "Density",
       fill = "Type", colour = "Type") +
  theme_minimal()

test_vector <- df_vectors$Heavy[df_vectors$.imp == 50]
densityplot(test_vector)
