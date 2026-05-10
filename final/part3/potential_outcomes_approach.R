library(mice)
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
mice_data <- mice(imp_data, m = 50, maxit = 10, method = "norm", seed = 123) #Keeps suggesting

#Analysis
library(ggplot2)
library(patchwork)

trt_levels <- c("Abstainer", "Moderate", "Heavy")
M          <- mice_data$m   # 50

#Density Plots
plot_list <- list()

for (treat in trt_levels) {
  obs_vals <- df$DirectChol[df$AlcoholCategory == treat]
  imp_vals <- as.vector(as.matrix(mice_data$imp[[treat]]))
  
  plot_df <- data.frame(
    value  = c(obs_vals, imp_vals),
    source = c(rep("Observed", length(obs_vals)),
               rep("Imputed",  length(imp_vals)))
  )
  
  plot_list[[treat]] <- ggplot(plot_df, aes(x = value, colour = source, fill = source)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_colour_manual(values = c("Observed" = "#2166ac", "Imputed" = "#d6604d")) +
    scale_fill_manual(  values = c("Observed" = "#2166ac", "Imputed" = "#d6604d")) +
    labs(title = treat, x = "DirectChol", y = "Density", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_list[["Abstainer"]] | plot_list[["Moderate"]] | plot_list[["Heavy"]]

#Mean Individual Average Causal Effect
# Build n x 3 matrix: row mean of 50 imputations per missing unit
mean_imp <- matrix(NA, nrow = nrow(df), ncol = length(trt_levels))
colnames(mean_imp) <- trt_levels

for (treat in trt_levels) {
  imp_mat                              <- mice_data$imp[[treat]]
  mean_imp[, treat]                    <- df[[treat]]
  missing_idx                          <- as.numeric(rownames(imp_mat))
  mean_imp[missing_idx, treat]         <- rowMeans(imp_mat)
}

# Individual effects vs Abstainer, then mean and percentiles across units
ace_unit <- data.frame(
  Comparison = character(), Mean_ACE = numeric(),
  Lower_2.5  = numeric(),   Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in c("Moderate", "Heavy")) {
  ind_fx <- mean_imp[, treat] - mean_imp[, "Abstainer"]
  ace_unit <- rbind(ace_unit, data.frame(
    Comparison = paste(treat, "- Abstainer"),
    Mean_ACE   = round(mean(ind_fx),             3),
    Lower_2.5  = round(quantile(ind_fx, 0.025),  3),
    Upper_97.5 = round(quantile(ind_fx, 0.975),  3)
  ))
}

print(ace_unit)

#Mean Individual Causal Effect
ace_imp <- data.frame(
  Comparison = character(), Mean_ACE = numeric(),
  Lower_2.5  = numeric(),   Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in c("Moderate", "Heavy")) {
  ate_vec <- numeric(M)
  for (m in 1:M) {
    comp        <- complete(mice_data, m)s
    ate_vec[m]  <- mean(comp[[treat]] - comp[["Abstainer"]], na.rm = TRUE)
  }
  ace_imp <- rbind(ace_imp, data.frame(
    Comparison = paste(treat, "- Abstainer"),
    Mean_ACE   = round(mean(ate_vec),             3),
    Lower_2.5  = round(quantile(ate_vec, 0.025),  3),
    Upper_97.5 = round(quantile(ate_vec, 0.975),  3)
  ))
}

print(ace_imp)

#Rubins Rules
rubin <- data.frame(
  Comparison = character(), ATE = numeric(), SE = numeric(),
  Lower_2.5  = numeric(),   Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in c("Moderate", "Heavy")) {
  ate_vec <- numeric(M)
  var_vec <- numeric(M)
  
  for (m in 1:M) {
    comp        <- complete(mice_data, m)
    ind_fx      <- comp[[treat]] - comp[["Abstainer"]]
    n           <- sum(!is.na(ind_fx))
    ate_vec[m]  <- mean(ind_fx, na.rm = TRUE)
    var_vec[m]  <- var(ind_fx,  na.rm = TRUE) / n
  }
  
  Q_bar  <- mean(ate_vec)
  V_W    <- mean(var_vec)
  V_B    <- var(ate_vec)
  V_T    <- V_W + (1 + 1/M) * V_B
  SE     <- sqrt(V_T)
  df_rb  <- (M - 1) * (1 + V_W / ((1 + 1/M) * V_B))^2
  t_crit <- qt(0.975, df = df_rb)
  
  rubin <- rbind(rubin, data.frame(
    Comparison = paste(treat, "- Abstainer"),
    ATE        = round(Q_bar,                3),
    SE         = round(SE,                   3),
    Lower_2.5  = round(Q_bar - t_crit * SE, 3),
    Upper_97.5 = round(Q_bar + t_crit * SE, 3)
  ))
}

print(rubin)

#Combining Everything
all_intervals <- rbind(
  data.frame(Comparison = ace_unit$Comparison, ATE = ace_unit$Mean_ACE,
             Lower = ace_unit$Lower_2.5, Upper = ace_unit$Upper_97.5,
             Method = "Unit percentile"),
  data.frame(Comparison = ace_imp$Comparison,  ATE = ace_imp$Mean_ACE,
             Lower = ace_imp$Lower_2.5,  Upper = ace_imp$Upper_97.5,
             Method = "Imputation percentile"),
  data.frame(Comparison = rubin$Comparison,    ATE = rubin$ATE,
             Lower = rubin$Lower_2.5,    Upper = rubin$Upper_97.5,
             Method = "Rubin's Rules")
)

ggplot(all_intervals, aes(x = ATE, y = Comparison, colour = Method, shape = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Unit percentile"        = "#2166ac",
                                 "Imputation percentile"  = "#4dac26",
                                 "Rubin's Rules"          = "#d6604d")) +
  scale_shape_manual(values  = c("Unit percentile"        = 16,
                                 "Imputation percentile"  = 17,
                                 "Rubin's Rules"          = 15)) +
  labs(title  = "95% intervals for ATE vs Abstainer",
       x      = "ATE (DirectChol, mmol/L)",
       y      = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")