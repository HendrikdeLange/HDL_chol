library(mice)
df <- read.csv("C:/Users/hendr/Downloads/fluoride.csv")

# Step 1 — treatment setup
obs_num    <- nrow(df)
df$trt     <- factor(df$treatmt, levels = c("W", "APF", "SF"))
dummies    <- model.matrix(~ df$trt - 1)
dummies[dummies == 0] <- NA
colnames(dummies) <- c("W", "APF", "SF")

# Step 2 — multiply responses
responses <- as.vector(df$after) * dummies
df        <- cbind(df, responses)

# Step 3 — build imputation dataset
data <- df[, c("age", "inst", "before", "W", "APF", "SF")]  

# Step 4 — MICE
mice_data <- mice(data, m = 50, maxit = 10, method = 'norm', seed = 10) 

#Bullet 1
library(mice)
library(ggplot2)
library(tidyr)
library(dplyr)

trt_levels <- c("W", "APF", "SF")

plot_list <- list()

for (treat in trt_levels) {
  
  # --- Observed values for this arm ---
  obs_vals <- df$after[df$trt == treat]
  
  # --- All 50 imputed draws for this arm ---
  imp_mat  <- mice_data$imp[[treat]]          # rows = missing units, cols = imputations
  imp_vals <- as.vector(as.matrix(imp_mat))   # stack all draws into one vector
  
  # --- Combine into one data frame ---
  plot_df <- data.frame(
    value  = c(obs_vals, imp_vals),
    source = c(rep("Observed",  length(obs_vals)),
               rep("Imputed",   length(imp_vals)))
  )
  
  plot_list[[treat]] <- ggplot(plot_df, aes(x = value, colour = source, fill = source)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_colour_manual(values = c("Observed" = "#2166ac", "Imputed" = "#d6604d")) +
    scale_fill_manual(  values = c("Observed" = "#2166ac", "Imputed" = "#d6604d")) +
    labs(
      title    = paste("Treatment:", treat),
      x        = "Response (after)",
      y        = "Density",
      colour   = NULL,
      fill     = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

# --- Print all three plots ---
library(patchwork)
plot_list[["W"]] | plot_list[["APF"]] | plot_list[["SF"]]




#BULLET 2
# --- Extract mean imputed value per unit per treatment arm ---
mean_imp <- sapply(trt_levels, function(treat) {
  
  imp_mat      <- mice_data$imp[[treat]]        # missing units × 50 imputations
  observed     <- df$after[df$trt == treat]     # observed values for this arm
  
  # Start with the full column of NAs from the original data
  full_vec           <- df$after
  full_vec[]         <- NA
  
  # Fill in observed values
  full_vec[df$trt == treat] <- observed
  
  # Fill in row means of imputed draws for missing units
  missing_idx        <- as.numeric(rownames(imp_mat))
  full_vec[missing_idx] <- rowMeans(imp_mat)
  
  full_vec
})

# mean_imp is now an obs_num × 3 matrix, one column per treatment arm
colnames(mean_imp) <- trt_levels

# --- Individual causal effects vs reference arm (W) ---
ace_output <- data.frame(
  Comparison = character(),
  Mean_ACE   = numeric(),
  Lower_2.5  = numeric(),
  Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in trt_levels[-1]) {      # APF and SF vs W
  
  ind_effects <- mean_imp[, treat] - mean_imp[, "W"]   # one effect per unit
  
  ace_output <- rbind(ace_output, data.frame(
    Comparison = paste(treat, "- W"),
    Mean_ACE   = round(mean(ind_effects),              3),
    Lower_2.5  = round(quantile(ind_effects, 0.025),  3),
    Upper_97.5 = round(quantile(ind_effects, 0.975),  3)
  ))
}

print(ace_output)

#bullet 3
M          <- mice_data$m       # 50
trt_levels <- c("W", "APF", "SF")

ace_output <- data.frame(
  Comparison = character(),
  Mean_ACE   = numeric(),
  Lower_2.5  = numeric(),
  Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in trt_levels[-1]) {   # APF vs W, SF vs W
  
  ate_vec <- numeric(M)
  
  for (m in 1:M) {
    
    # Complete dataset for imputation m
    comp <- complete(mice_data, m)
    
    # Individual causal effect per unit
    ind_effects <- comp[[treat]] - comp[["W"]]
    
    # ATE for this imputed dataset
    ate_vec[m] <- mean(ind_effects, na.rm = TRUE)
  }
  
  ace_output <- rbind(ace_output, data.frame(
    Comparison = paste(treat, "- W"),
    Mean_ACE   = round(mean(ate_vec),             3),
    Lower_2.5  = round(quantile(ate_vec, 0.025),  3),
    Upper_97.5 = round(quantile(ate_vec, 0.975),  3)
  ))
}

print(ace_output)

#Bullet 3
M          <- mice_data$m
trt_levels <- c("W", "APF", "SF")

rubin_output <- data.frame(
  Comparison = character(),
  ATE        = numeric(),
  SE         = numeric(),
  Lower_2.5  = numeric(),
  Upper_97.5 = numeric(),
  stringsAsFactors = FALSE
)

for (treat in trt_levels[-1]) {
  
  ate_vec <- numeric(M)
  var_vec <- numeric(M)
  
  for (m in 1:M) {
    
    comp         <- complete(mice_data, m)
    ind_effects  <- comp[[treat]] - comp[["W"]]
    n            <- sum(!is.na(ind_effects))
    
    # ATE and SE for this imputed dataset
    ate_vec[m]  <- mean(ind_effects, na.rm = TRUE)
    var_vec[m]  <- var(ind_effects,  na.rm = TRUE) / n    # SE² = sample var / n
  }
  
  # --- Rubin's Rules ---
  
  # Pooled ATE
  Q_bar <- mean(ate_vec)
  
  # Within-imputation variance (average of the M variances)
  V_W   <- mean(var_vec)
  
  # Between-imputation variance
  V_B   <- var(ate_vec)
  
  # Total variance
  V_T   <- V_W + (1 + 1/M) * V_B
  
  SE    <- sqrt(V_T)
  
  # Degrees of freedom (Barnard & Rubin 1999)
  df    <- (M - 1) * (1 + V_W / ((1 + 1/M) * V_B))^2
  
  t_crit <- qt(0.975, df = df)
  
  rubin_output <- rbind(rubin_output, data.frame(
    Comparison = paste(treat, "- W"),
    ATE        = round(Q_bar,           3),
    SE         = round(SE,              3),
    Lower_2.5  = round(Q_bar - t_crit * SE, 3),
    Upper_97.5 = round(Q_bar + t_crit * SE, 3)
  ))
}

print(rubin_output)

#Combine everything
library(ggplot2)

# --- Combine all interval estimates into one data frame ---

# From the percentile approach (previous question before Rubin's)
percentile_df <- data.frame(
  Comparison = ace_output$Comparison,
  ATE        = ace_output$Mean_ACE,
  Lower      = ace_output$Lower_2.5,
  Upper      = ace_output$Upper_97.5,
  Method     = "Percentile"
)

# From Rubin's Rules
rubin_df <- data.frame(
  Comparison = rubin_output$Comparison,
  ATE        = rubin_output$ATE,
  Lower      = rubin_output$Lower_2.5,
  Upper      = rubin_output$Upper_97.5,
  Method     = "Rubin's Rules"
)

plot_df <- rbind(percentile_df, rubin_df)

# --- Forest plot ---
ggplot(plot_df, aes(x = ATE, y = Comparison, colour = Method, shape = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height   = 0.15,
                 linewidth = 0.8,
                 position = position_dodge(width = 0.4)) +
  geom_point(size     = 3,
             position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Percentile"   = "#2166ac",
                                 "Rubin's Rules" = "#d6604d")) +
  scale_shape_manual(values  = c("Percentile"   = 16,
                                 "Rubin's Rules" = 17)) +
  labs(
    title   = "95% intervals for average treatment effects",
    x       = "ATE vs W (reference)",
    y       = NULL,
    colour  = "Method",
    shape   = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )