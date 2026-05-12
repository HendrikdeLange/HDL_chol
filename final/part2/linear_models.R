set.seed(123)
#Part 2
df <- readRDS("data/nhanes_subset.rds")
colnames(df)


#Model 4 (LR with model 3 factors)
model_4 <- lm(DirectChol ~ AlcoholCategory + Age + Gender + BMI, data = df)
summary(model_4)

#Model 5 (LR with all covariates)
model_5 <- lm(DirectChol ~ AlcoholCategory + Age + Gender + BMI 
              + BPSysAve + Diabetes + SmokeNow + PhysActive, data = df)
summary(model_5)

# QQ Plot of residuals for Model 5
qqnorm(residuals(model_5), 
       main = "Normal Q-Q Plot of Residuals - Model 5",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       pch  = 1,
       col  = "steelblue")
qqline(residuals(model_5), col = "red", lwd = 2)

# Histogram of Residuals for Model 5
hist(residuals(model_5),
     main = "Histogram of Residuals - Model 5",
     xlab = "Residuals",
     ylab = "Frequency",
     col  = "steelblue",
     border = "white",
     breaks = 30)

# Add a normal curve overlay
x <- seq(min(residuals(model_5)), max(residuals(model_5)), length = 200)
curve_vals <- dnorm(x, mean = mean(residuals(model_5)), sd = sd(residuals(model_5)))
lines(x, curve_vals * length(residuals(model_5)) * diff(hist(residuals(model_5), 
                                                             breaks = 30, plot = FALSE)$breaks)[1], 
      col = "red", lwd = 2)

legend("topright", legend = "Normal curve", col = "red", lwd = 2)
