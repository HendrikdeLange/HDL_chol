set.seed(123)
par(mfrow=c(2,2))
#Part 2
df <- readRDS("data/nhanes_subset.rds")
colnames(df)

#Model 4 (LR with model 3 factors)
model_4 <- lm(DirectChol ~ AlcoholCategory + Age + Gender + BMI, data = df)
summary(model_4)


#Model 5 (LR with all covariates) (NON-LOGGED)
df$DirectChol_exp <- exp(df$DirectChol)
model_5 <- lm(DirectChol_exp ~ AlcoholCategory + Age + Gender + BMI 
              + BPSysAve + Diabetes + SmokeNow + PhysActive, data = df)
summary(model_5)

residuals <- residuals(model_5)

hist(residuals, breaks = 15, main="Histogram of Residuals (Non-Logged HDL)", xlab = "Residuals")

qqnorm(residuals, main ="QQ (Non-Logged HDL)")
qqline(residuals)

#Model 5 (LR with all covariates)

model_5 <- lm(DirectChol ~ AlcoholCategory + Age + Gender + BMI 
              + BPSysAve + Diabetes + SmokeNow + PhysActive, data = df)
summary(model_5)

residuals <- residuals(model_5)

hist(residuals, breaks = 15, main = "Histogram of Residuals (Logged HDL)",xlab = "Residuals")

qqnorm(residuals, main ="QQ (Logged HDL)")
qqline(residuals)



