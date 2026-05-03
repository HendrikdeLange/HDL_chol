#Part 2
df <- readRDS("data/nhanes_subset.rds")
colnames(df)

#Model 4 (LR with model 3 factors)
model_4 <- lm(DirectChol ~ AlcoholCategory + Gender +PhysActive, data = df)
summary(model_4)

#Model 5 (LR with all covariates)
model_5 <- lm(DirectChol ~ AlcoholCategory + Age + Gender + BMI + Poverty +
Education + BPSysAve + Diabetes + SmokeNow + PhysActive, data = df)
summary(model_5)

#Model 6 (Objective Prior Bayesian Regression)

