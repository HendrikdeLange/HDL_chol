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

