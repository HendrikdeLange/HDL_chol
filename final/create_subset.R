library(NHANES)
library(dplyr)

# Constructing the subset
nhanes_subset <- NHANES %>%
  # 1. Recode Treatment as discussed
  mutate(AlcoholDay = ifelse(is.na(AlcoholDay), 0, AlcoholDay)) %>%
  mutate(AlcoholCategory = case_when(
    AlcoholDay == 0 ~ "Abstainer",
    AlcoholDay >= 1 & AlcoholDay <= 2 ~ "Moderate",
    AlcoholDay >= 3 ~ "Heavy"
  )) %>%
  # 2. Select Response, Treatment, and 7 Covariates
  select(
    DirectChol,      # Response
    AlcoholCategory, # Treatment
    Age, Gender, BMI, BPSysAve, Diabetes, SmokeNow, PhysActive
  ) %>%
  # 3. Handle missing values (important for Bayesian/Matrix-based models)

  na.omit()

nhanes_subset$DirectChol <- log(nhanes_subset$DirectChol)
# Check dimensions
dim(nhanes_subset)
summary(nhanes_subset)
n <- length(nhanes_subset$Age)#2946
n

#Save Data
saveRDS(nhanes_subset, "data/nhanes_subset.rds")

