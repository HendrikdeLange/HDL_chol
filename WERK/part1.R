#Part 1
df <- readRDS("data/nhanes_subset.rds")
colnames(df)

#Part1
#Model 1 (AOV)
model_1 <- aov(DirectChol ~factor(AlcoholCategory), data=df)
summary(model_1)

#post-hoc test
TukeyHSD(model_1)

#Model 2 (KW)
kruskal.test(DirectChol ~factor(AlcoholCategory), data=df)
#test which test is better
by(df$DirectChol, df$AlcoholCategory, shapiro.test)

#Model 3 (Multi-Way ANOVA)

library(dplyr)
library(ggplot2)

# ── Step 1: Check raw cell sizes ──────────────────────────────────────────────
table(df$AlcoholCategory, df$Gender, df$PhysActive)

# ── Step 2: Balance across all 12 cells ───────────────────────────────────────
min_n <- df %>%
  group_by(AlcoholCategory, Gender, PhysActive) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>%
  min()

cat("Subsampling", min_n, "per cell —", min_n * 12, "total observations\n")

set.seed(123) # for reproducibility
balanced <- df %>%
  group_by(AlcoholCategory, Gender, PhysActive) %>%
  slice_sample(n = min_n) %>%
  ungroup()

# Verify balance
table(balanced$AlcoholCategory, balanced$Gender, balanced$PhysActive)

# ── Step 3: Three-way ANOVA ────────────────────────────────────────────────────
model <- aov(DirectChol ~ factor(AlcoholCategory) * 
               factor(Gender) * 
               factor(PhysActive), 
             data = balanced)
summary(model)

# ── Step 4: Effect size (eta squared) ─────────────────────────────────────────
library(effectsize)
eta_squared(model)

# ── Step 5: Mean plot ─────────────────────────────────────────────────────────
library(jtools)
balanced %>%
  group_by(AlcoholCategory, Gender, PhysActive) %>%
  summarise(mean_chol = mean(DirectChol), .groups = "drop") %>%
  ggplot(aes(x = AlcoholCategory,
             y = mean_chol,
             color = Gender,
             group = Gender)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ PhysActive, 
             labeller = labeller(PhysActive = c("No" = "Not Physically Active",
                                                "Yes" = "Physically Active"))) +
  labs(x = "Alcohol Consumption Category",
       y = "Mean Direct HDL Cholesterol (mmol/L)",
       color = "Gender") +   # no title here
  theme_apa(legend.use.title = TRUE)
ggsave("figure1.png", width = 6.5, height = 4, dpi = 300, units = "in")