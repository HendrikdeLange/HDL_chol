library(ggplot2)
library(dplyr)
library(car)
df <- readRDS("data/nhanes_subset.rds")
colnames(df)
df %>% count(df$AlcoholCategory)
df$DirectChol <- exp(df$DirectChol)


#1 Continous Dependent Variable -CHECK
#2 Independence of Observatiions -CHECK (How the data was collected)
#Normality of Residuals
model_1 <- aov(DirectChol ~factor(AlcoholCategory), data=df)
residuals <- residuals(model_1, breaks = 30)
shapiro.test(residuals)
par(mfrow=c(2,2))
hist(residuals, main = "Histogram of Residuals Before Log")
qqnorm(residuals, main = "Histogram of Residuals Before Log")
qqline(residuals)


#RESIDULS ARE NOT NORMAL BUT ANOVA IS FAIRLY ROBUST IF N > 30 


#If we log DirectChol, the residuals are much more normal!!!!!!!!!!!!
df$DirectChol <- log(df$DirectChol)
model_1 <- aov(DirectChol ~factor(AlcoholCategory), data=df)
residuals <- residuals(model_1, breaks = 30)
shapiro.test(residuals)
hist(residuals,  main = "Histogram of Residuals After Log")
qqnorm(residuals, main = "Histogram of Residuals After Log")
qqline(residuals)


#Homogeneity of Variance
#Assumption is that variance amoung groups must be equal
print(leveneTest(DirectChol ~ factor(AlcoholCategory), data = df))
