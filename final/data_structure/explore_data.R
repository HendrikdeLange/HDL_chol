library(ggplot2)
library(corrplot)
df <- readRDS("data/nhanes_subset.rds")
colnames(df)

cont_df <- df[,sapply(df, is.numeric)]
factor_df <- df[,sapply(df, is.factor)]

corrplot(cor(cont_df), method = "number")
#BPSysAve and Age have a postive correlation of 0.41
#DirectChol and BMI have a negative correlation of -0.36


cont_vars <- c(colnames(cont_df))
factor_vars <- c(colnames(factor_df))
par(mfrow = c(2, 2))
for (i in 1:length(cont_vars)) {
  for (j in 1:length(factor_vars)) {
    boxplot(df[[cont_vars[i]]] ~ df[[factor_vars[j]]], 
            main = paste(cont_vars[i], "by", factor_vars[j]),
            xlab = factor_vars[j],
            ylab = cont_vars[i]) 
  }
}

#chisq.test(table(data$AlcoholConsumption, data$SmokeNow))
chisq.test(table(df$AlcoholCategory, df$SmokeNow))

abstainer_df <- subset(df, AlcoholCategory = "Abstainer")
ggplot(abstainer_df, aes(x = DirectChol, fill = SmokeNow)) +
  geom_density(alpha = 0.4)

moderate_df <- subset(df, AlcoholCategory == "Moderate")
ggplot(moderate_df, aes(x = DirectChol, fill= SmokeNow)) +
  geom_density(alpha = 0.4)


heavy_df <- subset(df, AlcoholCategory == "Heavy")
ggplot(heavy_df, aes(x = DirectChol, fill = SmokeNow)) +
  geom_density(alpha = 0.4)