#Part 1
library(ggplot2)
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




#Graphs for the treatment
ggplot(df, aes(x = DirectChol)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  facet_wrap(~ AlcoholCategory)
 
#KW is better
