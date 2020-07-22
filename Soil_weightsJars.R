setwd("~/Documents/Study/LaTrobe/Research/Methods/Experiment1")

data <- read.csv("Soil_weightsInJars.csv", header = TRUE, row.names=NULL)

require(rafalib)
mypar(3,1)
boxplot(T6~Treatment,data=data) 
aovfc <- aov(T2 ~ Treatment, data = data)
summary(aovfc)
aovfc


#another way with lm function
lmSpec <- lm(T2 ~ Treatment, data = data)
summary(lmSpec)
#library(multcomp)
#tukey_spec <- glht(lmSpec, linfct=mcp(group1="Tukey"))
#summary(tukey_spec)