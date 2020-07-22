setwd("~/Documents/Study/LaTrobe/Research/phD/Milestone3/R_M3")
water <- read.csv("water.csv", header = TRUE, row.names=NULL)
str(water)
#order factors
water$Treatment <- factor(water$Treatment, levels=c("Control","AutoclavedControl", "Lime", "PoultryLitter", "Biochar", 
                                                          "CitricAcid", "Inositol", "SodiumFormate", "SodiumAcetate",
                                                          "FumaricAcid", "H2O2"))



######################### boxplots WITH  JITTER #########################
require(ggplot2)
p <- ggplot(T16water, aes(x = Treatment, y = T16vgravi.H20)) +
     geom_boxplot(outlier.shape = NA) + geom_jitter() +
    theme_classic()
p

######################### stripchart #########################
stripchart(T16vgravi.H20 ~ Treatment, data = T16water, method = 'jitter' )



########## bargraph using SCIPLOT ##########

library(sciplot)
bargraph.CI(Treatment, T16vgravi.H20, data=T16water, legend=FALSE,
            ylab="Watercontent (gravi)", main = "Water by Treatments")


######################### ANOVA #########################
lm <- lm(T16vgravi.H20 ~ Treatment, data = T16water)
summary(lm) #no signifant difference
mean(T16water$T16vgravi.H20) #0.24
summary(T16water$T16vgravi.H20)



########################## GGPUBR #########################
require(ggpubr)
str(water)
boxp <- ggboxplot(water, x = "Treatment", y = c("T0_H20"),  add = "jitter",legend = "none",
                      outlier.shape = NA,
                      xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.34, method = "anova")+  # Add global annova p-value
 #geom_hline(yintercept = mean(T16water[T16water$Treatment %in% c("Control"),c("T16vgravi.H20")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                   ref.group = "Control", hide.ns = TRUE)
boxp

require(Rmisc)
avgT0 <- summarySE(water[c(-1),], measurevar="T0_H20", groupvars=c("Treatment")) #getting the standard error
avgT16 <- summarySE(water[c(-1),], measurevar="T16_H20", groupvars=c("Treatment")) #getting the standard error

require(reshape2)
water_melt <- melt(water, measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)
avg <- summarySE(water_melt[c(-1),], measurevar="value", groupvars=c("variable"))
avg
