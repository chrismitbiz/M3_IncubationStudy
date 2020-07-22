setwd("~/Documents/Study/LaTrobe/Research/phD/Milestone3/R_M3")

#packages
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/rstatix")

library(agricolae)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(nmle)
library(lme4)

##### loading all data ######
residues <- read.csv("Residues.csv", header = TRUE, row.names=NULL)
residues_na <- na.omit(residues)

residues <- residues[c(-2, -23,-31),]  #? minus the unusually low values in sodium formate (31), and poultry litter (23)
#should probably re-extract those 2 samples to confirm
residues <- residues[c(-2, -23),] 
#also sample#2 spilled during extraction so values are taken out as well. 

residues <- residues[c(-45,-46,-47,-48,-49,-50, -51),] 
#excluding airdried control
str(residues)
#order factors
residues$Treatment <- factor(residues$Treatment, levels=c("Control","AutoclavedControl", "Lime", "PoultryLitter", "Biochar", 
                                                          "CitricAcid", "Inositol", "SodiumFormate", "SodiumAcetate",
                                                          "FumaricAcid", "H2O2", "AirdriedSoil"))

#remove H202 and autoclaved control
residues <- residues[residues$Treatment %in% c("Control", "Lime", "PoultryLitter", "Biochar", 
                                               "CitricAcid", "Inositol", "SodiumFormate", "SodiumAcetate",
                                               "FumaricAcid"),]  
#order factors
residues$Treatment <- factor(residues$Treatment, levels=c("Control", "Lime", "PoultryLitter", "Biochar", 
                                                          "CitricAcid", "Inositol", "FumaricAcid", "SodiumFormate", "SodiumAcetate"))

##### loading trtms related to soil ag amendments + H202######
str(residues)
# subsetting based on variable values
residues_AG <- residues[residues$Treatment %in% c("Control", "Lime", "PoultryLitter", "Biochar", "H2O2","AutoclavedControl"),]          
residues_AG$Treatment <- factor(residues_AG$Treatment, levels=c("Control", "Lime", "PoultryLitter", "Biochar", "H2O2","AutoclavedControl"))

##### loading trtms related to carbon substrates ######
residues_C <- residues[residues$Treatment %in% c("Control","CitricAcid", "Inositol","FumaricAcid"),] 
residues_C$Treatment <- factor(residues_C$Treatment, levels=c("Control","CitricAcid", "Inositol", "FumaricAcid"))

##### loading trtms related to electron transport  ######
residues_E <- residues[residues$Treatment %in% c("Control", "SodiumFormate","SodiumAcetate"),] 
residues_E$Treatment <- factor(residues_E$Treatment, levels=c("Control", "SodiumFormate","SodiumAcetate"))

##### loading trtms related to electron transport  ######
residues <- residues[residues$Treatment %in% c("Control",  "CitricAcid", "Inositol","FumaricAcid", "SodiumFormate","SodiumAcetate"),] 
residues$Treatment <- factor(residues_C_E$Treatment, levels=c("Control", "CitricAcid", "Inositol","FumaricAcid", "SodiumFormate","SodiumAcetate"))

##### loading autoclaved and H202 treatments  ######
residues_Other <- residues[residues$Treatment %in% c("Control", "AutoclavedControl","H2O2"),] 
residues_Other$Treatment <- factor(residues_Other$Treatment, levels=c("Control", "AutoclavedControl","H2O2"))

residues$Treatment

####### basic boxplot ######
par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.1)
boxplot(Dieldrin ~ Treatment, data = residues, main = "Dieldrin")
boxplot(DDE ~ Treatment, data = residues, main = "DDE")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
boxplot(Dieldrin ~ Treatment, data = residues[c(-22,-30),], main = "Dieldrin") #without the unusually low values in sodium formate (31), and poultry litter (23)
boxplot(DDE ~ Treatment, data = residues, main = "DDE")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default



library(tidyverse)


###### GGPUBR box plot ##### 
library(ggpubr)
#mean of control, dieldrin
str(residues)
mean(residues[residues$Treatment %in% c("Control"),c("Diel_DStd_IncFront")]) #1.218


##### Carbon substrates #####
#Dieldrin
str(residues_C)
boxp1D <- ggboxplot(residues_C, x = "Treatment", y = "Diel_DStd_IncFront",  add = "jitter",legend = "none",
          color = "Treatment", palette = "uchicago",
          xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
                    ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of carbon substrates on Dieldrin")
boxp1D <- boxp1D + scale_y_continuous(breaks=seq(0,1.35,0.05))

boxp1D <- ggboxplot(residues_C, x = "Treatment", y = "Diel_DStd_Middle",  add = "jitter",legend = "none",
                    color = "Treatment", palette = "uchicago",
                    xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of carbon substrates on Dieldrin")
boxp1D <- boxp1D + scale_y_continuous(breaks=seq(0,1.35,0.05))

# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C, method = "t.test")
?compare_means
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid") )
#ggboxplot(residues_C, x = "Treatment", y = "Dieldrin",  add = "jitter",legend = "none",
#          outlier.shape = NA, color = "TreatmentGroup2")+
          #palette =c("black","darkgrey")) 
#  rotate_x_text(angle = 45)+ 
#  stat_compare_means(label.y = 1.5, method = "anova")+  # Add global annova p-value
#  geom_hline(yintercept = 1.218, linetype = 2)+
#  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#  ggtitle("Carbon substrates") # Add horizontal line at base mean
  #stat_compare_means(label = "p.signif", method = "t.test",   # Pairwise comparison against all
  #                   ref.group = ".all.", hide.ns = TRUE)


#DDE
  mean(residues[residues$Treatment %in% c("Control"),c("DDE")]) #0.31033
  
boxp1_DDE <- ggboxplot(residues_C, x = "Treatment", y = "DDE",  add = "jitter",legend = "none",
            color = "TreatmentGroup2", palette = "uchicago",
            xlab = FALSE) +
    rotate_x_text(angle = 45)+ 
    stat_compare_means(label.y = 0.32, method = "anova")+  # Add global annova p-value
    geom_hline(yintercept = 0.31033, linetype = 2)+ # Add horizontal line at base mean
    stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
                       ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of carbon substrates on DDE")
boxp1_DDE<- boxp1_DDE + scale_y_continuous(breaks=seq(0,0.32,0.005))

############### ggarrange D + DDE #########
ggarrange(boxp1D, boxp1_DDE, legend = "none", common.legend = TRUE, nrow = 2, ncol = 1, labels = "auto")







############ melt data to compare means of T0 and 16 #############
require(reshape2) #to get the table into the right format 
require(ggpubr)
residues_Donly <- residues[c(1:3,5:6)]
residues_melt <- melt(residues_Donly, measure.vars = c(3:5), na.rm = TRUE) #measure.vars = msr)
residues_melt <- melt(residues_C_E, measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)
#Dieldrin only
residues_melt <- melt(residues[c(1:2,5:6)], measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)
#DDE only
residues_melt <- melt(residues[c(1:4)], measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)
#env_melt$Region <- factor(env_melt$Region, levels = c("KV", "GE")) #changing the order of Region factor
str(residues_melt)
str(residues)



ggboxplot(residues_melt, x = "variable", y = "value",  add = "jitter",legend = "none",
                    color = "variable", palette = "uchicago", facet.by = "Treatment",
                    xlab = FALSE,   short.panel.labs = FALSE) + stat_compare_means(label = "p.format")
                                                                                
  
?ggboxplot


##### Carbon PLUS Energy substrates #####
#Dieldrin T16
boxp2D_0 <- ggboxplot(residues, x = "Treatment", y = "T16_Diel_DStd",  add = "jitter",legend = "none",
                     palette = "uchicago",outlier.shape = NA,
                    xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.3, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T16_Diel_DStd")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Residues after 16 weeks") + scale_y_continuous(name = "Dieldrin µg / g soil")
boxp2D_0 


#Dieldrin T0
boxp2D_T0 <- ggboxplot(residues_C_E, x = "Treatment", y = "T0_Diel_DStd",  add = "jitter",legend = "none",
                      color = "Treatment", palette = "uchicago",outlier.shape = NA,
                      xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.3, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T0_Diel_DStd")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Residues at beginning (T0)") + scale_y_continuous(name = "Dieldrin µg / g soil")
boxp2D_T0 

#Dieldrin T0 _all treatments
boxp2D_T0 <- ggboxplot(residues, x = "Treatment", y = "T0_Diel_DStd",  add = "jitter",legend = "none",
                       color = "Treatment", palette = "uchicago",outlier.shape = NA,
                       xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T0_Diel_DStd")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Residues at beginning (T0)") + scale_y_continuous(name = "Dieldrin µg / g soil")

boxp2D_T0 



#### plotting the residue loss  ####
str(residues_C_E)
residueloss <- (residues_C_E$T0_Dieldrin - residues_C_E$T16_Dieldrin) / residues_C_E$T0_Dieldrin 
residues_C_E <- cbind(residues_C_E, residueloss)
boxp2D_1 <- ggboxplot(residues_C_E, x = "Treatment", y = "residueloss",  add = "jitter",legend = "none",
                    color = "Treatment", palette = "uchicago",outlier.shape = NA,
                    xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.7, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("residueloss")]), linetype = 2)+ # Add horizontal line at base mean
  #stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
  #                   ref.group = "Control", hide.ns = TRUE)+
  ggtitle("% residue loss to control after 16 weeks") 
boxp2D_1 <- boxp2D_1 + scale_y_continuous(breaks=seq(0,1.4,0.1),name = "Loss (%) after 16 weeks")
boxp2D_1




#AG residues
boxp2D_3 <- ggboxplot(residues_AG, x = "Treatment", y = "T16_Dieldrin",  add = "jitter",legend = "none",
                   color = "Treatment", palette = "uchicago",outlier.shape = NA,
                  xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.3, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T16_Dieldrin")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Residues after 16 weeks (Middle Column)") + scale_y_continuous(name = "Dieldrin µg / g soil")
boxp2D_3



############### ggarrange 3x Dieldrin #########
ggarrange(boxp2D_0, boxp2D_1, boxp2D_3, legend = "none", common.legend = TRUE, nrow = 3, ncol = 1, labels = "auto")

?stat_compare_means
# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C, method = "t.test")
?compare_means
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid") )
#ggboxplot(residues_C, x = "Treatment", y = "Dieldrin",  add = "jitter",legend = "none",
#          outlier.shape = NA, color = "TreatmentGroup2")+
#palette =c("black","darkgrey")) 
#  rotate_x_text(angle = 45)+ 
#  stat_compare_means(label.y = 1.5, method = "anova")+  # Add global annova p-value
#  geom_hline(yintercept = 1.218, linetype = 2)+
#  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#  ggtitle("Carbon substrates") # Add horizontal line at base mean
#stat_compare_means(label = "p.signif", method = "t.test",   # Pairwise comparison against all
#                   ref.group = ".all.", hide.ns = TRUE)


#DDE
str(residues_C_E)
mean(residues[residues$Treatment %in% c("Control"),c("DDE")]) #0.31033

boxp2_DDE <- ggboxplot(residues_C_E, x = "Treatment", y = "DDE_OC4704",  add = "jitter",legend = "none",
                       color = "Treatment", palette = "uchicago",
                       xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.32, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = mean(residues_C_E[residues_C_E$Treatment %in% c("Control"),c("DDE_OC4704")]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of carbon substrates on DDE")
boxp2_DDE<- boxp2_DDE + scale_y_continuous(breaks=seq(0,0.32,0.005), name = "DDE µg / g soil") 

############### ggarrange D + DDE #########
ggarrange(boxp2D_1, boxp2_DDE, legend = "none", common.legend = TRUE, nrow = 2, ncol = 1, labels = "auto")
#dieldrin with values of middle column only





##### Soil Amendments  #####
#Dieldrin
boxp3D <- ggboxplot(residues_AG, x = "Treatment", y = "Dieldrin",  add = "jitter",legend = "none",
                    color = "TreatmentGroup2", palette = "uchicago",
                    xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of soil amendments on dieldrin")
boxp3D <- boxp3D + scale_y_continuous(breaks=seq(0,1.4,0.05),name = "Dieldrin µg / g soil")

?stat_compare_means
# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C, method = "t.test")
?compare_means
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid") )
#ggboxplot(residues_C, x = "Treatment", y = "Dieldrin",  add = "jitter",legend = "none",
#          outlier.shape = NA, color = "TreatmentGroup2")+
#palette =c("black","darkgrey")) 
#  rotate_x_text(angle = 45)+ 
#  stat_compare_means(label.y = 1.5, method = "anova")+  # Add global annova p-value
#  geom_hline(yintercept = 1.218, linetype = 2)+
#  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#  ggtitle("Carbon substrates") # Add horizontal line at base mean
#stat_compare_means(label = "p.signif", method = "t.test",   # Pairwise comparison against all
#                   ref.group = ".all.", hide.ns = TRUE)


#DDE
mean(residues[residues$Treatment %in% c("Control"),c("DDE")]) #0.31033

boxp3_DDE <- ggboxplot(residues_AG, x = "Treatment", y = "DDE",  add = "jitter",legend = "none",
                       color = "TreatmentGroup2", palette = "uchicago",
                       xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.32, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 0.31033, linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  ggtitle("Effect of soil amendments on DDE")
boxp3_DDE<- boxp3_DDE + scale_y_continuous(breaks=seq(0,0.32,0.005), name = "DDE µg / g soil") 

############### ggarrange D + DDE #########
ggarrange(boxp3D, boxp3_DDE, legend = "none", common.legend = TRUE, nrow = 2, ncol = 1, labels = "auto")


############### ggarrange 4 plots #########
ggarrange(boxp2D,  boxp3D, boxp2_DDE, boxp3_DDE, legend = "none", common.legend = TRUE, nrow = 2, ncol = 2, labels = "auto")










###### GGPUBR bar plot with labels ##### 
library(ggplot2)
library(agricolae)

# Bar plot of mean +/-se
#Colour palettes
#ggsci scientific journal palettes, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".

#  stat_compare_means() +      # Global p-value
#  stat_compare_means(ref.group = "0.5", label = "p.signif", label.y = c(22, 29))  # compare to ref.group


#### all treatments ######
ANOVA=aov(Dieldrin ~ Treatment, data = residues)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
comparison = HSD.test(ANOVA, "Paddock", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
target <- c("1", "2", "3", "4", 
            "5", "6", "7", "8",
            "9", "10", "11", "12")
groups1_ordered <- groups1[match(target, row.names(groups1)),]

str(residues)
# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C)
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid")
#                       , c("CitricAcid", "FumaricAcid"), c("Inositol", "FumaricAcid"))
ggbarplot(residues, x = "Treatment", y = "T16_Diel_DStd",legend = "none", 
                color = "Treatment", fill = "Treatment", palette = "uchicago",
                add = "mean_se", width = 0.25, label = groups1_ordered$groups, lab.vjust = -2, xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.6, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Soil amendments") # Add horizontal line at base mean




#### GGPUBR barplot with t test #####
str(residues)
#### DT0 ####
T0D <- ggbarplot(residues, x = "Treatment", y = "T0_Diel_DStd",legend = "none", 
          color = "black", fill = "Treatment", palette = "uchicago",
          add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30","gray45"))+
  stat_compare_means(label.y = 1.3, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = mean(na.omit(residues[residues$Treatment %in% c("Control"),c("T0_Diel_DStd")])), linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("T0 Dieldrin") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(0,1.25,0.25)) 



#### DT16 ####
T16D <- ggbarplot(residues, x = "Treatment", y = "T16_Diel_DStd",legend = "none", 
          color = "black", fill = "Treatment", palette = "uchicago",
          add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue","blue","royalblue4","gray45","gray45", 
                             "gray45","gray30", "gray30","gray30", "gray30","gray30","gray45"))+
  stat_compare_means(label.y = 1.3, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T16_Diel_DStd")]), linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Dieldrin concentrations after 4 months") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                   ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(0,1.25,0.25))  

T16D


#### DDE_T16 ####
T16DDE <- ggbarplot(residues, x = "Treatment", y = "T16_DDE",legend = "none", 
                  color = "black", fill = "Treatment", palette = "uchicago",
                  add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "DDE µg / g soil")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue","blue","royalblue4","gray45","gray45", 
                             "gray45","gray30", "gray30","gray30", "gray30","gray30","gray45"))+
  stat_compare_means(label.y = 0.1, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T16_DDE")]), linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("DDE concentrations after 4 months") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)
  #scale_y_continuous(breaks=seq(0,1.25,0.25))  

T16DDE


#### Dloss  ####
D_loss <- (residues$T0_Diel_DStd - residues$T16_Diel_DStd) / residues$T0_Diel_DStd * 100
residues2 <- cbind(residues, D_loss)

Dloss <- ggbarplot(residues2[c(-1),], x = "Treatment", y = "D_loss",legend = "none", 
                  color = "black", fill = "Treatment", palette = "uchicago",
                  add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin loss (%)") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30"))+
  stat_compare_means(label.y = 20, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Dieldrin loss (%)") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(-50,50,10))  
Dloss



#### DDE loss  ####
DDE_loss <- (residues$T0_DDE - residues$T16_DDE) / residues$T0_DDE * 100
residues2 <- cbind(residues2, DDE_loss)

DDEloss <- ggbarplot(residues2[c(-1),], x = "Treatment", y = "DDE_loss",legend = "none", 
                   color = "black", fill = "Treatment", palette = "uchicago",
                   add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "DDT (%)") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30"))+
  stat_compare_means(label.y = 20, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("DTT loss (%)") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(-50,50,10))  
DDEloss



#### D + DDE loss  ####
#use residues2 from above
D_loss <- (residues$T0_Diel_DStd - residues$T16_Diel_DStd) / residues$T0_Diel_DStd * 100
residues2 <- cbind(residues, D_loss)
DDE_loss <- (residues$T0_DDE - residues$T16_DDE) / residues$T0_DDE * 100
residues2 <- cbind(residues2, DDE_loss)

D_DDEloss <- ggbarplot(residues2[c(-1),], x = "Treatment", y = c("D_loss", "DDE_loss"),legend = "none", 
                     color = "black", fill = "Treatment", palette = "uchicago", combine = TRUE,
                     add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Residue Loss (%)") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30"
                             ))+
  stat_compare_means(label.y = -20, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Loss after 4 months (%)")+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)
D_DDEloss

require(reshape2) #to get the table into the right format 
require(ggpubr)
residues2_melt <- melt(residues2, measure.vars = c(7:8), na.rm = TRUE) #measure.vars = msr)
D_DDEloss <- ggbarplot(residues2_melt, x = "Treatment", y = "value",legend = "none", 
                       color = "black", fill = "variable", palette = "uchicago",
                       add = c("mean_se"), lab.vjust = -2, xlab = FALSE,ylab = "Residue Loss (%)")+
#palette =c("black","darkgrey")) 
rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30"
  ))+
  stat_compare_means(label.y = -20, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Loss after 4 months (%)")+ # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)

D_DDEloss


require(Rmisc)
dcc1 <- summarySE(residues2[c(-1),], measurevar="residueloss", groupvars=c("Treatment")) #getting the standard error


#### ggarrange dat thang ####
ggarrange(T0D,  T16D,Dloss,legend = "none", common.legend = TRUE, nrow = 2, ncol = 2, labels = "auto")


### ggbarplot T0,T16 comparison dieldrin ####
require(reshape2) #to get the table into the right format 
residues_Donly <- residues[c(1:2,5:6)]
residues_melt <- melt(residues_Donly, measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)


ggbarplot(residues_melt, x = "variable", y = "value",legend = "none", 
          palette = "uchicago", facet.by = "Treatment",   short.panel.labs = FALSE,
          add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil",
           color = "black", fill = "variable")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  #stat_compare_means(label.y = 1.3, method = "anova")+ # Add global annova p-value
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
   # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(0,1.25,0.25)) 




### ggbarplot T0,T16 comparison ddt ####
require(reshape2) #to get the table into the right format 
residues_DDTonly <- residues[c(1:2,3:4)]
residues_melt <- melt(residues_DDTonly, measure.vars = c(3:4), na.rm = TRUE) #measure.vars = msr)


ggbarplot(residues_melt, x = "variable", y = "value",legend = "none", 
          palette = "uchicago", facet.by = "Treatment",   short.panel.labs = FALSE,
          add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil",
          color = "black", fill = "variable")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.1, method = "anova")+ # Add global annova p-value
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE) + ggtitle("DDT")




str(residues)
#DDET0
DDET0 <- ggbarplot(residues, x = "Treatment", y = "T0_DDE",legend = "none", 
                 color = "black", fill = "Treatment", palette = "uchicago",
                 add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "DDE µg / g soil")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30","gray45"))+
  stat_compare_means(label.y = 0.11, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T0_DDE")]), linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("T0 DDE") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(0,0.1,0.01)) 
DDET0

#DDET16
DDET16 <- ggbarplot(residues, x = "Treatment", y = "T16_DDE",legend = "none", 
                  color = "black", fill = "Treatment", palette = "uchicago",
                  add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "DDE µg / g soil")+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue","blue","royalblue4","gray45","gray45", 
                             "gray45","gray30", "gray30","gray30", "gray30","gray30","gray45"))+
  stat_compare_means(label.y = 0.11, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = mean(residues[residues$Treatment %in% c("Control"),c("T16_DDE")]), linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("T16 DDE") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(0,0.1,0.01))  
DDET16

#DDEloss  
residuelossDDE <- (residues$T0_DDE - residues$T16_DDE) / residues$T0_DDE * 100
residues2 <- cbind(residues2, residuelossDDE)

DDEloss <- ggbarplot(residues2[c(-1),], x = "Treatment", y = "residuelossDDE",legend = "none", 
                   color = "black", fill = "Treatment", palette = "uchicago",
                   add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "DDE loss (%)") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  scale_fill_manual(values=c("blue", 
                             "gray30", "gray30","gray30", "gray30","gray30"))+
  stat_compare_means(label.y = 25, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("DDE loss (%)") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)+
  scale_y_continuous(breaks=seq(-50,50,10))  
DDEloss

require(Rmisc)
dcc2 <- summarySE(residues2[c(-1),], measurevar="residuelossDDE", groupvars=c("Treatment")) #getting the standard error


#### ggarrange dat thang ####
ggarrange(DDET0,  DDET16,DDEloss,legend = "none", common.legend = TRUE, nrow = 2, ncol = 2, labels = "auto")







#### Ag selection ######
ANOVA=aov(Dieldrin ~ Treatment, data = residues_AG)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
target <- c("Control", "Lime", "PoultryLitter", "Biochar")
groups1_ordered <- groups1[match(target, row.names(groups1)),]


# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C)
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid")
#                       , c("CitricAcid", "FumaricAcid"), c("Inositol", "FumaricAcid"))

ggbarplot(residues_AG, x = "Treatment", y = "Dieldrin",legend = "none", 
          color = "TreatmentGroup2", fill = "TreatmentGroup2", palette = "uchicago",
          add = "mean_se", width = 0.25, lab.vjust = -2, xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.6, method = "anova")+  # Add global annova p-value
  stat_compare_means(ref.group = "Control", label = "p.signif") +
  geom_hline(yintercept = 1.218, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Soil amendments") # Add horizontal line at base mean


##### Carbon substrates #####
#Tukeytest to add p values to graph
library(rstatix)
stat.testD <- aov(Dieldrin ~ Treatment, data = residues_C) %>% tukey_hsd()
stat.testDDE <- aov(DDE ~ Treatment, data = residues_C) %>% tukey_hsd()
#ANOVA=aov(Dieldrin ~ Treatment, data = residues_C)
#myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
#comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
#groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
#target <- c("Control", "CitricAcid", "Inositol","FumaricAcid")
#groups1_ordered <- groups1[match(target, row.names(groups1)),]

# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C, ref.group = "Control")
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid")
 #                       , c("CitricAcid", "FumaricAcid"), c("Inositol", "FumaricAcid"))

#Dieldrin
p1 <-ggbarplot(residues_C, x = "Treatment", y = "Dieldrin",legend = "none", 
          color = "TreatmentGroup2", fill = "TreatmentGroup2", palette = "uchicago",
          add = "mean_se", width = 0.25, label = c("","*","*",""), lab.nb.digits = 2, lab.vjust = -2,xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means( method = "anova")+  # Add global annova p-value
  #stat_compare_means(ref.group = "Control", label = "p.format", method = "t.test") + 
  geom_hline(yintercept = 1.218, linetype = 2)+
  #stat_pvalue_manual(stat.test[c(1:3),], label = "p.adj",y.position = c(1.3, 1.35, 1.4))+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Carbon substrates")  # Add horizontal line at base mean
  #ggtext(label.df, x = "Treatment", y = "Value", label = "Value")
p1 <- ggpar(p0,yticks.by = 0.25)
p1

#DDE
mean(residues[residues$Treatment %in% c("Control"),c("DDE")]) #0.31033
p1DDE <-ggbarplot(residues_C, x = "Treatment", y = "DDE",legend = "none", 
               color = "TreatmentGroup2", fill = "TreatmentGroup2", palette = "uchicago",
               add = "mean_se", width = 0.25, label = c("","*","*",""), lab.nb.digits = 2, lab.vjust = -2,xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.35, method = "anova")+  # Add global annova p-value
  #stat_compare_means(ref.group = "Control", label = "p.format", method = "t.test") + 
  geom_hline(yintercept = 0.31033, linetype = 2)+
  #stat_pvalue_manual(stat.test[c(1:3),], label = "p.adj",y.position = c(1.3, 1.35, 1.4))+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Carbon substrates")  # Add horizontal line at base mean
#ggtext(label.df, x = "Treatment", y = "Value", label = "Value")
p1DDE
ggpar(p1DDE,yticks.by = 0.05)



##### Electron chain #####
ANOVA=aov(Dieldrin ~ Treatment, data = residues_E)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
target <- c("Control", "SodiumFormate", "SodiumAcetate")
groups1_ordered <- groups1[match(target, row.names(groups1)),]

# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C)
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid")
#                       , c("CitricAcid", "FumaricAcid"), c("Inositol", "FumaricAcid"))

p3 <- ggbarplot(residues_E, x = "Treatment", y = "Dieldrin",legend = "none", 
          color = "TreatmentGroup2", fill = "TreatmentGroup2", palette = "uchicago",
          add = "mean_se", width = 0.20, lab.vjust = -2, xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.6, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Electron Transport substrates") # Add horizontal line at base mean



##### Autoclaved and H202 #####
ANOVA=aov(Dieldrin ~ Treatment, data = residues_Other)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
target <- c("Control", "AutoclavedControl","H2O2")
groups1_ordered <- groups1[match(target, row.names(groups1)),]

# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C)
# Visualize: Specify the comparisons you want
#my_comparisons <- list( c("Control", "CitricAcid"), c("Control", "Inositol"), c("Control", "FumaricAcid")
#                       , c("CitricAcid", "FumaricAcid"), c("Inositol", "FumaricAcid"))

p4 <- ggbarplot(residues_Other, x = "Treatment", y = "Dieldrin",legend = "none", 
                color = "TreatmentGroup2", fill = "TreatmentGroup2", palette = "uchicago",
                add = "mean_se", width = 0.20, lab.vjust = -2, xlab = FALSE)+
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 1.6, method = "anova")+  # Add global annova p-value
  geom_hline(yintercept = 1.218, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("Others") # Add horizontal line at base mean

############### ggarrange barplots #########
ggarrange(p1, p2, p3,  p4, legend = "none", common.legend = TRUE, nrow = 2, ncol = 2, labels = "auto")





##### GGPLOT barplots with HSD letters #####
library(ggplot2)
library(plotrix)
library(dplyr)
library(multcompView)
library(multcomp)
library(agricolae)

require(Rmisc)
require(tidyverse)
data2 <- summarySE(residues, measurevar="Dieldrin", groupvars=c("Treatment")) #getting the standard error
data2 <- data2 %>% arrange(desc(Dieldrin))

ANOVA=aov(Dieldrin ~ Treatment, data = residues)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)

comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
groups1 <- groups1 %>% arrange(desc(Dieldrin))

ggplot(data2,aes(x=reorder(Treatment, Dieldrin, FUN = median), y=data2$Dieldrin, width=.45)) + 
  geom_bar(position="dodge",stat="identity", colour = "black") + 
  theme_classic() +
  geom_errorbar(aes(ymin=data2$Dieldrin, ymax=data2$Dieldrin + data2$se)) +
  geom_text(aes(label=groups1$groups), hjust=-2)+ coord_flip() +
  theme(axis.text.x = element_text(size=8,vjust=0.55),
        axis.title.y=element_text(size=8,angle=90),
        title = element_text(size=8, face='bold', hjust=0.5))+
  labs(x="Treatments", y = "Dielrin µg / g soil")







##############  BARPLOT with ggplot and Rmisc for error bars ########

######### Dieldrin ##################
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
require(Rmisc)
require(ggplot2)
dcc1 <- summarySE(residues, measurevar="T16_Diel_DStd", groupvars=c("Treatment")) #getting the standard error
# Error bars represent standard error of the mean
theme_set(theme_bw()) 
p1 <- ggplot(dcc1, aes(x=Treatment, y=T16_Diel_DStd, fill=Treatment)) + 
  theme(legend.position="none") +
  geom_bar(position=position_dodge(0.9), stat="identity",col ="black", width = 0.75) +
  geom_errorbar(aes(ymin=T16_Diel_DStd-se, ymax=T16_Diel_DStd+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + scale_fill_brewer(palette = "Set1")+
  labs(x="Treatment", y=expression(Dieldrin~concentration~(µg~g^-1)), title='') + 
  scale_y_continuous(breaks=seq(0,2,0.25)) +  scale_fill_manual(values=c("blue","blue","royalblue4","gray45","gray45", 
                                                                          "gray45","gray45", "gray45","gray45", "gray45","gray45","gray45"))
p1 +  geom_point(data=residues,aes(Treatment)) 
                                                                              
                                                                             

write.csv(dcc1, "meanvaluesD.csv")



###### DDE #############
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
require(Rmisc)
require(ggplot2)
dcc2 <- summarySE(residues, measurevar="DDE", groupvars=c("Treatment")) #getting the standard error
# Error bars represent standard error of the mean
theme_set(theme_bw()) 
p2 <- ggplot(dcc2, aes(x=Treatment, y=DDE))+ #, fill=Treatment)) + 
  theme(legend.position="none") +
  geom_bar(position=position_dodge(0.9), stat="identity",colour="black", width = 0.75) +
  geom_errorbar(aes(ymin=DDE-se, ymax=DDE+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + scale_fill_brewer(palette = "Set1")+
  labs(x="Treatment", y=expression(DDE~concentration~(µg~g^-1)), title='') + 
  scale_y_continuous(breaks=seq(0,2,0.25)) 
p2


###### ANOVA  #########
#Dieldrin all 
ano_resid <- lm(Dieldrin ~ Treatment, data = residues)
summary(ano_resid)
library(car)
Anova(ano_resid, type = "III")
# 1. Homogeneity of variances
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1) 
plot(ano_resid, 1)
hist(residuals(ano_resid)) #normally distributed residuals
leveneTest(Dieldrin ~ Treatment, data = residues) #levene test, to test for homogenity of variance
# 2. Normality
plot(ano_resid, 2)
shapiro.test(residuals(ano_resid))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 


#Dieldrin with residues_AG
ano_residAG <- lm(Dieldrin ~ Treatment, data = residues_AG)
summary(ano_residAG)
library(car)
Anova(ano_residC, type = "III")
# 1. Homogeneity of variances
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1) 
plot(ano_residAG, 1)
hist(residuals(ano_residAG)) #normally distributed residuals
leveneTest(Dieldrin ~ Treatment, data = residues_AG) #levene test, to test for homogenity of variance
# 2. Normality
plot(ano_residAG, 2)
shapiro.test(residuals(ano_residAG))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 


#Dieldrin with residues_C
ano_residC <- lm(Dieldrin ~ Treatment, data = residues_C)
summary(ano_residC)
library(car)
Anova(ano_residC, type = "III")
# 1. Homogeneity of variances
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1) 
plot(ano_residC, 1)
hist(residuals(ano_residC)) #normally distributed residuals
leveneTest(Dieldrin ~ Treatment, data = residues_C) #levene test, to test for homogenity of variance
# 2. Normality
plot(ano_residC, 2)
shapiro.test(residuals(ano_residC))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 


#Dieldrin with residues_E
ano_residE <- lm(Dieldrin ~ Treatment, data = residues_E)
summary(ano_residE)
library(car)
Anova(ano_residE, type = "III")
# 1. Homogeneity of variances
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1) 
plot(ano_residE, 1)
hist(residuals(ano_residE)) #normally distributed residuals
leveneTest(Dieldrin ~ Treatment, data = residues_E) #levene test, to test for homogenity of variance
# 2. Normality
plot(ano_residC, 2)
shapiro.test(residuals(ano_residE))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 



#DDE
ano_residDDE <- lm(DDE ~ Treatment, data = residues)
summary(ano_residDDE)
library(car)
Anova(ano_residDDE, type = "III")
# 1. Homogeneity of variances
plot(ano_residDDE, 1)
hist(residuals(ano_residDDE)) #normally distributed residuals
leveneTest(DDE ~ Treatment, data = residues) #levene test, to test for homogenity of variance
# 2. Normality
plot(ano_residDDE, 2)
shapiro.test(residuals(ano_residDDE))



###### Pairwise T test ##########
pairwise.t.test(residues$Dieldrin, residues$Treatment, p.adj = "none")







############ TUKEY HSD ##########
library(agricolae)
ano_resid <- lm(Dieldrin ~ Treatment, data = residues)
#Dieldrin
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 
out  <- (HSD.test(ano_resid, "Treatment"))
plot(out,variation="SE")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 


ano_residAG <- lm(Dieldrin ~ Treatment, data = residues_AG)
outAG  <- (HSD.test(ano_residAG, "Treatment"))
plot(outAG,variation="SE")

ano_residC <- lm(Dieldrin ~ Treatment, data = residues_C)
outC  <- (HSD.test(ano_residC, "Treatment",unbalanced=TRUE))
plot(outC,variation="SE")
?HSD.test
ano_residE <- lm(Dieldrin ~ Treatment, data = residues_E)
outE  <- (HSD.test(ano_residE, "Treatment"))
plot(outE,variation="SE")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) 



#DDE
outDDE  <- (HSD.test(ano_residDDE, "Treatment"))
plot(outDDE,variation="SE")

outAG  <- (HSD.test(ano_residAG, "Treatment"))
plot(outAG,variation="SE")

outC  <- (HSD.test(ano_residC, "Treatment"))
plot(outC,variation="SE")




######## Multicomp #########
require(multcomp)

tukey_spec <- glht(ano_resid, linfct=mcp(Treatment="Tukey"))
summary(tukey_spec)
tuk.cld <- cld(tukey_spec) 

par(mfrow = c(1, 1), mai=c(1,1,1.5,1))
plot(tuk.cld)
?plot
#Note that “Tukey” here does not mean Tukey-adjusted comparisons.  It just sets up a matrix to compare each mean to each other mean.
#only sodium formate is different to other treatments but not the control
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
library(car)
crPlots(ano_resid)


ano_resid <- lm(Dieldrin ~ Treatment, data = residues[c(1:3,16:19),])
summary(ano_resid)
#Just comparing control vs Inositol, p-value: 0.02836

ano_resid <- lm(Dieldrin ~ Treatment, data = residues[c(1:3,8:11),])
summary(ano_resid)
#Just comparing control vs citric acid, p-value: 0.023


######## Multicomp view, labels on boxplot #########
library(plyr)
library(ggplot2)
library(multcompView)

ano_resid <- aov(Dieldrin ~ Treatment, data = residues)
tHSD <- TukeyHSD(ano_resid, ordered = FALSE, conf.level = 0.95)

generate_label_df <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$y)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}






########## LSD test ##########
library(agricolae) 
#using the model from the ANOVA
ano_resid <- aov(Dieldrin ~ Treatment, data = residues)
summary(ano_resid)

out <- LSD.test(ano_resid, "Treatment", p.adj = "holm", console = TRUE)
out$statistics$LSD

?LSD.test
#stargraph
# Variation range: max and min
plot(out)
plot(out,variation="IQR")
plot(out,variation="SE")

# Old version LSD.test()
df<-df.residual(ano_resid)
MSerror<-deviance(ano_resid)/df
out <- with(residues,LSD.test(Dieldrin,Treatment,df,MSerror))
out



#stargraph
# Variation interquartil range: Q75 and Q25
plot(out,variation="IQR")
#endgraph
out<-LSD.test(ano_resid,"Treatment",p.adj="hommel",console=TRUE)
plot(out,variation="SD") # variation standard deviation



####### calculating the % difference of mean treatments to control ########
require(Rmisc)
dcc1 <- summarySE(residues, measurevar="Dieldrin", groupvars=c("Treatment")) #getting the standard error

#Dieldrin
#control - inositol
(dcc1[c(1),]$Dieldrin - dcc1[c(6),]$Dieldrin) / dcc1[c(1),]$Dieldrin
#-11.1%, or -0.13.5 µg / g soil 

#control - citric acid
(dcc1[c(1),]$Dieldrin - dcc1[c(5),]$Dieldrin) / dcc1[c(1),]$Dieldrin
#-10.9% or -0.13.2 µg / g


#control - sodium formate
(dcc1[c(1),]$Dieldrin - dcc1[c(7),]$Dieldrin) / dcc1[c(1),]$Dieldrin
# +2.4% or 0.03 µg / g 

#DDE
require(Rmisc)
dcc2 <- summarySE(residues, measurevar="DDE", groupvars=c("Treatment")) #getting the standard error

#control - Fumaric Acid
(dcc2[c(1),]$DDE - dcc2[c(9),]$DDE) / dcc2[c(1),]$DDE
#-4.8%, or -0.015 µg / g soil 

#control - Sodium Acetete
(dcc2[c(1),]$DDE - dcc2[c(8),]$DDE) / dcc2[c(1),]$DDE
#-3.97%, or -0.012 µg / g soil 

#control - Sodium Formate
(dcc2[c(1),]$DDE - dcc2[c(7),]$DDE) / dcc2[c(1),]$DDE
#-6.4%, or -0.02 µg / g soil 

#control - Biochar
(dcc2[c(1),]$DDE - dcc2[c(4),]$DDE) / dcc2[c(1),]$DDE
#-6.55%, or -0.02 µg / g soil 

#control - Poultry manure
(dcc2[c(1),]$DDE - dcc2[c(3),]$DDE) / dcc2[c(1),]$DDE
#-5.59%, or -0.017 µg / g soil 






######## LME Linear mixed effect model ###########
library(lme4)
library(car)
library(visreg)
library(nlme)

#load data
residues <- read.csv("Residues.csv", header = TRUE, row.names=NULL)
residues <- residues[c(-2, -23,-31),] 
str(residues)

boxplot(Dieldrin ~ Treatment, residues)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default

M2<- gls(Dieldrin ~ Treatment, data = residues)
M3<- lme(Dieldrin ~ Treatment, random =~ 1 | Cleanup, data = residues)
M4<- lme(Dieldrin ~ Treatment, random =~ 1 | Cleanup, weights = varIdent(form =~ 1 | Treatment), data = residues)
anova(M2,M3,M4)
#M3
summary(M2)


#initial diagnostics
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
E <- resid(M2, type = 'normalized')
hist(E, xlab = "Residuals", main = "")
F2 <- fitted(M2)
plot(x = F2, y = E, xlab = "Fitted values", ylab = "Residuals")
plot(residues$Treatment, E, xlab = "Treatment", ylab = "Residuals")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default
shapiro.test(E)


#remove sodium acetate and lime, residuals are much wider than rest
residues <- residues[residues$Treatment %in% c("Control", "PoultryLitter", "Biochar", 
                                               "CitricAcid", "Inositol", "SodiumFormate",
                                               "FumaricAcid", "H2O2", "AutoclavedControl"),]  

M2.2 <- gls(Dieldrin ~ Treatment, data = residues)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
E <- resid(M2.2, type = 'normalized')
hist(E, xlab = "Residuals", main = "")
F2 <- fitted(M2.2)
plot(x = F2, y = E, xlab = "Fitted values", ylab = "Residuals")
plot(residues$Treatment, E, xlab = "Treatment", ylab = "Residuals")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default

plot(M2.2)
summary(M2.2)
library(car)
Anova(M2.2, type=c("III"))







