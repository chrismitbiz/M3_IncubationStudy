setwd("~/Documents/Study/LaTrobe/Research/phD/Milestone3/R_M3")

#require(BBmisc)
#standardsR <- read.csv("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Data_GCECD/StandardsR.csv",header = TRUE, row.names=NULL)
#standardsR_centered <- normalize(standardsR, method = "scale")
#write_csv(standardsR_centered,"~/Documents/Study/LaTrobe/Research/phD/Milestone3/Data_GCECD/standardsR_centered.csv")  
#packages
# Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/rstatix")

library(agricolae)
library(ggplot2)

library(dplyr)
library(Rmisc)
library(reshape2)

##### loading all data ######
data <- read.csv("~/Documents/Study/LaTrobe/Research/phD/Milestone3/R_M3/data.csv", header = TRUE, row.names=NULL)

mod <- lm(log1p(Dieldrin) ~ H20, data = data)
summary(mod)


#renaming
residues$Treatment <- as.character(residues$Treatment)
residues$Treatment[residues$Treatment == "AutoclavedControl"] <- "Heat"
#order factors
residues$Treatment <- factor(residues$Treatment, levels=c("Control", "Lime", "PoultryLitter", "Biochar", 
                                                                "CitricAcid", "Inositol", "SodiumFormate", "SodiumAcetate",
                                                                "FumaricAcid", "Heat","H2O2"))

mapping <- c("T0" = 0, "T16" = 16, "T39" = 39,"T52" = 52)
residues$Week_cont <- mapping[residues$Week]


residues_na <- na.omit(residues)
str(residues_na)
residues_na$Week


# some initial look
plot(Dieldrin ~ Week_cont, data = residues_na)
plot(Dieldrin ~ Week, data = residues)

boxplot(Dieldrin ~ Week_cont, data = residues_na)
stripchart(Dieldrin ~ Week_cont, data = residues_na,
           method = "jitter", vertical = TRUE, add = TRUE,pch = 20)
boxplot(residues_na$Dieldrin)
stripchart(residues_na$Dieldrin,
           method = "jitter", vertical = TRUE, add = TRUE,pch = 20)

boxplot(Dieldrin ~ Treatment, data = residues_na)
stripchart(Dieldrin ~ Treatment, data = residues_na,
           method = "jitter", vertical = TRUE, add = TRUE,pch = 20)

shapiro.test(residues_na$Dieldrin)
hist(residues_na$Dieldrin)
densityplot(residues_na$Dieldrin)


shapiro.test(log(residues_na$Dieldrin+1))
hist(log(residues_na$Dieldrin+1))


plot(Dieldrin ~ H20, data = residues_na)

boxplot(H20 ~ Week_cont, data = residues_na)
stripchart(H20 ~ Week_cont, data = residues_na,
           method = "jitter", vertical = TRUE, add = TRUE,pch = 20)
shapiro.test(residues_na$H20)
histogram(residues_na$H20)



residues %>% 
ggplot(data = ., aes(x = Week, y = Dieldrin)) +
  geom_boxplot( outlier.shape  = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "", y = "Dieldrin\n") +
  #facet_wrap(~ OTU, scales = "free") +
  ggtitle("All treatments")




###### GGPUBR box plot ##### 

#mean of control, dieldrin

#Dieldrin overview across all weeks
#subsettign by week if wanted
#subset(residues_na, Week %in% c("T0"))
residues %>%
ggboxplot(data = ., x = "Week", y = "Dieldrin",  add = "jitter",
          color = "Category", palette = "simpsons",xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(aes(group = Week),label = "..p.signif..", paired = TRUE) +
  #stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  #stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
  #                  ref.group = "Control", hide.ns = TRUE)+
  ggtitle("")
boxpplot 


residues %>% select(Samples:Dieldrin,Week_cont) %>%
  ggboxplot(data = ., x = "Week_cont", y = "Dieldrin",  add = "jitter",
            color = "Category", palette = "simpsons",xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  stat_compare_means(aes(group = Week),label = "..p.signif..", paired = TRUE)



#boxp1D <- boxp1D + scale_y_continuous(breaks=seq(0,1.35,0.05))
boxpplot <- ggboxplot(residues_na, x = "Week", y = "Dieldrin",  add = "jitter",
                      color = "Week",xlab = FALSE, label = "Treatment"
                      ) +
  rotate_x_text(angle = 45)+ 
  #stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  #stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
  #                  ref.group = "Control", hide.ns = TRUE)+
  ggtitle("")

# Perorm pairwise comparisons
#compare_means(Dieldrin ~ Treatment,  data = residues_C, method = "t.test")
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






##### GGPUBR line plots ####
#subset(residues_na, Treatment %in% c("Control", "Lime"))

#to get continuous x axis (numeric value)
residues_na$Treatment
str(residues_na)
lineplot1 <- ggline(subset(residues_na, Treatment %in% c("Control", "AutoclavedControl", "Lime", "Biochar", "PoultryLitter")), x = "Week", y = "Dieldrin",  
                   add = c("mean_se", "dotplot"),
                    color = "Treatment", palette = "uchicago",xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  #stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  #stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
  #                  ref.group = "Control", hide.ns = TRUE)+
  ggtitle("")
lineplot1 

lineplot2 <- ggline(subset(residues_na, Treatment %in% c("Control", "CitricAcid", "Inositol","AutoclavedControl"
                                                         )), x = "Week", y = "Dieldrin",  
                    add = c("mean_se", "dotplot"),
                    color = "Treatment", palette = "uchicago",xlab = FALSE) +
  rotate_x_text(angle = 45)+ 
  #stat_compare_means(label.y = 1.4, method = "anova")+  # Add global annova p-value
  #stat_compare_means(label = "p.format", method = "t.test",   # Pairwise comparison against all
  #                  ref.group = "Control", hide.ns = TRUE)+
  ggtitle("")
lineplot2





##### GGPUBR scatter plots ####
scatterplot <- ggscatter(residues_na , x = c("Week_cont"), y = "Dieldrin",
          add = "loess", facet.by = "Treatment", palette = "simpsons",
          col = "colours", repel = TRUE, ylab = "Dieldrin µg / g",xlab = "Week",
          add.params = list(color = "black", size = 0.5), # Customize reg. line
          conf.int = FALSE, legend = "none")




ggscatter(residues , x = c("Week"), y = "Dieldrin",
                         add = "loess", facet.by = "Treatment", palette = "simpsons",
                          repel = TRUE, ylab = "Dieldrin µg / g",xlab = "Week",
                         add.params = list(color = "black", size = 0.5), # Customize reg. line
                         conf.int = FALSE, legend = "none", col = "Category") 

ggscatter(data , x = c("H20"), y = "Dieldrin",
                         add = "loess",  palette = "simpsons",
                         repel = TRUE, ylab = "Dieldrin µg / g", xlab = "Watercontent",
                         add.params = list(color = "black", size = 0.5), # Customize reg. line
                         conf.int = FALSE, legend = "top", color = "Week") 




#scatter with average dotted line
Davg <- residues %>% 
  select(Dieldrin, Week) %>% 
  filter(Week == "T0") %>% 
  summarize(mean(Dieldrin))

DavgT52 <- residues %>% 
  select(Dieldrin, Week) %>% 
  filter(Week == "T52")%>%  filter(Dieldrin > 0) %>% 
  summarize(mean(Dieldrin))

residues %>% select(Samples:Dieldrin,Week_cont) %>%
ggscatter(data = . , x = c("Week_cont"), y = "Dieldrin",
          add = "loess", facet.by = "Treatment", palette = "simpsons", repel = TRUE, ylab = "Dieldrin µg / g",xlab = "Week",
          add.params = list(color = "black", size = 0.5), # Customize reg. line
          conf.int = FALSE, legend = "none", col = "Category") +
  geom_hline(yintercept = Davg$`mean(Dieldrin)`, linetype = "dotted")+
geom_hline(yintercept = DavgT52$`mean(Dieldrin)`, linetype = "dotted")


ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss_scatter_averagedottedlines.pdf", 
       height=5, width=6, units='in', dpi=600)







###### GGPUBR bar plot  ##### 
#Dieldrin T0
        #add labels
        avg <- summarySE(residues_na %>% filter(Week == "T0"), measurevar="Dieldrin", groupvars=c("Treatment"))
      ANOVA=aov(log(Dieldrin+1) ~ Treatment, data = residues_na %>% filter(Week == "T0"))
      myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
      comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
      #comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
      groups1 <- as.vector(comparison$groups)
      #order groups if needed
      #groups1 <- groups1 %>% arrange(desc(Dieldrin))
      target <- c("Control",  "Lime", "PoultryLitter", 
                  "Biochar", "CitricAcid", "Inositol", "SodiumFormate",
                  "SodiumAcetate", "FumaricAcid", "Heat","H2O2")
      groups1_ordered <- groups1[match(target, row.names(groups1)),]
      label <- as.vector(groups1_ordered$groups)
      
barplotT0 <- ggbarplot(residues_na %>% filter(Week == "T0"), x = "Treatment", y = "Dieldrin",
                           color = "black", fill = "colours", palette = "simpsons", label = FALSE,
                           add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil") +
        rotate_x_text(angle = 45) +
        #  stat_compare_means(aes(group = Treatment),label = "..p.signif..", paired = TRUE, method = "anova") +  # Pairwise comparison against all
        stat_compare_means(label.y = 1.5, method = "anova")+ # Add global annova p-value
      #geom_hline(yintercept = mean(residues_na$Dieldrin), linetype = 2) # Add horizontal line at base mean
annotate("text", x = c(avg$Treatment),
         y=c(1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3), 
         label = label)
      barplotT0


      
      
      #Dieldrin T39
      #add labels
      avg <- summarySE(residues_na %>% filter(Week == "T39"), measurevar="Dieldrin", groupvars=c("Treatment"))
      ANOVA=aov(log(Dieldrin+1) ~ Treatment, data = residues_na %>% filter(Week == "T39"))
      myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
      comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
      #comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
      groups1 <- as.vector(comparison$groups)
      #order groups if needed
      #groups1 <- groups1 %>% arrange(desc(Dieldrin))
      target <- c("Control",  "Lime", "PoultryLitter", 
                  "Biochar", "CitricAcid", "Inositol", "SodiumFormate",
                  "SodiumAcetate", "FumaricAcid", "Heat","H2O2")
      groups1_ordered <- groups1[match(target, row.names(groups1)),]
      label <- as.vector(groups1_ordered$groups)
      
      barplotT39 <- ggbarplot(residues_na %>% filter(Week == "T39"), x = "Treatment", y = "Dieldrin",
                             color = "black", fill = "colours", palette = "simpsons", label = FALSE,
                             add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil") +
        rotate_x_text(angle = 45) +
        #  stat_compare_means(aes(group = Treatment),label = "..p.signif..", paired = TRUE, method = "anova") +  # Pairwise comparison against all
        stat_compare_means(label.y = 1.5, method = "anova")+ # Add global annova p-value
        #geom_hline(yintercept = mean(residues_na$Dieldrin), linetype = 2) # Add horizontal line at base mean
        annotate("text", x = c(avg$Treatment),
                 y=c(1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3), 
                 label = label)
      barplotT39

      
      
      
      
barplotT0T39 <- ggbarplot(residues[-c(45:88),] , x = "Treatment", y = "Dieldrin",
                              color = "Week", fill = "Week", palette = c("cornflowerblue", "gray45"), label = FALSE,
                              add = c("mean_se"),error.plot = "upper_errorbar", lab.vjust = -2, xlab = FALSE,ylab = "Dieldrin µg / g soil",
                          position = position_dodge(0.9)) +
        rotate_x_text(angle = 45)     
            
barplotT0T39    


# Bar plot of mean +/-se
#Colour palettes
#ggsci scientific journal palettes, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".

#  stat_compare_means() +      # Global p-value
#  stat_compare_means(ref.group = "0.5", label = "p.signif", label.y = c(22, 29))  # compare to ref.group


#### Dloss  ####
#residueT0 <- residues_na[c(-23,-25,-28, -39,-42),] %>% filter(Week == "T0")
#residueT16 <- residues_na %>% filter(Week == "T16")
#residueT39 <- residues_na[c(-102,-104,-107,-110,-111,-120,-123),] %>% filter(Week == "T39")
#lossT16 <- ((residueT0$Dieldrin - residueT16$Dieldrin) / residueT0$Dieldrin)
#lossT39 <- ((residueT0$Dieldrin - residueT39$Dieldrin) / residueT0$Dieldrin)
#totalloss <- lossT16+lossT39

#for comparing T0 an T39 no need to remove rows. only when including T16 as sample 2 is missing. 
#residueT0 <- residues_na[c(-2),] %>% filter(Week == "T0")
residueT0 <- residues %>% filter(Week == "T0")
#residueT16 <- residues_na %>% filter(Week == "T16")
#residueT39 <- residues_na[c(-89),] %>% filter(Week == "T39")
residueT39 <- residues %>% filter(Week == "T39")
residueT52 <- residues %>% filter(Week == "T52")
#lossT16 <- ((residueT0$Dieldrin - residueT16$Dieldrin) / residueT0$Dieldrin)
lossT39 <- ((residueT0$Dieldrin - residueT39$Dieldrin) / residueT0$Dieldrin * 100)
lossT52 <- ((residueT0$Dieldrin - residueT52$Dieldrin) / residueT0$Dieldrin * 100)
lossT39T52 <- ((residueT39$Dieldrin - residueT52$Dieldrin) / residueT39$Dieldrin * 100)
#lossdata <- cbind(lossT16, lossT39, residueT0)
lossdataT39 <- cbind(lossT39, residueT0)
lossdataT52 <- cbind(lossT52, residueT0)
lossdataT52 <- lossdataT52 %>% filter(lossT52 >0)
lossdataT39T52 <- cbind(lossT39T52, residueT0)
lossdataT39T52 <- lossdataT39T52 %>% filter(lossT52 >0)


lossdata_mlt <- melt(lossdata, measure.vars = c(1:2))

#loss from T0 to T39

#boxplot
#loss from T0 to T39
data.frame(Samples = residues$Samples, dloss = (residueT0$Dieldrin - residueT39$Dieldrin) / residueT0$Dieldrin * 100) %>% 
  left_join(residueT0) %>% 
  ggboxplot(data = ., x = "Treatment", y = "dloss",  add = "jitter",
            color = "Category", palette = "simpsons",xlab = FALSE, title = "loss from T0 to T39") +
  rotate_x_text(angle = 45) 


#boxplot
#loss from T0 to T52
data.frame(Samples = residues$Samples, dloss = (residueT0$Dieldrin - residueT52$Dieldrin) / residueT0$Dieldrin * 100) %>% 
  left_join(residueT0) %>% 
  ggboxplot(data = ., x = "Treatment", y = "dloss",  add = "jitter",
            color = "Category", palette = "simpsons",xlab = FALSE,title = "loss from T0 to T52") +
  rotate_x_text(angle = 45)


#boxplot
#loss from T39 to T52
data.frame(Samples = residues$Samples, dloss = (residueT39$Dieldrin - residueT52$Dieldrin) / residueT39$Dieldrin * 100) %>% 
  left_join(residueT0) %>% 
  ggboxplot(data = ., x = "Treatment", y = "dloss",  add = "jitter",
            color = "Category", palette = "simpsons",xlab = FALSE,title = "loss from T39 to T52") +
  rotate_x_text(angle = 45)


#loss barplot T0-T52
ggbarplot(lossdataT52, x = "Treatment", y = c("lossT52"),legend = "none", 
          color = "black", fill = "Category",palette = "simpsons",
          add = c("mean_se"), lab.vjust = -2, xlab = FALSE,
          ylab = "") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 25,label.x = 2, method = "kruskal")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Dieldrin loss until T52") +  # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                       ref.group = "Control", hide.ns = TRUE)
ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss_total.pdf", 
       height=3.5, width=4, units='in', dpi=600)


#loss barplot T0-T39
ggbarplot(lossdataT39, x = "Treatment", y = c("lossT39"),legend = "none", 
          color = "black", fill = "Category",palette = "simpsons",
          add = c("mean_se"), lab.vjust = -2, xlab = FALSE,
          ylab = "") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 25, label.x = 2,method = "kruskal")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Dieldrin loss until T39 ") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)
ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss_T0_T39.pdf", 
       height=3.5, width=4, units='in', dpi=600)



#loss barplot T39 - T52
ggbarplot(lossdataT39T52, x = "Treatment", y = c("lossT39T52"),legend = "none", 
          color = "black", fill = "Category",palette = "simpsons",
          add = c("mean_se"), lab.vjust = -2, xlab = FALSE,
          ylab = "") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 25,label.x = 2, method = "kruskal")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Dieldrin loss from T39 to T52")+  # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                      ref.group = "Control", hide.ns = TRUE)

ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss_T39_T52.pdf", 
       height=3.5, width=4, units='in', dpi=600)









#add labels
avg <- summarySE(lossdata, measurevar="lossT39", groupvars=c("Treatment"))
ANOVA=aov(log(lossT39+1) ~ Treatment, data = lossdata)
myTukey=TukeyHSD(x=ANOVA, conf.level = 0.95, return=TRUE)
comparison = HSD.test(ANOVA, "Treatment", group=TRUE)
#comparison = LSD.test(ANOVA, "Treatment", group=TRUE)
groups1 <- as.vector(comparison$groups)
#order groups if needed
#groups1 <- groups1 %>% arrange(desc(Dieldrin))
target <- c("Control",  "Lime", "PoultryLitter", 
            "Biochar", "CitricAcid", "Inositol", "SodiumFormate",
            "SodiumAcetate", "FumaricAcid", "Heat","H2O2")
groups1_ordered <- groups1[match(target, row.names(groups1)),]
label <- as.vector(groups1_ordered$groups)
#need to figure out how to display the label vector in correct order

ggbarplot(lossdata, x = "Treatment", y = c("lossT39"),legend = "none", 
                    color = "black", fill = "colours",palette = "simpsons",
                    add = c("mean_se", "dotplot"), lab.vjust = -2, xlab = FALSE,
                    ylab = "") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 20, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Dieldrin loss (%)") + # Add horizontal line at base mean
#  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
#                     ref.group = "Control", hide.ns = FALSE)
  annotate("text", x = c(avg$Treatment),
           y=c(17,17,17,17,17,17,17,17,17,17,17), 
           label = label)





#using melted data to get a stacked T16, T39 loss
Dloss2 <- ggbarplot(lossdata_mlt, x = "Treatment", y = "value", 
                   color = "variable", fill = "variable",palette = "simpsons",
                   add = c("mean", "dotplot"), lab.vjust = -2, 
                   xlab = FALSE,ylab = "Dieldrin loss (%)") +
  #palette =c("black","darkgrey")) 
  rotate_x_text(angle = 45)+ 
  stat_compare_means(label.y = 0.2, method = "anova")+ # Add global annova p-value
  geom_hline(yintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0.104, linetype = 2)+
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  ggtitle("") + # Add horizontal line at base mean
  stat_compare_means(label = "..p.signif..", method = "t.test",   # Pairwise comparison against all
                     ref.group = "Control", hide.ns = TRUE)
# Change the default order of items
#order = c("D2", "D1", "D0.5"))




Dloss3 <- ggboxplot(lossdata_mlt, x = "variable", y = "value",  add = "jitter",
                      color = "Treatment", palette = "simpsons",xlab = FALSE,
                    facet.by = "Treatment", legend = "none") +
  rotate_x_text(angle = 45)+ ylab("Dieldrin loss (%)") +
  ggtitle("Dieldrin loss (%)")

Dloss3 



### ggarrange #### 
ggarrange( Dloss1, scatterplot, barplotT0T39,ncol = 2, nrow = 2,common.legend = FALSE,
         heights = c(1, 1))
ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss_overview.pdf", 
       height=8, width=10, units='in', dpi=600)

ggarrange(scatterplot,Dloss1,common.legend = FALSE,
          heights = c(1.4, 1))

ggsave("~/Documents/Study/LaTrobe/Research/phD/Milestone3/Figures/loss.pdf", 
       height=4.5, width=10, units='in', dpi=600)


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



######## Mixed effect model - repeated measures ###########
library(lme4)
library(car)
library(visreg)
library(nlme)
library(pbkrtest)
library(lmerTest)
library(sjPlot)

boxplot(Dieldrin ~ Treatment, residues)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) #default

#from ch 6 -  Mixed effect models and extension in ecology with R (zuur et al)
M0 <- gls(Dieldrin ~ Week_cont+Treatment,
          na.action = na.omit, data = residues_na)
summary(M0)
#visually inspect independance of residuals (as it is the same subject it is not independant)
E <- residuals(M0, type = "normalized")
#I1 <- !is.na(residues$Dieldrin)
#Efull <- vector(length = length(residues$Dieldrin))
#Efull <- NA
#Efull[I1] <- E
acf(E, na.action = na.pass,
      main = "Auto-correlation plot for residuals")
plot(residues_na$Treatment, E, xlab = "Treatment", ylab = "Residuals")
plot(residues_na$Week, E, xlab = "Treatment", ylab = "Residuals")
plot(residues_na$Dieldrin, E, xlab = "concentrations", ylab = "Residuals", type = "n")
text(residues_na$Dieldrin, E, labels = residues_na$Samples)
hist(E, xlab = "Residuals", main = "")
shapiro.test(E)



#now compare AIC,BIC of model with and without autocorrelation structure
M1 <- gls(Dieldrin ~ Week_cont + Treatment,
          na.action = na.omit, data = residues_na,
          correlation = corCompSymm(form =~ Samples))
anova(M0,M1) #no improvement in model
plot(M0)
qqnorm(resid(M0))
#plot_grid(plot_model(M0, type = "diag")) #works for lme but not gls

# Explore the fixed structure
M0 <- lme(Dieldrin ~ Week+Treatment, random =~ 1 | Cleanup,
          data = residues_na, method = "ML")
M1<- lme(Dieldrin ~ Week*Treatment, random =~ 1 | Cleanup, 
         data = residues_na, method = "ML")
anova(M0,M1)

# Explore the random structure
M0 <- gls(Dieldrin ~ Treatment, data = residues_na)
M1<- lme(Dieldrin ~ Treatment, random =~ Week | Samples, data = residues_na)
M2<- lme(Dieldrin ~ Treatment, random =~ 1 | Samples, data = residues_na)
M3<- lme(Dieldrin ~ Treatment, random =~ 1 | Week/Treatment, data = residues_na)
M4<- lme(Dieldrin ~ Treatment, random =~ 1 | Week, data = residues_na)
M5<- lme(Dieldrin ~ Treatment, random =~ 1 | Week/Samples, data = residues_na)
M6<- lme(Dieldrin ~ Treatment, random =~ 1 | Cleanup, data = residues_na)
M7<- lme(Dieldrin ~ Treatment, random =~ 1 | Grinding, data = residues_na)
M8<- lme(Dieldrin ~ Treatment, random =~ Week + Treatment | Grinding, data = residues_na)
anova(M0,M1,M2,M3,M4,M5,M6,M7)

summary(M6)
plot(M6)
plot(fitted(M4), residuals(M4))
qqPlot(residuals(M4))
Anova(M4, test="F")
plot_grid(plot_model(M6, type = "diag")) #works for lme but not gls
#check interactions
with(residues_na, interaction.plot(Treatment, Cleanup, Dieldrin))
with(residues_na, interaction.plot(Treatment, Week, Dieldrin))

#follow more on the tutorial below to add correlation structure etc

plot(M2.2)
summary(M2.2)
library(car)
Anova(M2.2, type=c("III"))



#http://www.flutterbys.com.au/stats/tut/tut9.3a.html#mjx-eqn-tabrandomizedBlockAnova-correlationStructures
#Section 8
#follow tutprial
ggplot(residues_na) + geom_boxplot(aes(y = Dieldrin, x = as.factor(Week)))+
  geom_jitter(aes(y = Dieldrin, x = as.factor(Week), colour = Treatment))








### linear regression to get half life #####
#https://chemicalstatistician.wordpress.com/2013/03/24/inference-in-simple-linear-regression-with-logarithmic-transformation-an-illustration-with-exponential-decay-of-ddt-in-trout/


plot(log(Dieldrin) ~ Week_cont, data = residues_na)
reg <- lm (log(Dieldrin) ~ Week_cont, data = residues_na)
summary(reg)
decayrate <- reg$coefficients["Week_cont"]
#halflife 
halflife <- log(0.5)/decayrate
#368.9 weeks or 7 years
halflife/52
# across all 11 treatments the halflife is 7 years


#Halflife control
reg <- lm(log(Dieldrin) ~ Week_cont, data = subset(residues_na, Treatment %in% c("Control")))
plot(log(Dieldrin) ~ Week_cont, data = subset(residues_na, Treatment %in% c("Control")))
summary(reg)
decayrate <- reg$coefficients["Week_cont"]
#halflife 
halflife <- log(0.5)/decayrate
halflife/52
#4.8 years

#Halflife Acetate
reg <- lm(log(Dieldrin) ~ Week_cont, data = subset(residues_na, Treatment %in% c("SodiumAcetate")))
plot(log(Dieldrin) ~ Week_cont, data = subset(residues_na, Treatment %in% c("SodiumAcetate")))
summary(reg)
decayrate <- reg$coefficients["Week_cont"]
#halflife 
halflife <- log(0.5)/decayrate
halflife/52
#3.4 years