# set user specific working directory
setwd('./strains_github/results/ExtendedData_Figure1/')
getwd()

################################################################################
# lme test results; panel 1b
################################################################################
library(tidyverse)
library("lmerTest")
library("lme4")
library("sjPlot")
library("lattice")
library("multcomp")
library("nlme")
library("glmmTMB")

df1 <-
  read_tsv("./data_extended_figure1c_weights.tsv") %>% rename("StartDate" = 6) %>%
  pivot_longer(-(1:6), names_to = "Week", values_to = "Mass") %>%
  filter(
    Genotype != "NCoR_KO",
    Genotype != "SPRET/eiJ",
    Group != "SPRET/eiJ",
    Week != "32",
    Week != "34"
  )

df1 <- 
  df1 %>% mutate(
    Genotype = factor(Genotype),
    Diet = factor(Diet),
    Group = factor(Group),
    Week = as.numeric(Week)
  )

logdf1 <- 
  df1 %>% mutate(
    Genotype = factor(Genotype),
    Diet = factor(Diet),
    Group = factor(Group),
    Week = as.numeric(Week),
    Mass = log(Mass)
  )

str(df1)

# model used in manuscript
fm8 <- lmer(Mass ~ 1 + Week * Diet * Genotype + (1 + Week | Mouse), data = df1, REML = F)
summary(fm8)
aov8 <- anova(fm8)
aov8 %>% knitr::kable(format = "simple")

plot_model(fm8, type = "diag")
plot_model(fm8, type = "est")
plot_model(fm8, type = "pred")
plot_model(fm8, type = "int")

aj <- filter(df1, Genotype == "A/J")
aj_fm1 <- lmer(Mass ~ 1 + Week + (1 + Week | Mouse), data = aj, REML = F)
aj_fm2 <- lmer(Mass ~ 1 + Week * Diet + (1 + Week | Mouse), data = aj, REML = F)
anova(aj_fm1, aj_fm2)
anova(aj_fm1)
anova(aj_fm2)
summary(aj_fm2)

balb <- filter(df1, Genotype == "BALB/cJ")
balb_fm1 <- lmer(Mass ~ 1 + Week + (1 + Week | Mouse), data = balb, REML = F)
balb_fm2 <- lmer(Mass ~ 1 + Week * Diet + (1 + Week | Mouse), data = balb, REML = F)
anova(balb_fm1, balb_fm2)
summary(balb_fm2)

c57 <- filter(df1, Genotype == "C57BL/6J")
c57_fm1 <- lmer(Mass ~ 1 + Week + (1 + Week | Mouse), data = c57, REML = F)
c57_fm2 <- lmer(Mass ~ 1 + Week * Diet + (1 + Week | Mouse), data = c57, REML = F)
anova(c57_fm1, c57_fm2)
anova(c57_fm2)
summary(c57_fm2)

##################### comparisions of unused models #####################
fm1 <- lmer(Mass ~ 1 + (1 + Week | Mouse), data = df1, REML = F)
fm2 <- lmer(Mass ~ 1 + Week + (1 + Week | Mouse), data = df1, REML = F)
fm3 <- lmer(Mass ~ 1 + Diet + (1 + Week | Mouse), data = df1, REML = F)
fm4 <- lmer(Mass ~ 1 + Genotype + Diet + (1 + Week | Mouse), data = df1, REML = F)
fm5 <- lmer(Mass ~ 1 + Genotype + Diet + Week + (1 + Week | Mouse), data = df1, REML = F)
fm6 <- lmer(Mass ~ 1 + Genotype + Diet * Week + (1 + Week | Mouse), data = df1, REML = F)
fm7 <- lmer(Mass ~ 1 + Genotype * Week + Diet + (1 + Week | Mouse), data = df1, REML = F)
fm8 <- lmer(Mass ~ 1 + Week * Diet * Genotype + (1 + Week | Mouse), data = df1, REML = F)
anova(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8) %>% knitr::kable()

anova(fm5, fm6, fm7, fm8)

summary(fm7)
aov7 <- anova(fm7)
aov7 %>% knitr::kable(format = "simple")


################################################################################
# panel 1c; weight change plot
################################################################################
library(RColorBrewer)
library(cowplot)
library(patchwork)

data <- read_tsv("./data_extended_figure1c_weights.tsv") %>%
  pivot_longer(cols = -(1:6))

colnames(data) <- c("Genotype", "Mouse", "Diet", "Group", "DOB", "StartDate", "Week", "Mass")

notSpret <- filter(data, Genotype != "NCoR_KO", Genotype != "SPRET/eiJ", Week != "32", Week != "34")
notSpret$Week = as.numeric(as.character(notSpret$Week))

notSpret %>% filter(Mass != "NA") %>%
  group_by(Genotype, Diet, Week) %>% 
  summarize(n = n()) %>%
  pivot_wider(names_from = c("Genotype", "Diet"), values_from = c("n")) %>%
  write_csv("./extended_figure1bc_sample_size.csv")

functionBoxplot <- function(df){
  df %>%
    ggplot(aes(
      x = factor(Week),
      y = Mass,
      fill = Group
    )) +
    geom_boxplot(outlier.shape = NA, lwd = 0.25) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.1, alpha = 0.5, shape = "bullet") +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      line = element_line(size = 0.5, colour = "black"),
      rect = element_rect(size = 0.5, colour = "black"),
      axis.text.x  = element_text(angle = 0, vjust = 0.75),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 6, colour = "black"),
      axis.ticks = element_line(size = 0.5, colour = "black"),
      legend.position = "none", legend.key.width = unit(0.2, "cm"),
      strip.background = element_blank(), strip.text = element_blank()) + 
    labs(x = "", y = "Grams") 
}

functionBoxplot(notSpret) + facet_wrap(~ Genotype, nrow = 3) +
  scale_fill_manual(values = c("#a50f15", "#fcbba1", "#08519c","#c6dbef", "#006d2c","#c7e9c0"))
ggsave(filename = "extended_figure1c.pdf", height = 3.65, width = 2)

################################################################################
# scoring; panel 1e
################################################################################

#Steatosis score
a = c(1,0,1,1,2)
b = c(0,0,0,0,0)
c = c(3,3,3,3,3)
data1 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y1 = kruskal.test(data1)

#Lobular inflammation score
a = c(2,0,1,1,1)
b = c(0,0,0,0,0)
c = c(1,1,1,1,1)
data2 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y2 = kruskal.test(data2)

#NASH CRN score
a = c(3,0,2,2,3)
b = c(0,0,0,0,0)
c = c(5,5,4,4,4)
data3 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y3 = kruskal.test(data3)

#Fibrosis score
a = c(0,0,0,0,0)
b = c(0,0,0,0,0)
c = c(3,2,1,3,1)
data4 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y4 = kruskal.test(data4)

#Steatosis score
a = c(1,1,1,1,2)
b = c(0,0,0,0,0)
c = c(3,3,3,3,3)
data5 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y5 = kruskal.test(data5)

#Lobular inflammation score
a = c(1,3,1,2,1)
b = c(0,0,0,0,0)
c = c(1,1,2,1,1)
data6 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y6 = kruskal.test(data6)

#NASH CRN score
a = c(2,4,2,3,3)
b = c(0,0,0,0,0)
c = c(4,4,5,5,4)
data7 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y7 = kruskal.test(data7)

#Fibrosis score
a = c(0,1,1,0,0)
b = c(0,0,0,0,0)
c = c(2,0,2,2,3)
data8 = list('A/J'=a, 'BALB/cJ'=b, 'C57BL/6'=c)
y8 = kruskal.test(data8)

pdf("strainsAMLNScoring.pdf", width = 9.5, height = 4.5, pointsize = 8)
par(mfrow = c(2,4), mar=c(2, 2, 2, 1) + 0.1)#it goes c(bottom, left, top, right) 
stripchart(data1, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Steatosis 20 Wk\np = ",y1$p.value), pch = 21, col = 'black', bg='grey', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data2, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Inflammatory Activity 20 Wk\np = ",y2$p.value), pch = 21, col = 'black', bg='grey', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data3, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("NASH CRN 20 Wk\np = ",y3$p.value), pch = 21, col = 'black', bg='grey', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data4, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Fibrosis 20 Wk\np = ",y4$p.value), pch = 21, col = 'black', bg='grey', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data5, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Steatosis 30 Wk\np = ",y5$p.value), pch = 21, col = 'black', bg='blue', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data6, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Inflammatory Activity 30 Wk\np = ",y6$p.value), pch = 21, col = 'black', bg='blue', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data7, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("NASH CRN 30 Wk\np = ",y7$p.value), pch = 21, col = 'black', bg='blue', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
stripchart(data8, method = "jitter", vertical = T, las = T, ylab = "Score", 
           main = paste("Fibrosis 30 Wk\np = ",y8$p.value), pch = 21, col = 'black', bg='blue', 
           cex = 2, cex.main = 1, cex.axis = 1, cex.axis = 0.7)
dev.off()

