# --------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-05 20:42:50
# @Last Modified by:
# @Last Modified time: 2023-03-02 20:52:23
# @Description: DNA methylation for 3 Gens
# --------------------------------------------------------------------------#
# Clear memory
rm(list = ls())
gc()

# Library packages as needed
library(broom)
library(car)
library(data.table)
library(ggplot2)
library(ggpubr)
library(MASS)
library(phia) # for interaction
library(purrr)
library(readxl)
library(dplyr)

# Set working directory
setwd("G:/我的坚果云/002项目/2022_Chenlin/投稿期刊plants/data summary")
getwd()

chen_data_01gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "1st_gen_20210305")
chen_data_02gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "2nd_gen_20210305")
chen_data_03gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "3rd_gen_20210305")

# Modify data format
chen_data_01gen <- chen_data_01gen %>%
      dplyr::mutate(X1gen_treatment = as.factor(X1gen_treatment))

chen_data_02gen <- chen_data_02gen %>%
      dplyr::mutate(
            X1gen_treatment = as.factor(X1gen_treatment),
            X1gen_organ_origin = as.factor(X1gen_organ_origin)
      )

chen_data_03gen <- chen_data_03gen %>%
      dplyr::mutate(
            X1gen_treatment = as.factor(X1gen_treatment),
            X3gen_treatment = as.factor(X3gen_treatment),
            X1gen_organ_origin = as.factor(X1gen_organ_origin)
      )

str(chen_data_01gen)
str(chen_data_02gen)
str(chen_data_03gen)

# Check the homogeneity of variance of data
leveneTest.01gen <- leveneTest(X5mc ~ X1gen_treatment, data = chen_data_01gen)
leveneTest.01gen

leveneTest.02gen <- leveneTest(X5mc ~ X1gen_treatment * X1gen_organ_origin, data = chen_data_02gen)
leveneTest.02gen

leveneTest.03gen <- leveneTest(X5mc ~ X1gen_treatment * X1gen_organ_origin, data = chen_data_03gen)
leveneTest.03gen

# t.test for G1
t.test.01gen <- t.test(X5mc ~ X1gen_treatment, var.equal = TRUE, conf.level = 0.95, data = chen_data_01gen)
t.test.01gen

# two-way ANOVA for G2
type3 <- list(X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)
lm_5mc_2gen <- lm(X5mc ~ X1gen_treatment * X1gen_organ_origin, data = chen_data_02gen, contrasts = type3)
x2aov_5mc <- Anova(lm_5mc_2gen, type = "III")
x2aov_5mc

# two-way ANOVA for G3
type3 <- list(X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)
lm_5mc_3gen <- lm(X5mc ~ X1gen_treatment * X1gen_organ_origin, data = chen_data_03gen, contrasts = type3)
x3aov_5mc <- Anova(lm_5mc_3gen, type = "III")
x3aov_5mc

# G1
# Group the data by X1gen_organ_origin and summarize the variables
est_5mc_1gen <- chen_data_01gen %>%
      filter(!is.na(X5mc)) %>%
      dplyr::group_by(X1gen_treatment) %>%
      dplyr::summarise(
            N    = length(X5mc),
            mean = mean(X5mc, na.rm = TRUE),
            sd   = sd(X5mc, na.rm = TRUE),
            se   = sd(X5mc, na.rm = TRUE) / sqrt(N)
      )

# Generate the plot for G1
plot_5mc_1gen <- ggplot(
      est_5mc_1gen,
      aes(x = as.factor(X1gen_treatment), y = mean, group = 1, color = "blue")
) +
      geom_errorbar(width = .1, aes(ymin = mean - se, ymax = mean + se), size = 1) +
      geom_line(size = 1) +
      geom_point(shape = 16, size = 5) +
      theme(
            panel.background = element_rect(fill = NA),
            axis.line = element_blank(),
            legend.position = "none",
            axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
            axis.text = element_text(family = "serif", colour = "black", size = 16),
            axis.title = element_text(family = "serif", colour = "black", size = 16)
      ) +
      scale_y_continuous(limits = c(0, 15)) +
      scale_color_manual(values = c("black")) +
      geom_rangeframe(
            data = data.frame(
                  X1gen_treatment = c("Control", "Herbivory"),
                  mean = c(0, 15)
            ),
            size = 1.5,
            color = "black"
      ) +
      geom_text(aes(x = 1.5, y = 14, label = " * "), family = "serif", size = 7) +
      labs(x = "G1 herbivory", y = "DNA methylation level \nof G1 (%)")

plot_5mc_1gen

# G2
# Group the data by X1gen_organ_origin and X1gen_treatment and summarize the variables
pd <- position_dodge(width = 0.2)

est_5mc_2gen <- chen_data_02gen %>%
      filter(!is.na(X5mc)) %>%
      dplyr::group_by(X1gen_organ_origin, X1gen_treatment) %>%
      dplyr::summarise(
            N    = length(X5mc),
            mean = mean(X5mc, na.rm = TRUE),
            sd   = sd(X5mc, na.rm = TRUE),
            se   = sd(X5mc, na.rm = TRUE) / sqrt(N)
      )

#
est_5mc_2gen$X1gen_organ_origin <- relevel(est_5mc_2gen$X1gen_organ_origin, "Primary root")
levels(est_5mc_2gen$X1gen_organ_origin)

# Generate the plot for G2
plot_5mc_2gen <- ggplot(
      est_5mc_2gen,
      aes(x = as.factor(X1gen_organ_origin), y = mean, group = factor(X1gen_treatment))
) +
      geom_errorbar(width = 0.1, aes(ymin = mean - se, ymax = mean + se), size = 1, position = pd) +
      geom_line(aes(linetype = factor(X1gen_treatment)), size = 1, position = pd) +
      geom_point(shape = 16, size = 5, position = pd) +
      theme(
            panel.background = element_rect(fill = NA),
            legend.title = element_text(family = "serif", colour = "black", size = 16),
            legend.text = element_text(family = "serif", colour = "black", size = 16),
            legend.position = c(0.4, 0.8),
            # legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
            axis.text = element_text(family = "serif", colour = "black", size = 16),
            axis.title = element_text(family = "serif", colour = "black", size = 16)
      ) +
      scale_y_continuous(limits = c(0, 10)) +
      geom_rangeframe(
            data = data.frame(
                  X1gen_treatment = c("Control", "Herbivory"),
                  X1gen_organ_origin = c("Primary root", "Secondary root"),
                  mean = c(0, 10)
            ),
            size = 1.5,
            color = "black"
      ) +
      labs(x = "G1 root order", y = "DNA methylation level \nof G2 (%)", linetype = "G1 herbivory")

plot_5mc_2gen

# G3
# Group the data by X1gen_organ_origin and X1gen_treatment and summarize the variables
pd <- position_dodge(width = 0.2)

est_5mc_3gen <- chen_data_03gen %>%
      filter(!is.na(X5mc)) %>%
      dplyr::group_by(X1gen_organ_origin, X1gen_treatment) %>%
      dplyr::summarise(
            N    = length(X5mc),
            mean = mean(X5mc, na.rm = TRUE),
            sd   = sd(X5mc, na.rm = TRUE),
            se   = sd(X5mc, na.rm = TRUE) / sqrt(N)
      )

#
est_5mc_3gen$X1gen_organ_origin <- relevel(est_5mc_3gen$X1gen_organ_origin, "Primary root")
levels(est_5mc_3gen$X1gen_organ_origin)

# Generate the plot for G3
plot_5mc_3gen <- ggplot(
      est_5mc_3gen,
      aes(x = as.factor(X1gen_organ_origin), y = mean, group = factor(X1gen_treatment))
) +
      geom_errorbar(width = .1, aes(ymin = mean - se, ymax = mean + se), size = 1, position = pd) +
      geom_line(aes(linetype = factor(X1gen_treatment)), size = 1, position = pd) +
      geom_point(shape = 16, size = 5, position = pd) +
      theme(
            panel.background = element_rect(fill = NA),
            legend.title = element_text(family = "serif", colour = "black", size = 16),
            legend.text = element_text(family = "serif", colour = "black", size = 16),
            legend.position = c(0.4, 0.8),
            # legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
            axis.text = element_text(family = "serif", colour = "black", size = 16),
            axis.title = element_text(family = "serif", colour = "black", size = 16)
      ) +
      scale_y_continuous(limits = c(0, 10)) +
      geom_rangeframe(
            data = data.frame(
                  X1gen_treatment = c("Control", "Herbivory"),
                  X1gen_organ_origin = c("Primary root", "Secondary root"),
                  mean = c(0, 10)
            ),
            size = 1.5,
            color = "black"
      ) +
      labs(x = "G1 root order", y = "DNA methylation level \nof G3 (%)", linetype = "G1 herbivory")

plot_5mc_3gen

# Create fig.6
fig.6 <- ggarrange(plot_5mc_1gen,
      plot_5mc_2gen,
      plot_5mc_3gen,
      labels = c("A", "B", "C"),
      ncol = 1,
      nrow = 3,
      align = "hv"
)

ggexport(fig.6,
      filename = "./Output/fig.6.png",
      width = 1500,
      height = 3000,
      pointsize = 14,
      res = 300
)
