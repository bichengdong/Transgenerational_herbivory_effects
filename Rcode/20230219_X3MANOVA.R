#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-06 10:42:08
# @Last Modified by:
# @Last Modified time: 2023-03-02 20:31:48
# @Description: MAnovas for G3
#--------------------------------------------------------------------------#
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

# Load data
chen_data_03gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "3rd_gen_20210305")

chen_data_03gen <- data.frame(chen_data_03gen)
glimpse(chen_data_03gen)
View(chen_data_03gen)

# Modify data format
chen_data_03gen <- chen_data_03gen %>%
      dplyr::mutate(
            X1gen_treatment = as.factor(X1gen_treatment),
            X3gen_treatment = as.factor(X3gen_treatment),
            X1gen_organ_origin = as.factor(X1gen_organ_origin)
      )

str(chen_data_03gen)

# ==========================================================================
# MANOVA of growth variables
# ==========================================================================
growth.vars <- cbind(sqrt(chen_data_03gen$leafmass),
                     sqrt(chen_data_03gen$stemmass),
                     chen_data_03gen$rootmass,
                     chen_data_03gen$numnode,
                     chen_data_03gen$numleaf,
                     chen_data_03gen$stolonlen)

# Plots Mahalanobis distances against chi-squared distribution
coord <- qqplot(qchisq(ppoints(nrow(growth.vars)), df = ncol(growth.vars)), mahalanobis(growth.vars, colMeans(growth.vars), cov(growth.vars)))
abline(a = 0, b = 1)

# Performs Box's M test for homogeneity of covariance matrices
biotools::boxM(growth.vars, chen_data_03gen$X3gen_treatment)
biotools::boxM(growth.vars, chen_data_03gen$X1gen_treatment)
biotools::boxM(growth.vars, chen_data_03gen$X1gen_organ_origin)

# Plots the location of observations in a 3D space defined by the first three principal components
mvoutlier::aq.plot(physiology.vars)

# Specifies the contrasts for the ANOVA model
type3 <- list(X3gen_treatment = contr.sum, X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)

# Fits the linear model using growth variables and treatment as the predictor variable
lm_growth <-lm(cbind(sqrt(leafmass), sqrt(stemmass), rootmass, numnode, numleaf, stolonlen) ~ X1gen_treatment * X1gen_organ_origin * X3gen_treatment,
          data = chen_data_03gen,
          contrasts = type3)

# Prints the summary of the ANOVA model
summary(lm_growth)

# Performs multivariate analysis of variance using the Pillai test statistic
x3manova_growth <- car::Manova(lm_growth, test.statistic = "Pillai", type = "III")
x3manova_growth

# ==========================================================================
# MANOVA of physiological variables
# ==========================================================================
chen_data_03gen01 <- chen_data_03gen %>% tidyr::drop_na()

physiology.vars <- cbind(chen_data_03gen01$phenolic,
                       chen_data_03gen01$sugar,
                       chen_data_03gen01$starch,
                       chen_data_03gen01$carbohydrate)


# Plots Mahalanobis distances against chi-squared distribution
coord <- qqplot(qchisq(ppoints(nrow(physiology.vars)), df = ncol(physiology.vars)), mahalanobis(physiology.vars, colMeans(physiology.vars), cov(physiology.vars)))
abline(a = 0, b = 1)

# Performs Box's M test for homogeneity of covariance matrices
biotools::boxM(physiology.vars, chen_data_03gen01$X3gen_treatment)
biotools::boxM(physiology.vars, chen_data_03gen01$X1gen_treatment)
biotools::boxM(physiology.vars, chen_data_03gen01$X1gen_organ_origin)

# Plots the location of observations in a 3D space defined by the first three principal components
type3 <- list(X3gen_treatment = contr.sum, X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)

# Fits the linear model using physiological variables and treatment as the predictor variable
lm_physiology <-lm(cbind(phenolic, sugar, starch, carbohydrate) ~ X1gen_treatment * X1gen_organ_origin * X3gen_treatment,
          data = chen_data_03gen,
          contrasts = type3)

# Prints the summary of the ANOVA model
summary(lm_physiology)

# Performs multivariate analysis of variance using the Pillai test statistic
x3manova_physiology <- car::Manova(lm_physiology, type = "III")
x3manova_physiology
