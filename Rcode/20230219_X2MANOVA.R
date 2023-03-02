#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-06 10:42:08
# @Last Modified by:
# @Last Modified time: 2023-03-02 20:52:47
# @Description: MANOVAs for G2
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
chen_data_02gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "2nd_gen_20210305")

# Modify data format
chen_data_02gen <- chen_data_02gen %>%
      dplyr::mutate(
            X1gen_treatment = as.factor(X1gen_treatment),
            X1gen_organ_origin = as.factor(X1gen_organ_origin)
      )

str(chen_data_02gen)

# ==========================================================================
# MANOVA of growth variables
# ==========================================================================
growth.vars <- cbind(
      chen_data_02gen$leafmass,
      chen_data_02gen$stemmass,
      chen_data_02gen$rootmass,
      chen_data_02gen$numnode,
      chen_data_02gen$numleaf,
      chen_data_02gen$stolonlen
)

# Plots Mahalanobis distances against chi-squared distribution
coord <- qqplot(qchisq(ppoints(nrow(growth.vars)), df = ncol(growth.vars)), mahalanobis(growth.vars, colMeans(growth.vars), cov(growth.vars)))
abline(a = 0, b = 1)

# Performs Box's M test for homogeneity of covariance matrices
biotools::boxM(growth.vars, chen_data_02gen$X1gen_treatment)
biotools::boxM(growth.vars, chen_data_02gen$X1gen_organ_origin)

# Plots the location of observations in a 3D space defined by the first three principal components
mvoutlier::aq.plot(growth.vars)

# Specifies the contrasts for the ANOVA model
type2 <- list(X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)

# Fits the linear model using growth variables and treatment as the predictor variable
lm_growth <- lm(cbind(leafmass, stemmass, rootmass, numnode, numleaf, stolonlen) ~ X1gen_treatment * X1gen_organ_origin,
      data = chen_data_02gen,
      contrasts = type2
)

# Prints the summary of the ANOVA model
summary(lm_growth)

# Performs multivariate analysis of variance using the Pillai test statistic
x2manova_growth <- car::Manova(lm_growth, test.statistic = "Pillai", type = "III")
x2manova_growth

# ==========================================================================
# MANOVA of physiological variables
# ==========================================================================
chen_data_02gen01 <- chen_data_02gen %>% tidyr::drop_na()

physiology.vars <- cbind(
      chen_data_02gen01$phenolic,
      chen_data_02gen01$sugar,
      chen_data_02gen01$starch,
      chen_data_02gen01$carbohydrate
)


# Plots Mahalanobis distances against chi-squared distribution
coord <- qqplot(qchisq(ppoints(nrow(physiology.vars)), df = ncol(physiology.vars)), mahalanobis(physiology.vars, colMeans(physiology.vars), cov(physiology.vars)))
abline(a = 0, b = 1)

# Performs Box's M test for homogeneity of covariance matrices
biotools::boxM(physiology.vars, chen_data_02gen01$X1gen_treatment)
biotools::boxM(physiology.vars, chen_data_02gen01$X1gen_organ_origin)

# Plots the location of observations in a 3D space defined by the first three principal components
mvoutlier::aq.plot(physiology.vars)

# Specifies the contrasts for the ANOVA model
type2 <- list(X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)

# Fits the linear model using physiological variables and treatment as the predictor variable
lm_physiology <- lm(cbind(phenolic, sugar, starch, carbohydrate) ~ X1gen_treatment * X1gen_organ_origin,
      data = chen_data_02gen,
      contrasts = type2
)

# Prints the summary of the ANOVA model
summary(lm_physiology)

# Performs multivariate analysis of variance using the Pillai test statistic
x2manova_physiology <- car::Manova(lm_physiology, type = "III")
x2manova_physiology
