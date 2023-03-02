#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-05 20:42:50
# @Last Modified by:
# @Last Modified time: 2023-03-02 20:25:22
# @Description: t.test for G1
#--------------------------------------------------------------------------#
# Clear memory
rm(list = ls())
gc()

# Load packages as needed
library(broom)
library(car)
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(purrr)
library(dplyr)

# Set working directory
setwd("G:/我的坚果云/002项目/2022_Chenlin/投稿期刊plants/data summary")
getwd()

# Load data
chen_data_01gen <- read_excel("./Data/20230227data_Chenlin.xlsx", sheet = "1st_gen_20210305")

# Modify data format
chen_data_01gen <- chen_data_01gen %>%
      dplyr::mutate(X1gen_treatment = as.factor(X1gen_treatment))

# define the response variable names
response_vars <- c(
      "leafmass", "stemmass", "rootmass", "numnode", "numleaf", "stolonlen",
      "totalmass", "rsratio", "phenolic", "sugar", "starch", "carbohydrate"
)

# Levene's test results
levene_results <- response_vars %>%
      purrr::map(~ leveneTest(as.formula(paste(.x, "~ X1gen_treatment")), data = chen_data_01gen)) %>%
      purrr::set_names(response_vars) %>%
      map_dfr(broom::tidy, .id = "variable")

levene_results

# use map to fit the linear models
t.test_list <- response_vars %>%
      purrr::map(~ t.test(as.formula(paste(.x, "~ X1gen_treatment ")), var.equal = TRUE, conf.level = 0.95, data = chen_data_01gen)) %>%
      purrr::set_names(response_vars)

# covert list to data.frame
t.test_df <- t.test_list %>%
      map_dfr(broom::tidy, .id = "variable") %>%
      mutate(
            statistic = round(statistic, 2),
            p.value = round(p.value, 3)
      ) %>%
      mutate(symbol = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 & p.value > 0.001 ~ "**",
            p.value < 0.05 & p.value > 0.01 ~ "*",
            p.value > 0.05 ~ "ns",
            is.na(p.value) ~ ""
      ))


# Save results
save(levene_results, t.test_list, t.test_df, file = "./Output/t_test_results_1gen.Rdata")

# ==========================================================================
# plot
# ==========================================================================
# Create a function to generate the plots for different variables
# usage of double bracket
# relevant ref
# https://blog.csdn.net/weixin_42514521/article/details/112438066

# Define the function that generates a plot for a variable
generate_plot <- function(data, var, ylabel, ylimits) {
      if (is.null(data) || !("X1gen_treatment" %in% names(data))) {
            stop("Input data is empty or does not have X1gen_treatment variable")
      }

      # Group the data by X1gen_treatment and summarize the variables
      est_var <- data %>%
            group_by(X1gen_treatment) %>%
            dplyr::summarise(
                  N = n(),
                  mean = mean(.data[[var]], na.rm = TRUE),
                  sd = sd(.data[[var]], na.rm = TRUE),
                  se = sd(.data[[var]], na.rm = TRUE) / sqrt(N)
            )

      # Generate the plot
      ggplot(
            est_var,
            aes(
                  x = as.factor(X1gen_treatment),
                  y = mean,
                  group = 1,
                  color = "black"
            )
      ) +
            geom_errorbar(aes(
                  ymin = mean - se,
                  ymax = mean + se
            ),
            width = 0.1,
            size = 1
            ) +
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
            scale_y_continuous(limits = ylimits) +
            scale_color_manual(values = c("black")) +
            labs(x = "", y = ylabel) +
            geom_rangeframe(
                  data = data.frame(
                        X1gen_treatment = c("Control", "Herbivory"),
                        mean = c(ylimits[1], ylimits[2])
                  ),
                  size = 1.5,
                  color = "black"
            )
}


# plot features for each variable
var.df <- list(
      totalmass = data.frame(ylabs = "Total mass (g)", ylimits01 = 0, ylimits02 = 20),
      numnode = data.frame(ylabs = "Number of nodes", ylimits01 = 0, ylimits02 = 200),
      leafmass = data.frame(ylabs = "Leaf mass (g)", ylimits01 = 0, ylimits02 = 6),
      numleaf = data.frame(ylabs = "Number of leaves", ylimits01 = 0, ylimits02 = 500),
      stemmass = data.frame(ylabs = "Stem mass (g)", ylimits01 = 0, ylimits02 = 6),
      stolonlen = data.frame(ylabs = "Stolon length (cm)", ylimits01 = 0, ylimits02 = 1000),
      rootmass = data.frame(ylabs = "Root mass (g)", ylimits01 = 0, ylimits02 = 4),
      rsratio = data.frame(ylabs = "Root to shoot ratio", ylimits01 = 0, ylimits02 = 1),
      phenolic = data.frame(ylabs = "Total phenolics (%)", ylimits01 = 0, ylimits02 = 4),
      sugar = data.frame(ylabs = "Soluble sugar (%)", ylimits01 = 0, ylimits02 = 8),
      starch = data.frame(ylabs = "Starch (%)", ylimits01 = 0, ylimits02 = 8),
      carbohydrate = data.frame(ylabs = "Total NSC (%)", ylimits01 = 0, ylimits02 = 15)
)

var.df <- map(names(var.df), ~ mutate(var.df[[.x]], symbol = t.test_df[t.test_df$variable == .x, ]$symbol)) %>%
      set_names(names(var.df))

# Generate plots using the generate_plot function for each variable in the df
plots <- map(names(var.df), ~ generate_plot(data = chen_data_01gen, var = .x, ylabel = var.df[[.x]]$ylabs, ylimits = c(var.df[[.x]]$ylimits01, var.df[[.x]]$ylimits02))) %>%
      set_names(names(var.df))

plots

# Add the t-test results to the plots and create figures
plots02 <- map(names(var.df), ~ plots[[.x]] + geom_text(data = var.df[[.x]], aes(x = 1.5, y = ylimits02, label = symbol), family = "serif", size = 6)) %>%
      set_names(names(var.df))

plots02

# Create appendix.fig.2
appendix.fig.2 <- ggarrange(
      plots02$totalmass,
      plots02$numnode,
      plots02$leafmass,
      plots02$numleaf,
      plots02$stemmass,
      plots02$stolonlen,
      plots02$rootmass,
      plots02$rsratio,
      labels = c("A", "E", "B", "F", "C", "G", "D", "H"),
      ncol = 2,
      nrow = 4,
      align = "hv"
)

appendix.fig.2

# Export appendix.fig.2
ggexport(appendix.fig.2,
      filename = "./Output/appendix.fig.2.png",
      width = 3000,
      height = 4000,
      pointsize = 14,
      res = 300
)


#  Create appendix.fig.3
appendix.fig.3 <- ggarrange(
      plots02$phenolic,
      plots$starch,
      plots$sugar,
      plots02$carbohydrate,
      labels = c("A", "C", "B", "D"),
      ncol = 2,
      nrow = 4,
      align = "hv"
)

appendix.fig.3

# Export appendix.fig.3
ggexport(appendix.fig.3,
      filename = "./Output/appendix.fig.3.png",
      width = 3000,
      height = 4000,
      pointsize = 14,
      res = 300
)
