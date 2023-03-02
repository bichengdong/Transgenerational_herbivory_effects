#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-06 10:42:08
# @Last Modified by:
# @Last Modified time: 2023-03-02 20:53:34
# @Description: three-way ANOVA for G3
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
library(ggthemes)
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

# define the response variable names
response_vars <- c(
      "sqrt(leafmass)", "sqrt(stemmass)", "rootmass", "numnode", "numleaf", "stolonlen",
      "sqrt(totalmass)", "rsratio", "phenolic", "sugar", "starch", "carbohydrate"
)

# Levene's test results
levene_results <- response_vars %>%
      purrr::map(~ leveneTest(as.formula(paste(.x, "~ X3gen_treatment * X1gen_treatment * X1gen_organ_origin")), data = chen_data_03gen)) %>%
      purrr::set_names(response_vars) %>%
      map_dfr(broom::tidy, .id = "variable")

levene_results


# two-way ANOVAs
#
# Since Anova may report errors, it needs to be adjusted
# https://www.r-bloggers.com/2011/03/anova-%e2%80%93-type-iiiiii-ss-explained/
# https://rdoodles.rbind.io/2020/10/type-3-anova-in-r-an-easy-way-to-publish-wrong-tables/#:~:text=In%20R%2C%20so-called%20%E2%80%9CType%20I%20sums%20of%20squares%E2%80%9D,cannot%20simply%20specify%20car%3AAnova%20%28m1%2C%20type%20%3D%203%29.

type3 <- list(X3gen_treatment = contr.sum, X1gen_treatment = contr.sum, X1gen_organ_origin = contr.sum)

# use map to fit the linear models
lm_list <- response_vars %>%
      purrr::map(~ lm(as.formula(paste(.x, "~ X3gen_treatment * X1gen_treatment * X1gen_organ_origin")), data = chen_data_03gen, contrasts = type3)) %>%
      purrr::set_names(response_vars)

# use map to extract the ANOVA tables
x3aov_list <- lm_list %>%
      purrr::map(~ Anova(.x, type = "III"))

# covert list to data.frame
x3aov_df <- x3aov_list %>%
      map_dfr(broom::tidy, .id = "variable") %>%
      filter(!(term %in% c("(Intercept)", "Residuals"))) %>%
      mutate(
            sumsq = round(sumsq, 2),
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

x3aov_df

# Save results
save(levene_results,
      lm_list,
      x3aov_list,
      x3aov_df,
      file = "./Output/x3aov_results_3gen.Rdata"
)

write.csv(x3aov_df, "./Output/x3aov_df.csv", row.names = FALSE)

# -----------------------------------------------------------
# interaction
interaction_list <- lm_list %>%
      purrr::map(~ testInteractions(.x, fixed = "X1gen_organ_origin", across = "X1gen_treatment", adjustment = "none"))

interaction_df <- interaction_list %>%
      purrr::map(~ broom::tidy(.x)) %>%
      purrr::set_names(response_vars) %>%
      data.table::rbindlist(idcol = "variable") %>%
      mutate(
            Value = round(Value, 3),
            sumsq = round(sumsq, 2),
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

interaction_df

# --------------------------------------------------
# plot
# Define the function that generates a plot for a variable
generate_plot <- function(data, var, ylabel, ylimits) {

      # Group the data by X1gen_treatment and summarize the variables
      est_var <- data %>%
            group_by(X1gen_treatment, X1gen_organ_origin, X3gen_treatment) %>%
            dplyr::summarise(
                  N = n(),
                  mean = mean(.data[[var]], na.rm = TRUE),
                  sd = sd(.data[[var]], na.rm = TRUE),
                  se = sd(.data[[var]], na.rm = TRUE) / sqrt(N)
            ) %>%
            mutate(X1gen_organ_origin = relevel(X1gen_organ_origin, "Primary root"))

      pd <- position_dodge(width = 0.2)

      # Generate the plot
      ggplot(
            est_var,
            aes(x = as.factor(X3gen_treatment), y = mean, group = factor(X1gen_treatment))
      ) +
            geom_errorbar(width = 0.1, aes(ymin = mean - se, ymax = mean + se), size = 1, position = pd) +
            geom_line(aes(linetype = factor(X1gen_treatment)), size = 1, position = pd) +
            geom_point(shape = 16, size = 5, position = pd) +
            facet_wrap(~X1gen_organ_origin) +
            theme(
                  panel.background = element_rect(fill = NA),
                  legend.title = element_text(family = "serif", colour = "black", size = 16),
                  legend.text = element_text(family = "serif", colour = "black", size = 16),
                  legend.position = "none",
                  strip.text.x = element_text(family = "serif", colour = "black", size = 16),
                  axis.line = element_blank(),
                  axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                  axis.text = element_text(family = "serif", colour = "black", size = 16),
                  axis.title = element_text(family = "serif", colour = "black", size = 16)
            ) +
            scale_y_continuous(limits = ylimits) +
            scale_color_manual(values = c("black")) +
            labs(x = "G3 herbivory", y = ylabel, linetype = "G1 herbivory") +
            geom_rangeframe(
                  data = data.frame(est_var[, 1:3],
                        mean = c(ylimits[1], ylimits[2], ylimits[1], ylimits[1], ylimits[1], ylimits[2], ylimits[1], ylimits[1])
                  ),
                  size = 1.5,
                  color = "black"
            )
}

# Create features for each variable in data.frame
var.df <- data.frame(
      vars = c("totalmass", "numnode", "leafmass", "numleaf", "stemmass", "stolonlen", "rootmass", "rsratio", "phenolic", "sugar", "starch", "carbohydrate"),
      ylabs = c("Total mass (g)", "Number of nodes", "Leaf mass (g)", "Number of leaves", "Stem mass (g)", "Stolon length (cm)", "Root mass (g)", "Root to shoot ratio", "Total phenolics (%)", "Soluble sugar (%)", "Starch (%)", "Total NSC (%)"),
      ylimits01 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
      ylimits02 = c(15, 250, 5, 600, 8, 1000, 4, 1.5, 4, 6, 8, 10)
)

# Generate plots using the generate_plot function for each variable in the df
plots <- set_names(
      map(var.df$vars, ~ generate_plot(data = chen_data_03gen, var = .x, ylabel = var.df$ylabs[var.df$vars %in% .x], ylimits = c(var.df$ylimits01[var.df$vars %in% .x], var.df$ylimits02[var.df$vars %in% .x]))),
      var.df$vars
)

plots

# Add the legends to the plots
# plots$totalmass <- plots$totalmass + theme(legend.position = c(0.3, 0.2))
# plots$phenolic <- plots$phenolic + theme(legend.position = c(0.3, 0.8))
legends <- tibble(vars = c("totalmass", "phenolic"), position = list(c(0.4, 0.8), c(0.4, 0.8)))
plots[legends$vars] <- map2(plots[legends$vars], legends$position, ~ .x + theme(legend.position = .y))


# Create fig.4
fig.4 <- ggarrange(
      plots$totalmass,
      plots$numnode,
      plots$leafmass,
      plots$numleaf,
      plots$stemmass,
      plots$stolonlen,
      plots$rootmass,
      plots$rsratio,
      labels = c("A", "E", "B", "F", "C", "G", "D", "H"),
      ncol = 2,
      nrow = 4,
      align = "hv"
)

fig.4

# Export fig.4
ggexport(fig.4,
      filename = "./Output/fig.4.png",
      width = 3000,
      height = 4000,
      pointsize = 14,
      res = 300
)


#  Create fig.5
fig.5 <- ggarrange(
      plots$phenolic,
      plots$starch,
      plots$sugar,
      plots$carbohydrate,
      labels = c("A", "C", "B", "D"),
      ncol = 2,
      nrow = 4,
      align = "hv"
)

fig.5

# Export fig.5
ggexport(fig.5,
      filename = "./Output/fig.5.png",
      width = 3000,
      height = 4000,
      pointsize = 14,
      res = 300
)
