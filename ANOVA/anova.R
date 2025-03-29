library(stringr)
library(purrr)
library(janitor)
library(MuMIn)
library(broom)
library(forcats)
library(tidyverse)
library(arm)
library(vegan)
library(ggvegan)
library(dplyr)
library(visreg)
library(car)
library(reshape2)
library(MASS)

# Load data
data <- read.csv("acr_all.csv")
str(data)

## what assumptions should we be checking?
# 1. Residuals are random - no relationship to our treatment
# 2. Homogeneity of residuals across our treatment groups
# 3. Residuals are normally distributed
# 4. Residuals have a mean of 0)
# Fit the ANOVA model


# Assuming 'data' is your dataset and you want to model the interaction between location and genotype
aov_model <- aov(
  asparagine ~ location + genotype + location * genotype,
  data = data)

# Summary of the ANOVA model
summary(aov_model)

#Randomness of residuals
plot(residuals(aov_model) ~ fitted(aov_model), xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
#Homogeneity of Residuals (Homoscedasticity)
leveneTest(residuals(aov_model), group = data$genotype)

#Normality of Residuals
qqnorm(residuals(aov_model))
qqline(residuals(aov_model), col = "red")
shapiro.test(residuals(aov_model))
mean(residuals(aov_model))
