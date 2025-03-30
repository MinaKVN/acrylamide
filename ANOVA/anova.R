library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcompView)
library(RVAideMemoire)
library(dplyr)


# Load data
data <- read.csv("acr_all.csv")
str(data)

data$id <- as.factor(data$id)
data$location    <- as.factor(data$location)
data$asparagine <- as.numeric(data$asparagine)

## what assumptions should we be checking?
# 1. Residuals are random - no relationship to our treatment
# 2. Homogeneity of residuals across our treatment groups
# 3. Residuals are normally distributed
# 4. Residuals have a mean of 0)
# Fit the ANOVA model

## asparagine
#  want to model the interaction between location and genotype
asa_model <- aov(
  asparagine ~ location + id + location * id,
  data = data)

# Summary of the ANOVA model
summary(asa_model)
anova(asa_model)

# Residual Plots 
par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(asa_model)

# Shapiro Wilk 
shapiro.test(rstandard(asa_model))

# Tukey's test
tukey <- TukeyHSD(asa_model)
print(tukey)
cld <- multcompLetters4(asa_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk <- group_by(data, id) %>%
  summarise(mean=mean(asparagine), quant = quantile(asparagine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk$cld <- cld$Letters

print(Tk)

ggplot(data, aes(id, asparagine)) + 
  geom_boxplot() +
  labs(x="genotype", y="asparagine") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = id, y = quant, label = cld), size = 3, vjust=-1, hjust =-1)
