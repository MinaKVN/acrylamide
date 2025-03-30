
#if you hav not started using git and r studio yet run this first
library(usethis)
use_git_config(user.name = "MinaKVN", user.email = "minakavyani@gmail.com")

#loading packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcompView)
library(RVAideMemoire)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggplot2)


# Load data
data <- read.csv("acr_all.csv")
str(data)
data$id <- as.factor(data$id)
data$location    <- as.factor(data$location)
data$asparagine <- as.numeric(data$asparagine)


# Update the values in the 'location' column based on the old values
data <- data %>%
  mutate(location = case_when(
    location == "ea" ~ "Elora",
    location == "oa" ~ "Ottawa",
    location == "pa" ~ "Palmerston",
    location == "wr" ~ "Winchester",
    TRUE ~ location  # Default to keep original value if no match
  ))

 
## what assumptions should we be checking?
# 1. Residuals are random - no relationship to our treatment
# 2. Homogeneity of residuals across our treatment groups
# 3. Residuals are normally distributed
# 4. Residuals have a mean of 0)
# Fit the ANOVA model
## asparagine

#  want to model the interaction between location and genotype
asa_model <- aov(asparagine ~ location + id + location * id, data = data)

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

# Plot by genotype_ASA with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(asparagine),
            quant = quantile(asparagine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = asparagine, fill = id)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.2,
    color = "black",
    size = 1.3,
    alpha = 0.2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 3,
    color = "red3"
  ) +  # Adding mean as a point
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size = 11)) +
  xlab("genotype") +
  ylab("asparagine concentration (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -5,
    hjust = 0.5)

ggsave(
  "asaxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location ASA

# Tukey's test
tukey <- TukeyHSD(asa_model)
print(tukey)
cld <- multcompLetters4(asa_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(asparagine),
            quant = quantile(asparagine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = asparagine, fill = location)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.2,
    color = "black",
    size = 1.3,
    alpha = 0.2
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 3,
    color = "red3"
  ) +  # Adding mean as a point
  scale_fill_viridis(discrete = TRUE,
                     alpha = 0.6,
                     option = "turbo") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size = 11)) + #labs(fill = "Location")+
  xlab("location") + ylab ("asparagine concentration (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -5,
    hjust = -0.2
  )

# saving the graph in Tiff format. change the name accordingly

ggsave(
  "asaxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype asa
# Tukey's test

tukey <- TukeyHSD(asa_model)
print(tukey)
cld <- multcompLetters4(asa_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(asparagine),
            quant = quantile(asparagine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = asparagine, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("asparagine concentration (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -3.5,
    hjust = 0.5
  )

# saving the graph in Tiff format. change the name accordingly

ggsave(
  "glycinexGxE.tiff",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)

#https://statdoe.com/barplot-for-two-factors-in-r/

