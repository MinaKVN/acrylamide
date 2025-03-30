
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
    location == "pn" ~ "Palmerston",
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
  ylab("asparagine (μg/g)") +
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
  xlab("location") + ylab ("asparagine (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -5,
    hjust = -0.2
  )

# saving the graph in png format. change the name accordingly

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
  xlab("location") + ylab ("asparagine (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -3.5,
    hjust = 0.5
  )

# saving the graph in png format. change the name accordingly

ggsave(
  "asaxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)
#this is a good website for visulaization of ANOVA results
#https://statdoe.com/barplot-for-two-factors-in-r/

####################################
#####glutamine anova and vis########
####################################

str(data)
glu_model <- aov(glutamine ~ location + id + location * id, data = data)

# Summary of the ANOVA model
summary(glu_model)
anova(glu_model)

# Residual Plots

par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(glu_model)

# Shapiro Wilk

shapiro.test(rstandard(glu_model))

# Tukey's test
tukey <- TukeyHSD(glu_model)
print(tukey)
cld <- multcompLetters4(glu_model, tukey)
print(cld)

# Plot by genotype_glu with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(glutamine),
            quant = quantile(glutamine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = glutamine, fill = id)) +
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
  ylab("glutamine (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -5,
    hjust = -1)

ggsave(
  "gluxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location glu

# Tukey's test
tukey <- TukeyHSD(glu_model)
print(tukey)
cld <- multcompLetters4(glu_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(glutamine),
            quant = quantile(glutamine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = glutamine, fill = location)) +
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
  xlab("location") + ylab ("glutamine (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -5,
    hjust = -0.2
  )

# saving the graph in png format. change the name accordingly

ggsave(
  "gluxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype glu
# Tukey's test

tukey <- TukeyHSD(glu_model)
print(tukey)
cld <- multcompLetters4(glu_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(glutamine),
            quant = quantile(glutamine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = glutamine, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("glutamine (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -1.65,
    hjust = 0.58
  )

# saving the graph in png format. change the name accordingly

ggsave(
  "gluxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)

####################################
#####glycine anova and vis########
####################################

str(data)
gly_model <- aov(glycine ~ location + id + location * id, data = data)

# Summary of the ANOVA model
summary(gly_model)
anova(gly_model)

# Residual Plots

par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(gly_model)

# Shapiro Wilk

shapiro.test(rstandard(gly_model))

# Tukey's test
tukey <- TukeyHSD(gly_model)
print(tukey)
cld <- multcompLetters4(gly_model, tukey)
print(cld)

# Plot by genotype_gly with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(glycine),
            quant = quantile(glycine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = glycine, fill = id)) +
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
  ylab("glycine (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -2,
    hjust = -1)

ggsave(
  "glyxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location gly

# Tukey's test
tukey <- TukeyHSD(gly_model)
print(tukey)
cld <- multcompLetters4(gly_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(glycine),
            quant = quantile(glycine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = glycine, fill = location)) +
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
  xlab("location") + ylab ("glycine (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -2.5,
    hjust = -0.3
  )

# saving the graph in png format. change the name accordingly

ggsave(
  "glyxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype gly
# Tukey's test

tukey <- TukeyHSD(gly_model)
print(tukey)
cld <- multcompLetters4(gly_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(glycine),
            quant = quantile(glycine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = glycine, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("glycine (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -1.60,
    hjust = -0.2
  )

# saving the graph in png format. change the name accordingly

ggsave(
  "glyxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)

####################################
#####lysine anova and vis########
####################################

str(data)
#making Lysine to lysine
data<- data %>% rename(lysine = Lysine)

lys_model <- aov(
  lysine ~ location + id + location * id, data = data)

# Summary of the ANOVA model
summary(lys_model)
anova(lys_model)

# Residual Plots

par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(lys_model)

# Shapiro Wilk

shapiro.test(rstandard(lys_model))

# Tukey's test
tukey <- TukeyHSD(lys_model)
print(tukey)
cld <- multcompLetters4(lys_model, tukey)
print(cld)

# Plot by genotype_lys with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(lysine),
            quant = quantile(lysine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = lysine, fill = id)) +
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
  ylab("lysine (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -2,
    hjust = -1)

ggsave(
  "lysxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location lys

# Tukey's test
tukey <- TukeyHSD(lys_model)
print(tukey)
cld <- multcompLetters4(lys_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(lysine),
            quant = quantile(lysine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = lysine, fill = location)) +
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
  xlab("location") + ylab ("lysine (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -2.5,
    hjust = -0.3
  )

# saving the graph in png format. change the name accordinlys

ggsave(
  "lysxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype lys
# Tukey's test

tukey <- TukeyHSD(lys_model)
print(tukey)
cld <- multcompLetters4(lys_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(lysine),
            quant = quantile(lysine, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = lysine, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("lysine (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -1.5,
    hjust = -0.1
  )

# saving the graph in png format. change the name accordinlys

ggsave(
  "lysxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)


####################################
#####proline anova and vis########
####################################

str(data)
prol_model <- aov(proline ~ location + id + location * id, data = data)

# Summary of the ANOVA model
summary(prol_model)
anova(prol_model)

# Residual Plots

par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(prol_model)

# Shapiro Wilk

shapiro.test(rstandard(prol_model))

# Tukey's test
tukey <- TukeyHSD(prol_model)
print(tukey)
cld <- multcompLetters4(prol_model, tukey)
print(cld)

# Plot by genotype_prol with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(proline),
            quant = quantile(proline, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = proline, fill = id)) +
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
  ylab("proline (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -2,
    hjust = -1)

ggsave(
  "prolxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location prol

# Tukey's test
tukey <- TukeyHSD(prol_model)
print(tukey)
cld <- multcompLetters4(prol_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(proline),
            quant = quantile(proline, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = proline, fill = location)) +
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
  xlab("location") + ylab ("proline (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -2.5,
    hjust = -0.3
  )

# saving the graph in png format. change the name accordinprol

ggsave(
  "prolxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype prol
# Tukey's test

tukey <- TukeyHSD(prol_model)
print(tukey)
cld <- multcompLetters4(prol_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(proline),
            quant = quantile(proline, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = proline, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("proline (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -1.2,
    hjust = -0.35
  )

# saving the graph in png format. change the name accordinprol

ggsave(
  "prolxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)

####################################
#####protein anova and vis########
####################################

str(data)
prot_model <- aov(protein ~ location + id + location * id, data = data)

# Summary of the ANOVA model
summary(prot_model)
anova(prot_model)

# Residual Plots

par(mfrow = c(2, 2)) # Split the plotting panel into a 2 x 2 grid plot(model1)
plot(prot_model)

# Shapiro Wilk

shapiro.test(rstandard(prot_model))

# Tukey's test
tukey <- TukeyHSD(prot_model)
print(tukey)
cld <- multcompLetters4(prot_model, tukey)
print(cld)

# Plot by genotype_prot with mean line

Tk_G <- group_by(data, id) %>%
  summarise(mean = mean(protein),
            quant = quantile(protein, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$id)
Tk_G$cld <- cld$Letters
print(Tk_G)

data %>%
  ggplot(aes(x = id, y = protein, fill = id)) +
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
  ylab("protein (μg/g)") +
  geom_text(
    data = Tk_G,
    aes(x = id, y = quant, label = cld),
    size = 5,
    color = "blue4", 
    vjust = -5,
    hjust = 0.5)

ggsave(
  "protxG.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location prot

# Tukey's test
tukey <- TukeyHSD(prot_model)
print(tukey)
cld <- multcompLetters4(prot_model, tukey)
print(cld)


# table with factors and 3rd quantile
Tk_E <- group_by(data, location) %>%
  summarise(mean = mean(protein),
            quant = quantile(protein, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$location)
Tk_E$cld <- cld$Letters
print(Tk_E)
data %>%
  ggplot(aes(x = location, y = protein, fill = location)) +
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
  xlab("location") + ylab ("protein (μg/g)") +
  geom_text(
    data = Tk_E,
    aes(x = location, y = quant, label = cld),
    size = 5,
    color = "blue4",
    vjust = -3,
    hjust = -0.3
  )

# saving the graph in png format. change the name accordinprot

ggsave(
  "protxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 800
)

# Plot by location and genotype prot
# Tukey's test

tukey <- TukeyHSD(prot_model)
print(tukey)
cld <- multcompLetters4(prot_model, tukey)
print(cld)
# table with factors and 3rd quantile
Tk_GE <- data %>%
  group_by(location, id) %>%
  summarise(mean = mean(protein),
            quant = quantile(protein, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- as.data.frame.list(cld$'location:id')
Tk_GE$cld <- cld$Letters
print(Tk_GE)
data %>%
  ggplot(aes(x = location, y = protein, fill = id)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 1)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(size = 11)) + labs(fill = "genotype") +
  xlab("location") + ylab ("protein (μg/g)") +
  geom_text(
    data = Tk_GE,
    aes(x = location, y = quant, label = cld),
    position = position_dodge(1),
    size = 4,
    color = "blue4",
    vjust = -2,
    hjust = -0.1
  )

# saving the graph in png format. change the name accordinprot

ggsave(
  "protxGxE.png",
  units = "in",
  width = 5.5,
  height = 3.75,
  dpi = 700
)





