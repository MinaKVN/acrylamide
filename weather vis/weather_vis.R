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

#grain filling period

data %>%
  ggplot(aes(x = location, y = gfp, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.2,
    color = "black",
    size = 1.3,
    alpha = 0.1
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
  theme(
    legend.position = "none",
    axis.title.x = element_text(
      size = 11,
      color = "grey10",
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 11,
      color = "grey10",
      face = "bold"
    ),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  ) +
  xlab("environment") +
  ylab("growing filling period (days)")

# Save the plot
ggsave("locationxGFP.png", units="in", width=5.5, height=3.75, dpi=800)

#growing degree days

data %>%
  ggplot(aes(x = location, y = total_GDD, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.2,
    color = "black",
    size = 1.3,
    alpha = 0.1
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
  theme(
    legend.position = "none",
    axis.title.x = element_text(
      size = 11,
      color = "grey10",
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 11,
      color = "grey10",
      face = "bold"
    ),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  ) +
  xlab("environment") +
  ylab("growing degree days")

# Save the plot
ggsave("locationxGDD.png", units="in", width=5.5, height=3.75, dpi=800)
