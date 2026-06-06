library(MASS)
library(survival)
library(car)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(tidyverse)

life_data <- read.table("MK_LRS_genome_data.csv", h=T, sep=",")

life_data$het_scaled  <- as.vector(scale(life_data$het))
life_data$roh_scaled <- as.vector(scale(life_data$roh))
life_data$snpeff_scaled  <- as.vector(scale(life_data$hom_load_SnpEff))
life_data$cadd_scaled  <- as.vector(scale(life_data$hom_load_CADD))
model <- glm.nb(LRS ~ het_scaled + breeding_years + sex_genome, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
summary(model)
model <- glm.nb(LRS ~ roh_scaled + breeding_years + sex_genome, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
summary(model)
model <- glm.nb(LRS ~ snpeff_scaled + breeding_years + sex_genome, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
summary(model)
model <- glm.nb(LRS ~ cadd_scaled + breeding_years + sex_genome, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
summary(model)

plot_prediction <- function(glm_model, genetics) {
newdata <- life_data %>%
  summarise(min_het = min(.data[[genetics]]), max_het = max(.data[[genetics]]), median_years = median(breeding_years)) %>%
  rowwise() %>%
  do(data.frame(value = seq(.$min_het, .$max_het, length.out = 100), breeding_years = .$median_years))
names(newdata)[names(newdata) == "value"] <- genetics

# Predict on the new data
newdata$predicted <- predict(glm_model, newdata = newdata, type = "response")

# For confidence intervals, predict on link scale then transform
pred_link <- predict(glm_model, newdata = newdata, type = "link", se.fit = TRUE)
newdata <- newdata %>%
  mutate(
    fit_link = pred_link$fit,
    se_link = pred_link$se.fit,
    lower = exp(fit_link - 1.96 * se_link),
    upper = exp(fit_link + 1.96 * se_link)
  )

# Plot using ggplot2
ggplot(newdata, aes(x = .data[[genetics]], y = predicted)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey70", color = NA) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )
}

model_het <- glm.nb(LRS ~ het + breeding_years, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
model_roh <- glm.nb(LRS ~ roh + breeding_years, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
model_snpeff <- glm.nb(LRS ~ hom_load_SnpEff + breeding_years, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)
model_cadd <- glm.nb(LRS ~ hom_load_CADD + breeding_years, data = life_data[which(life_data$bred_in_2024==0),], na.action = na.omit)

plot_prediction(model_het, "het") +
  labs(
    x = "Heterozygosity",
    y = "Predicted Lifetime Reproductive Success",
    title = "(a) Heterozygosity"
  )

plot_prediction(model_roh, "roh") +
  labs(
    x = "Runs of homozygosity",
    y = "Predicted Lifetime Reproductive Success",
    title = "(b) Runs of Homozygosity"
  )

plot_prediction(model_snpeff, "hom_load_SnpEff") +
  labs(
    x = "Homozygous load (SnpEff)",
    y = "Predicted Lifetime Reproductive Success",
    title = "(c) Homozygous load (SnpEff)"
  )

plot_prediction(model_cadd, "hom_load_CADD") +
  labs(
    x = "Homozygous load (CADD)",
    y = "Predicted Lifetime Reproductive Success",
    title = "Homozygous load (CADD)"
  )
