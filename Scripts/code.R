# EM HANSSON April 2022

# Load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(arm)
library(pbkrtest)

# Data import
df <- read_csv("Data/df.csv") 

df_mod <- df %>%
  mutate(exp_day2 = exp_day - mean(exp_day)) %>% # Center day variable for model fit
  filter(exp_day > 7) %>% # Remove establishment phase for better fit
  filter(!(test == "D" & chamber == "13")) %>% # Remove chamber with leak
  filter(exp_day < 51)  # Remove post-contamination phase

# Models
## Straight line through data, exp_day just a random effect, no effect of rpm or of test
mod_rpm <- lmer(log(cell_count) ~ rpm + ( 1 | chamberID) + (1 | exp_day2), data = df_mod, 
                    control = lmerControl(optimizer="Nelder_Mead"))

mod_test <- lmer(log(cell_count) ~ test + ( 1 | chamberID) + (1 | exp_day2), data = df_mod, 
                             control = lmerControl(optimizer="Nelder_Mead"))

# Summaries and residuals
summary(mod_rpm)
Anova(mod_rpm)
plot(mod_rpm)

summary(mod_test)
Anova(mod_test)
plot(mod_test)

# Reduced model and signficance testing for robustness
mod_red <- lmer(log(cell_count) ~ ( 1 | chamberID) + (1 | exp_day2), data = data_relaxed, 
                                     control = lmerControl(optimizer="Nelder_Mead"))

anova(mod_rpm, mod_red)
PBmodcomp(mod_rpm, mod_red, nsim = 1000)

anova(mod_test, mod_red)
PBmodcomp(mod_test, mod_red, nsim = 1000)

# Plots included in chapter

### Import data
# Population density
df_density_sum <- df %>%
  filter(!(test == "D" & chamber == "13")) %>% # Remove chamber with leak
  filter(exp_day < 51) %>%  # Remove post-contamination phase
  group_by(test, rpm, exp_day) %>%
  summarise(
    mean_cell_count = mean(cell_count),
    se_cell_count = (sd(cell_count)/sqrt(n()))) %>%
  ungroup()

ggplot(df %>% mutate(mean_cell_count = cell_count)  %>% filter(!(test == "D" & chamber == "13")) %>%
         filter(exp_day < 51), aes(x = exp_day, y = log(mean_cell_count), colour = test, linetype = rpm)) +
  geom_point() +
  geom_smooth(data = df_density_sum, aes(ymax = log(mean_cell_count + se_cell_count), ymin = log(mean_cell_count - se_cell_count)), stat = "identity") +
  scale_x_continuous(breaks = seq(1, 40, 2), labels = seq(1, 40, 2), limits = c(1, 40), name = "Days since inoculation") +
  scale_y_continuous(name = "Population density log(cells/ml)", limits = c(13.2, 15.8)) +
  scale_linetype_manual(name = "Dilution rate",
                        labels = c("1.25RPM" = "0.15/day",
                                   "2.5RPM" = "0.3/day"),
                        values = c("1.25RPM" = "solid",
                                   "2.5RPM" = "31")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA)))


# Leak
df_leak <- df %>%
  filter(test == "D") %>%
  filter(chamber == "4" | chamber == "6" | chamber == "8" | chamber == "9" | chamber == "10" | chamber == "13")

df_leak_sum <- df_leak %>%
  group_by(chamber, exp_day) %>%
  summarise(
    mean_cell_count = mean(cell_count),
    se_cell_count = (sd(cell_count)/sqrt(n()))) %>%
  ungroup()

ggplot(df_leak %>% mutate(mean_cell_count = cell_count), aes(x = exp_day, y = log(mean_cell_count), colour = factor(chamber))) +
  geom_point() +
  geom_smooth(data = df_leak_sum, aes(ymax = log(mean_cell_count + se_cell_count), ymin = log(mean_cell_count - se_cell_count)), stat = "identity") +
  scale_x_continuous(limits = c(30, 45), labels = seq(30, 45, 1), breaks = seq(30, 45, 1), name = "Days since inoculation") +
  scale_y_continuous(name = "Population density log(cells/ml)", limits = c(13.7, 15.3)) +
  theme_bw() +
  geom_segment(aes(x = 33, xend = 34.5, y = 14.1, yend = 14.1), arrow = arrow(length = unit(0.5, "cm")), size = 0.5, colour = "grey50") +
  annotate("text", x = 32.4, y = 14.1, label = "Leak", size = 8, colour = "grey35") +
  guides(color = guide_legend(override.aes = list(fill = NA)))

# Contamination
df_contamination <- df %>%
  filter(test == "D") %>%
  filter(chamber == "4" | chamber == "6" | chamber == "8" | chamber == "9" | chamber == "10" | chamber == "13") %>%
  mutate(contamination = ifelse(chamber == "4" | chamber == "13", "N", "Y"))

df_contamination_sum <- df_contamination %>%
  group_by(contamination, exp_day) %>%
  summarise(
    mean_cell_count = mean(cell_count),
    se_cell_count = (sd(cell_count)/sqrt(n()))) %>%
  ungroup()

ggplot(df_contamination %>% mutate(mean_cell_count = cell_count), aes(x = exp_day, y = log(mean_cell_count), colour = contamination, linetype = contamination)) +
  geom_point() +
  geom_smooth(data = df_contamination_sum, aes(ymax = log(mean_cell_count + se_cell_count), ymin = log(mean_cell_count - se_cell_count)), stat = "identity") +
  scale_x_continuous(limits = c(40, 76), labels = seq(40, 76, 2), breaks = seq(40, 76, 2), name = "Days since inoculation") +
  scale_y_continuous(name = "Population density log(cells/ml)", limits = c(13.2, 15.8)) +
  scale_linetype_manual(name = "Contamination",
                        labels = c("Y" = "Yes",
                                   "N" = "No"),
                        values = c("Y" = "solid",
                                   "N" = "21")) +
  theme_bw() +
  geom_segment(aes(x = 54, xend = 54, y = 15.5, yend = 15.1), arrow = arrow(length = unit(0.5, "cm")), size = 0.5, colour = "grey50", show.legend = FALSE) +
  annotate("text", x = 54, y = 15.6, label = "Possible contamination event", size = 8, colour = "grey35") +
  geom_segment(aes(x = 62, xend = 62, y = 13.5, yend = 14.1), arrow = arrow(length = unit(0.5, "cm")), size = 0.5, colour = "grey50", show.legend = FALSE) +
  annotate("text", x = 62, y = 13.40, label = "Contamination detected", size = 8, colour = "grey35") +
  guides(color = guide_legend(override.aes = list(fill = NA)))
