# Odds ratio for 10 (low), else (10 > x < 90), and high (>90)
# figure 1 and 2 (old ones) into a unique panel
pacman::p_load(data.table, ggplot2, ggthemr, envalysis, tidyverse, broom, patchwork)
options(scipen = 999) # disable scientific notation

# MAKE VERY SIMILAR TO THE https://www.nature.com/articles/s41380-023-02293-8
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

fem_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(females, IID, PRS), by = "IID")

man_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(males, IID, PRS), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4 + gender, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS = new_PRS)

shapiro.test(fem_pcs$PRS)
new_PRS_fem <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = fem_pcs))
females <- cbind(select(females, -PRS), PRS = new_PRS_fem)

shapiro.test(man_pcs$PRS)
new_PRS_man <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = man_pcs))
males <- cbind(select(males, -PRS), PRS = new_PRS_man)

# plot object
database_long <-
  select(database, W0, W1, W2, PRS) %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")), # that way the intercept is the low
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

females_long <-
  select(database, gender, W0, W1, W2, PRS) %>%
  filter(gender == "Female") %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")),
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

males_long <-
  select(database, gender, W0, W1, W2, PRS) %>%
  filter(gender == "Male") %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")), # that way the intercept is the low
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

# Testing if the inclusion of wave into the model changes anything
# Type 1: 90th x 10th
full_model <- glm(diagnosis ~ type1 * wave, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type1, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: no significance

# Type 2: 10th x else
full_model <- glm(diagnosis ~ type2 * wave, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type2, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: with significance

# Type 3: 90th x else
full_model <- glm(diagnosis ~ type3 * wave, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type3, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: with significance

# We want to see wave-specific odds ratios of diagnosis (ORs) for each risk group
# so we do stratify this

################### ALL SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
  database_long %>%
  filter(type1 %in% c("high", "low")) %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type1, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "low", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# Low risk vs. else stratified by wave
type2_results <-
  database_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type2, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "low"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "low") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# High risk vs. else stratified by wave
type3_results <-
  database_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type3, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

fp1 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "all"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "all"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "all"))

################### FEMALE SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
  females_long %>%
  filter(type1 %in% c("high", "low")) %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type1, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "low", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# Low risk vs. else stratified by wave
type2_results <-
  females_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type2, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "low"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "low") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# High risk vs. else stratified by wave
type3_results <-
  females_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type3, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

fp2 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "females"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "females"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "females"))

################### MALES SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
  males_long %>%
  filter(type1 %in% c("high", "low")) %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type1, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "low", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# Low risk vs. else stratified by wave
type2_results <-
  males_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type2, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "low"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "low") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# High risk vs. else stratified by wave
type3_results <-
  males_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type3, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

fp3 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "males"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "males"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "males"))

wd <-
  rbind(fp1, fp2, fp3) %>%
  mutate(
    strata = factor(strata, levels = c("all", "males", "females")),
    comparison = factor(comparison, levels = c("90th vs. 10th", "Else vs. 10th", "90th vs. Else")),
    signif = ifelse(p.value < 0.05, "*", ""))

final <- mutate(wd,
    estimate = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high),
    p.value = format.pval(p.value, digits = 3, eps = .001)) %>%
  select(strata, comparison, wave, estimate, p.value) %>%
  rename(
    Sample = strata,
    Comparison = comparison,
    Wave = wave,
    `OR (95% CI)` = estimate,
    `P-value` = p.value) %>%
  flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::align(part = "all", align = "center") %>%
  flextable::theme_booktabs() %>%
  flextable::autofit()
	
flextable::save_as_docx(
  "Supplementary Table S15" = final,
  path = "sup_tab15.docx")
