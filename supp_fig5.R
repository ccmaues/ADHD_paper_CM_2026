pacman::p_load(ggplot2, patchwork, dplyr, survival, survminer, tidyr, envalysis)

# Data
data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS), by = "IID") %>%
  select(-FID)

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <-
	cbind(select(database, -PRS), PRS = new_PRS) %>%
  mutate(
    risk = ntile(PRS, 100),
    percentile = case_when(
      risk >= 90 ~ "90th",
      risk <= 10 ~ "10th",
      TRUE ~ "else"),
    percentile = factor(
      percentile,
      levels = c("10th", "else", "90th"))) %>%
  select(IID, site, W0, W1, W2, age_W0, age_W1, age_W2, percentile, gender, site)

# keep controls the same
without_entry <-
  filter(database, W2 == 0) %>%
  select(IID, age_W2, W2) %>%
  rename(time = 2, status = 3)

# with the first occurance
temp1 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, W0, W1, W2) %>%
  pivot_longer(
    cols = starts_with("W"),
    names_to = "wave",
    values_to = "diagnosis")

# Age data
temp2 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, age_W0, age_W1, age_W2) %>%
  pivot_longer(
    cols = starts_with("age_W"),
    names_to = "wave",
    values_to = "age") %>%
  mutate(wave = gsub("age_", "", wave))

with_entry <-
  inner_join(temp1, temp2, by = c("IID", "wave")) %>%
  filter(diagnosis == 1) %>% # any time diagnosis
  group_by(IID) %>%
  filter(age == min(age)) %>%
  ungroup() %>%
  select(-wave) %>%
  rename(status = 2, time = 3)

# Make survival object
survival_data <-
  rbind(with_entry, without_entry) %>%
  inner_join(., select(database, IID, site, percentile, gender, site), by = "IID") %>%
  select(-IID) %>%
  data.frame()
str(survival_data)

# IMPORTANTE: O N não tá igual aos últimos... PQ? deveria ser 1558 e não 1553. REFAZER DPS.

## Cox model testing
cox <- coxph(Surv(time, status) ~ strata(percentile) + site + gender, data = survival_data)

# Panel A: Schoenfeld test to check the propotional harzard assumption
p1 <- ggcoxzph(cox.zph(cox), ggtheme = theme_publish())

# Panel B: check if any group is driving the model disproportionately
p2 <- ggcoxdiagnostics(cox, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_publish())
# Panel C: Likelihood Ratio Tests for Nested Models
model_null <- coxph(Surv(time, status) ~ strata(percentile), data = survival_data)
model_no_gender <- coxph(Surv(time, status) ~ strata(percentile) + site, data = survival_data)
model_no_site <- coxph(Surv(time, status) ~ strata(percentile) + gender, data = survival_data)
model_full <- coxph(Surv(time, status) ~ strata(percentile) + gender + site, data = survival_data)
anova(model_null, model_full, model_no_gender, model_no_site, test = "LRT")

# Note: the survival curves themselves serve as a
# verirfication for categorical variables  in
# cox model (overlap shows possible violation of
# proportional hazard assumption). This check is
# shown at the Main figure 2 of the letter.
# Martingale or deviance residuals are for continuous
# variables in the plot.

# Final plot
p1 / p2
ggpubr::as_ggplot(p1)
