# Proportional Hazard ratio for cox model for each sex (Schoenfield Test)
pacman::p_load(dplyr, tidyr, patchwork, survminer, survival, gtsummary)
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

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

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

# PRS correction
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
  inner_join(., select(hist, IID, any_hist), by = "IID") %>%
  select(IID, site, W0, W1, W2, any_hist, age_W0, age_W1, age_W2, percentile, gender, site)

## keep controls the same
without_entry <-
  filter(database, W2 == 0) %>%
  select(IID, age_W2, W2) %>%
  rename(time = 2, status = 3)
str(without_entry)

## with the first occurance
temp1 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, W0, W1, W2) %>%
  pivot_longer(
    cols = starts_with("W"),
    names_to = "wave",
    values_to = "diagnosis")
str(temp1)

## Age data
temp2 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, age_W0, age_W1, age_W2) %>%
  pivot_longer(
    cols = starts_with("age_W"),
    names_to = "wave",
    values_to = "age") %>%
  mutate(wave = gsub("age_", "", wave))
str(temp2)

with_entry <-
  inner_join(temp1, temp2, by = c("IID", "wave")) %>%
  filter(diagnosis == 1) %>% # any time diagnosis
  group_by(IID) %>%
  filter(age == min(age)) %>%
  ungroup() %>%
  select(-wave) %>%
  rename(status = 2, time = 3)
str(with_entry)

survival_data_fem <-
  rbind(with_entry, without_entry) %>%
  inner_join(., select(database, IID, any_hist, site, percentile, gender, site), by = "IID") %>%
  filter(gender == "Female") %>%
  select(-IID, -gender) %>%
  data.frame()

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = survival_data_fem)
fit <- survfit(cox)
ph_fem <- cox.zph(cox)

survival_data_man <-
  rbind(with_entry, without_entry) %>%
  inner_join(., select(database, IID, any_hist, site, percentile, gender, site), by = "IID") %>%
  filter(gender == "Male") %>%
  select(-IID, -gender) %>%
  data.frame()

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = survival_data_man)
fit <- survfit(cox)
ph_man <- cox.zph(cox)

tab_fem <-
  data.frame(
    Variable = rownames(ph_fem$table),
    Chisq_F = ph_fem$table[, "chisq"],
    P_F = ph_fem$table[, "p"])

tab_man <-
  data.frame(
    Variable = rownames(ph_man$table),
    Chisq_M = ph_man$table[, "chisq"],
    P_M = ph_man$table[, "p"])

final <-
  left_join(tab_fem, tab_man, by = "Variable") %>%
  mutate(
    Chisq_F = sprintf("%.2f", Chisq_F),
    P_F = ifelse(P_F < 0.001, "<0.001", sprintf("%.3f", P_F)),
    Chisq_M = sprintf("%.2f", Chisq_M),
    P_M = ifelse(P_M < 0.001, "<0.001", sprintf("%.3f", P_M))) %>%
  flextable::flextable() %>%
  flextable::autofit()

flextable::save_as_docx(
  "Supplementary Table S11" = final,
  path = "sup_tab11.docx")