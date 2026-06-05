# Figure 5 - part 2 - table
pacman::p_load(dplyr, tidyr, DescTools, flextable)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
  inner_join(., data$PCA_all_samples, by = "IID")

# Family history
hist <- data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
# apply te correction AT the glm of R2
# shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))

database <-
  cbind(database, PRS_for_risk = new_PRS) %>%
  mutate(risk = ntile(PRS_for_risk, 5)) %>%
  inner_join(., hist, by = "IID")

# Pseudo-R2 longitudinally
# I could do a ANOVA for checking the R2 diffs
# per risk
ggthemr("fresh")

# add familty history
fp1 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + gender + any_hist + age_W0, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + gender + any_hist + age_W0, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + gender + any_hist + age_W0, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + gender + any_hist + age_W0, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + gender + any_hist + age_W0, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + gender + any_hist + age_W1, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + gender + any_hist + age_W1, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + gender + any_hist + age_W1, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + gender + any_hist + age_W1, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + gender + any_hist + age_W1, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + gender + any_hist + age_W2, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + gender + any_hist + age_W2, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + gender + any_hist + age_W2, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + gender + any_hist + age_W2, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + gender + any_hist + age_W2, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke")))

fp2 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke")))

fp3 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + any_hist + age_W0, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + any_hist + age_W1, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + any_hist + age_W2, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke")))

tab_all <-
  fp1 %>%
  mutate(cohort = "All")

tab_male <-
  fp2 %>%
  mutate(cohort = "Male")

tab_female <-
  fp3 %>%
  mutate(cohort = "Female")

final <-
  bind_rows(tab_all, tab_male, tab_female) %>%
  mutate(R2 = round(R2 * 100, 2)) %>%
  pivot_wider(names_from = risk, values_from = R2, names_prefix = "Q") %>%
  relocate(cohort, wave) %>%
  flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::align(part = "all", align = "center") %>%
  flextable::autofit()
	
flextable::save_as_docx(
  "Supplementary Table S16" = final,
  path = "sup_tab16.docx")
