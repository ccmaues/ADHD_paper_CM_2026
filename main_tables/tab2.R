# PRS table characteristics
pacman::p_load(dplyr, flextable, officer)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")
# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	inner_join(., hist, by = "IID")
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
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4 + gender, family = "gaussian", data = all_pcs))
new_PRS_fem <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = fem_pcs))
new_PRS_man <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = man_pcs))

wd <- list(
	all = cbind(select(database, -PRS), PRS = new_PRS),
	fem = cbind(select(females, -PRS), PRS = new_PRS_fem),
	man = cbind(select(males, -PRS), PRS = new_PRS_man)) %>%
	lapply(., mutate, percentile = ntile(PRS, 100)) %>%
	lapply(., mutate, risk_strata = case_when(
		percentile >= 90 ~ "High",
		percentile >= 10 ~ "Low",
		TRUE ~ "Intermediate"))

section_rows <- which(wd$Overall == "")

final <- wd$all %>%
  group_by(risk_strata) %>%
  summarise(
    N = n(),
    Female = sprintf("%d (%.1f%%)", sum(gender == "Female"), 100 * mean(gender == "Female")),
    SP = sprintf("%d (%.1f%%)", sum(site == "SP"), 100 * mean(site == "SP")),
    Family_history = sprintf("%d (%.1f%%)", sum(any_hist == 1), 100 * mean(any_hist == 1)),
    ADHD_W0 = sprintf("%d (%.1f%%)", sum(W0 == 1), 100 * mean(W0 == 1)),
    ADHD_W1 = sprintf("%d (%.1f%%)", sum(W1 == 1), 100 * mean(W1 == 1)),
    ADHD_W2 = sprintf("%d (%.1f%%)", sum(W2 == 1), 100 * mean(W2 == 1)),
  flextable() %>%
  bold(part = "header") %>%
  bold(i = section_rows, bold = TRUE) %>%
  bg(i = section_rows, bg = "#F2F2F2") %>%
  hline(
    i = section_rows,
    border = fp_border(width = 1)) %>%
  autofit() %>%
	add_footer_lines(
    values = c(
      "Categorical variables are reported as n (%).",
      "PRS strata were defined as Low (≤20th percentile), Intermediate (21st–79th percentile), and High (≥80th percentile)."))

save_as_docx(
  "Table S2. Cohort charfacteristiscs per risk strata." = final,
  path = "tab2.docx")

