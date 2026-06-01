# Multivariate regression of raw PRS
pacman::p_load(dplyr, broom, flextable, gtsummary)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	select(IID, site, gender, PRS, W0, W1, W2, age_W0, age_W1, age_W2)

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  select(IID, PC1, PC2, PC3, PC4)

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# Working data
wd <- plyr::join_all(
	list(database, all_pcs, hist),
	by = "IID", type = "inner") %>%
	select(-IID)

str(wd)

final <-
	glm(PRS ~ ., family = "gaussian", data = select(wd, !c(age_W0, age_W1, age_W2))) %>%
	tbl_regression(intercept = TRUE) %>%
	as_flex_table()

save_as_docx(
  "Supplementary Table S5" = final,
  path = "sup_tab5.docx")
