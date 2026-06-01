# Odds ratio for each PRS decile per wave
pacman::p_load(dplyr, broom, flextable, gtsummary)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	select(IID, gender, W0, W1, W2)

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  select(IID, PC1, PC2, PC3, PC4) %>%
	inner_join(., select(data$proband_data, IID, PRS), by = "IID")

new_PRS <- residuals(glm(
	PRS ~ PC1 + PC2 + PC3 + PC4,
	family = "gaussian",
	data = all_pcs)) %>%
	as.data.frame() %>%
	cbind(data$proband_data$IID, .) %>%
	rename(IID = 1, PRS = 2)

# Family history
hist <- data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# Working data
wd <-
	plyr::join_all(
		list(database, new_PRS, hist),
		by = "IID", type = "inner") %>%
  mutate(decile = ntile(PRS, 10)) %>%
	select(!c(PRS, IID))

# Models
mod_w0 <-
	glm(W0 ~ factor(decile) + gender + any_hist, family = binomial, data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w1 <-
	glm(W1 ~ factor(decile) + gender + any_hist, family = binomial, data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w2 <-
	glm(W2 ~ factor(decile) + gender + any_hist, family = binomial, data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

final <-
	tbl_merge(
  		tbls = list(mod_w0, mod_w1, mod_w2),
  		tab_spanner = c("**W0**","**W1**","**W2**")) %>%
  	bold_labels() %>%
  	italicize_levels() %>%
  	as_flex_table()

save_as_docx(
  "Supplementary Table S7" = final,
  path = "sup_tab7.docx")
