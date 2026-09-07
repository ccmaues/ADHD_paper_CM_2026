# Odds ratio for each PRS decile per wave for males
pacman::p_load(dplyr, broom, flextable, gtsummary)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Male data
males <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	select(
		IID,
		gender,
		PRS,
		W0, W1, W2, W3,
		age_W0, age_W1, age_W2, age_W3) %>%
	filter(gender == "Male")

# PCA
sex_stratified_pcs_M <-
	data$PCA_by_sex %>%
	inner_join(
		select(males, IID, PRS),
		by = "IID")

# PRS correction
new_PRS <-
	sex_stratified_pcs_M %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = sex_stratified_pcs_M))) %>%
	select(
		IID,
		PRS = PRS_adj)

# Family history
hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(starts_with("parent_"), ~ . == 1),
			1,
			0)) %>%
	select(IID, any_hist)

# Working data
wd <-
	males %>%
	select(-PRS) %>%
	inner_join(new_PRS, by = "IID") %>%
	inner_join(hist, by = "IID") %>%
	mutate(
		decile = ntile(PRS, 10)) %>%
	select(-PRS, -IID)

# Models
mod_w0 <-
	glm(
		W0 ~ factor(decile) + any_hist + age_W0,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w1 <-
	glm(
		W1 ~ factor(decile) + any_hist + age_W1,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w2 <-
	glm(
		W2 ~ factor(decile) + any_hist + age_W2,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w3 <-
	glm(
		W3 ~ factor(decile) + any_hist + age_W3,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

# Final table
final <-
	tbl_merge(
		tbls = list(
			mod_w0,
			mod_w1,
			mod_w2,
			mod_w3),
		tab_spanner = c(
			"**W0**",
			"**W1**",
			"**W2**",
			"**W3**")) %>%
	bold_labels() %>%
	italicize_levels() %>%
	as_flex_table()

save_as_docx(
	"Supplementary Table S9" = final,
	path = "sup_tab9.docx")