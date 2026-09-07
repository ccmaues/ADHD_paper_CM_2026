# Odds ratio of ADHD diagnosis per wave predicted by the calculated PRS
# all samples
pacman::p_load(dplyr, broom, flextable, gtsummary)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Proband data
database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	select(
		IID,
		gender,
		W0, W1, W2, W3,
		age_W0, age_W1, age_W2, age_W3,
		site)

# PCA + raw PRS
all_pcs <-
	data$PCA_all_samples %>%
	select(
		IID,
		PC1, PC2, PC3, PC4) %>%
	inner_join(
		select(data$proband_data, IID, PRS),
		by = "IID")

# PRS correction
new_PRS <-
	all_pcs %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = all_pcs))) %>%
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
	plyr::join_all(
		list(
			database,
			new_PRS,
			hist),
		by = "IID",
		type = "inner") %>%
	select(-IID)

str(wd)

# Models
mod_w0 <-
	glm(
		W0 ~ gender + age_W0 + PRS + any_hist + site,
		family = binomial,
		data = wd)

mod_w1 <-
	glm(
		W1 ~ gender + age_W1 + PRS + any_hist + site,
		family = binomial,
		data = wd)

mod_w2 <-
	glm(
		W2 ~ gender + age_W2 + PRS + any_hist + site,
		family = binomial,
		data = wd)

mod_w3 <-
	glm(
		W3 ~ gender + age_W3 + PRS + any_hist + site,
		family = binomial,
		data = wd)

# Final table
final <-
	tbl_merge(
		tbls = list(
			tbl_regression(
				mod_w0,
				exponentiate = TRUE,
				intercept = TRUE),

			tbl_regression(
				mod_w1,
				exponentiate = TRUE,
				intercept = TRUE),

			tbl_regression(
				mod_w2,
				exponentiate = TRUE,
				intercept = TRUE),

			tbl_regression(
				mod_w3,
				exponentiate = TRUE,
				intercept = TRUE)),

		tab_spanner = c(
			"**W0**",
			"**W1**",
			"**W2**",
			"**W3**")) %>%
	as_flex_table()

save_as_docx(
	"Supplementary Table S6" = final,
	path = "sup_tab6.docx")