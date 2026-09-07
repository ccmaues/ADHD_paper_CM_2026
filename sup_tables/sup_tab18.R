# PRS interaction verification with Odds Ratio model
pacman::p_load(
	dplyr,
	flextable,
	tidyr,
	broom,
	purrr)

# Article dataset
data <- readRDS(
	"C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# --------------------------------------------------
# Family history
# --------------------------------------------------

hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(
				starts_with("parent_"),
				~ . == 1),
			1,
			0)) %>%
	select(
		IID,
		any_hist)

# --------------------------------------------------
# Diagnosis standardization
# --------------------------------------------------

database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			as.factor)) %>%
	inner_join(
		hist,
		by = "IID")

# --------------------------------------------------
# PCA
# --------------------------------------------------

all_pcs <-
	data$PCA_all_samples %>%
	inner_join(
		select(
			database,
			IID,
			PRS),
		by = "IID") %>%
	select(
		IID,
		PC1,
		PC2,
		PC3,
		PC4,
		PRS)

# PCA-adjusted PRS
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

# --------------------------------------------------
# Working data
# --------------------------------------------------

wd <-
	database %>%
	select(-PRS) %>%
	inner_join(
		new_PRS,
		by = "IID") %>%
	rename(
		diagnosis_W0 = W0,
		diagnosis_W1 = W1,
		diagnosis_W2 = W2,
		diagnosis_W3 = W3) %>%
	select(
		IID,
		gender,
		site,
		starts_with("diagnosis"),
		starts_with("age_"),
		any_hist,
		PRS)

# --------------------------------------------------
# PRS × gender interaction
# --------------------------------------------------

m_w0 <-
	glm(
		diagnosis_W0 ~
			PRS * gender +
			age_W0 +
			site +
			any_hist,
		data = wd,
		family = "binomial")

m_w1 <-
	glm(
		diagnosis_W1 ~
			PRS * gender +
			age_W1 +
			site +
			any_hist,
		data = wd,
		family = "binomial")

m_w2 <-
	glm(
		diagnosis_W2 ~
			PRS * gender +
			age_W2 +
			site +
			any_hist,
		data = wd,
		family = "binomial")

m_w3 <-
	glm(
		diagnosis_W3 ~
			PRS * gender +
			age_W3 +
			site +
			any_hist,
		data = wd,
		family = "binomial")

mods <- list(
	W0 = m_w0,
	W1 = m_w1,
	W2 = m_w2,
	W3 = m_w3)

# --------------------------------------------------
# Extract interaction terms
# --------------------------------------------------

interactions <-
	imap_dfr(
		mods,
		~ tidy(
			.x,
			exponentiate = TRUE,
			conf.int = TRUE) %>%
			filter(
				term == "PRS:genderMale") %>%
			mutate(
				wave = .y)) %>%
	select(
		wave,
		estimate,
		conf.low,
		conf.high,
		p.value) %>%
	mutate(
		estimate = sprintf(
			"%.2f",
			estimate),

		conf.low = sprintf(
			"%.2f",
			conf.low),

		conf.high = sprintf(
			"%.2f",
			conf.high),

		p.value = format.pval(
			p.value,
			digits = 3,
			eps = .001))

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	interactions %>%
	rename(
		Wave = wave,
		`Interaction OR` = estimate,
		`95% CI lower` = conf.low,
		`95% CI upper` = conf.high,
		`P-value` = p.value) %>%
	flextable() %>%
	bold(
		part = "header") %>%
	align(
		part = "all",
		align = "center") %>%
	autofit()

save_as_docx(
	"Supplementary Table S18" = final,
	path = "sup_tab18.docx")