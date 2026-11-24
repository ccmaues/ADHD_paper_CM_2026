# Figure 5 - part 1 - table
# ADHD prevalence by PRS quintile

pacman::p_load(
	dplyr,
	tidyr,
	purrr,
	flextable)

# --------------------------------------------------
# Database and common cohort preparation (S16 / S17)
# Same eligibility and quintiles as Fig5_unified_audit.R
# --------------------------------------------------

data <- readRDS(
	"C:/Users/Cass/Documents/work/cass_23112026_ARTICLE.rds")

waves <- c("W0", "W1", "W2", "W3")
pcs <- paste0("PC", 1:4)

# Proband data
proband <-
	data$proband_data %>%
	mutate(
		across(
			all_of(waves),
			~ ifelse(.x == 2, 1, .x)))

# Canonical family history from expanded processing
hist <-
	data$family_history %>%
	select(IID, any_hist)

# Verify unique identifiers and binary outcomes
check_id <- function(df) {
	stopifnot(!anyNA(df$IID), !anyDuplicated(df$IID))
}

invisible(lapply(
	list(proband, hist, data$PCA_all_samples, data$PCA_by_sex),
	check_id))

stopifnot(
	all(unlist(proband[waves]) %in% c(0, 1, NA)),
	all(proband$gender %in% c("Male", "Female")))

# PCA-adjusted PRS and fixed cohort-specific quintiles
make_cohort <- function(df, pca, cohort) {

	joined <-
		df %>%
		select(-any_of(c(pcs, "any_hist"))) %>%
		inner_join(
			select(pca, IID, all_of(pcs)),
			by = "IID") %>%
		inner_join(hist, by = "IID") %>%
		drop_na(PRS, all_of(pcs)) %>%
		arrange(IID)

	if (nrow(joined) != nrow(df)) {
		warning(cohort, ": retained ", nrow(joined), " of ", nrow(df),
			" participants; review eligibility.")
	}

	stopifnot(nrow(joined) >= 5)

	joined %>%
		mutate(
			PRS_raw = PRS,
			PRS = residuals(
				glm(
					PRS ~ PC1 + PC2 + PC3 + PC4,
					family = gaussian,
					data = joined,
					na.action = na.fail)),
			risk = ntile(PRS, 5))
}

# All samples: pooled PCA; sex subsets: sex-specific PCA
cohorts <- list(
	All = make_cohort(
		proband,
		data$PCA_all_samples,
		"All"),
	Male = make_cohort(
		filter(proband, gender == "Male"),
		data$PCA_by_sex,
		"Male"),
	Female = make_cohort(
		filter(proband, gender == "Female"),
		data$PCA_by_sex,
		"Female"))

# Check cohort and quintile sizes
print(imap_dfr(cohorts, ~ count(.x, risk) %>%
	mutate(cohort = .y, .before = 1)))

# --------------------------------------------------
# Prevalence by fixed risk quintile and wave
# --------------------------------------------------

prev <-
	imap_dfr(cohorts, ~ mutate(.x, cohort = .y)) %>%
	pivot_longer(
		all_of(waves),
		names_to = "wave",
		values_to = "diagnosis") %>%
	group_by(cohort, risk, wave) %>%
	summarise(
		N = n(),
		N_cases = sum(diagnosis == 1, na.rm = TRUE),
		N_controls = sum(diagnosis == 0, na.rm = TRUE),
		N_missing = sum(is.na(diagnosis)),
		N_observed = N_cases + N_controls,
		prevalence = if_else(
			N_observed > 0,
			100 * N_cases / N_observed,
			NA_real_),
		.groups = "drop")

stopifnot(all(prev$N == prev$N_observed + prev$N_missing))

# Check sample sizes and case counts
print(prev %>%
	group_by(cohort, wave) %>%
	summarise(
		N = sum(N),
		cases = sum(N_cases),
		.groups = "drop"))

# --------------------------------------------------
# Final table
# --------------------------------------------------

s17 <-
	prev %>%
	select(cohort, risk, wave, prevalence) %>%
	pivot_wider(
		names_from = wave,
		values_from = prevalence) %>%
	# Compute percentage-point change before rounding.
	mutate(
		delta = W3 - W0,
		across(c(all_of(waves), delta), ~ round(.x, 2))) %>%
	arrange(match(cohort, names(cohorts)), risk) %>%
	mutate(
		ntile = factor(
			risk,
			levels = 1:5,
			labels = c("1st", "2nd", "3rd", "4th", "5th"))) %>%
	select(cohort, ntile, all_of(waves), delta) %>%
	rename(`Δ W0–W3` = delta)

final <-
	s17 %>%
	flextable() %>%
	bold(part = "header") %>%
	align(part = "all", align = "center") %>%
	add_header_row(
		values = c("", "", "Prevalence (%)", ""),
		colwidths = c(1, 1, 4, 1)) %>%
	autofit()

save_as_docx(
	"Supplementary Table S17" = final,
	path = "sup_tab17.docx")

write.csv(s17, "sup_tab17.csv", row.names = FALSE)
