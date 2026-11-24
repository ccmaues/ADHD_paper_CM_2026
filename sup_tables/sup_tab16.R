# Figure 5 - part 2 - table
# Nagelkerke pseudo-R² by PRS quintile and wave

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
# Pseudo-R² function (full model; matched null model)
# --------------------------------------------------

calc_r2 <- function(df, wave, risk_group, cohort) {

	covars <- c("any_hist", paste0("age_", wave))
	if (cohort == "All") covars <- c("gender", covars)

	temp <-
		df %>%
		filter(risk == risk_group) %>%
		select(all_of(c(wave, "PRS", covars))) %>%
		drop_na()

	out <- tibble(
		cohort = cohort,
		wave = wave,
		risk = risk_group,
		R2 = NA_real_)

	label <- paste(cohort, wave, paste0("Q", risk_group))

	if (length(unique(temp[[wave]])) != 2 ||
		!is.finite(sd(temp$PRS)) || sd(temp$PRS) == 0) {
		warning(label, ": insufficient outcome or PRS variation.")
		return(out)
	}

	# Match the unified audit: omit constant covariates only if necessary.
	keep_covars <- covars[vapply(temp[covars], function(x)
		length(unique(x)) > 1, logical(1))]
	if (length(keep_covars) != length(covars)) {
		warning(label, ": constant covariate omitted; review model.")
	}

	null <- glm(
		reformulate("1", response = wave),
		family = binomial,
		data = temp,
		na.action = na.fail)

	model <- glm(
		reformulate(c("PRS", keep_covars), response = wave),
		family = binomial,
		data = temp,
		na.action = na.fail)

	if (!null$converged || !model$converged || model$boundary ||
		any(fitted(model) < 1e-8 | fitted(model) > 1 - 1e-8)) {
		warning(label, ": review convergence or boundary fitted values.")
	}

	# Nagelkerke pseudo-R² in percent, same formula as unified audit
	ll0 <- as.numeric(logLik(null))
	ll <- as.numeric(logLik(model))
	n <- nrow(temp)
	out$R2 <- 100 * (-expm1(2 * (ll0 - ll) / n)) /
		(-expm1(2 * ll0 / n))

	out
}

# --------------------------------------------------
# Analysis
# --------------------------------------------------

r2 <-
	imap_dfr(cohorts, function(df, cohort) {
		map_dfr(waves, function(wave) {
			map_dfr(1:5, function(r) {
				calc_r2(df, wave, r, cohort)
			})
		})
	})

# --------------------------------------------------
# Final table
# --------------------------------------------------

s16 <-
	r2 %>%
	mutate(
		wave = factor(wave, levels = waves),
		risk = factor(risk, levels = 1:5),
		R2 = round(R2, 2)) %>%
	pivot_wider(
		names_from = risk,
		values_from = R2,
		names_prefix = "Q") %>%
	relocate(cohort, wave)

final <-
	s16 %>%
	flextable() %>%
	bold(part = "header") %>%
	align(part = "all", align = "center") %>%
	autofit()

save_as_docx(
	"Supplementary Table S16" = final,
	path = "sup_tab16.docx")

write.csv(s16, "sup_tab16.csv", row.names = FALSE)
