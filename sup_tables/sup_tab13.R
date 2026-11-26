pacman::p_load(
	dplyr,
	tidyr,
	purrr,
	flextable,
	officer)

# --------------------------------------------------
# Database and cohort preparation
# --------------------------------------------------

data <- readRDS(
	"C:/Users/Cass/Documents/work/cass_23112026_ARTICLE.rds")

waves <- c("W0", "W1", "W2", "W3")
cohort_levels <- c("All", "Male", "Female")
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

# --------------------------------------------------
# Model comparison on matched complete cases
# --------------------------------------------------

calc_r2 <- function(df, wave, r, cohort) {
	covars <- c("any_hist", paste0("age_", wave))
	if (cohort == "All") covars <- c("gender", covars)

	sub <- df %>% filter(risk == r)
	temp <- sub %>% select(all_of(c(wave, "PRS", covars))) %>% drop_na()
	y <- temp[[wave]]

	out <- tibble(
		cohort = cohort, risk = r, wave = wave,
		N_model = nrow(temp), N_cases_model = sum(y == 1),
		N_controls_model = sum(y == 0),
		PRS_SD_cohort = sd(df$PRS), PRS_SD_stratum = sd(sub$PRS),
		PRS_SD = sd(temp$PRS), PRS_IQR = IQR(temp$PRS),
		PRS_mean = if (nrow(temp)) mean(temp$PRS) else NA_real_,
		PRS_SD_ratio = sd(temp$PRS) / sd(df$PRS),
		PRS_case_control_diff = if (length(unique(y)) == 2)
			mean(temp$PRS[y == 1]) - mean(temp$PRS[y == 0]) else NA_real_,
		R2_PGS = NA_real_, R2_base = NA_real_, R2_full = NA_real_,
		delta_R2 = NA_real_, LRT_p = NA_real_, converged = NA,
		boundary = NA, note = "OK")

	if (length(unique(y)) != 2) {
		out$note <- "No outcome variation"
		return(out)
	}

	# Constant covariates are uninformative and cannot be fit as factors.
	keep_covars <- covars[vapply(temp[covars], function(x)
		length(unique(x)) > 1, logical(1))]
	if (length(keep_covars) != length(covars))
		out$note <- "Constant covariate omitted; review model"
	if (!is.finite(sd(temp$PRS)) || sd(temp$PRS) == 0) {
		out$note <- "No PGS variation"
		return(out)
	}

	fit <- function(x) glm(
		reformulate(if (length(x)) x else "1", response = wave),
		family = binomial, data = temp, na.action = na.fail)

	tryCatch({
		null <- fit(character(0))
		pgs <- fit("PRS")
		base <- fit(keep_covars)
		full <- fit(c("PRS", keep_covars))
		# Explicit common null likelihood for Nagelkerke pseudo-R2.
		ll0 <- as.numeric(logLik(null))
		n <- nrow(temp)
		nr2 <- function(m) {
			ll <- as.numeric(logLik(m))
			100 * (-expm1(2 * (ll0 - ll) / n)) / (-expm1(2 * ll0 / n))
		}
		out$R2_PGS <- nr2(pgs)
		out$R2_base <- nr2(base)
		out$R2_full <- nr2(full)
		out$delta_R2 <- out$R2_full - out$R2_base
		if (full$rank - base$rank != 1L)
			stop("PGS does not add one estimable coefficient")
		out$LRT_p <- anova(base, full, test = "LRT")[2, "Pr(>Chi)"]
		out$converged <- all(vapply(list(null, pgs, base, full),
									function(m) m$converged, logical(1)))
		out$boundary <- any(vapply(list(pgs, base, full),
								  function(m) m$boundary ||
									  any(fitted(m) < 1e-8 | fitted(m) > 1 - 1e-8),
								  logical(1)))
		out
	}, error = function(e) {
		out$note <- conditionMessage(e)
		out
	})
}

r2 <-
	imap_dfr(cohorts, function(df, cohort)
		map_dfr(waves, function(w)
			map_dfr(1:5, function(r) calc_r2(df, w, r, cohort))))

# --------------------------------------------------
# Analysis checks
# --------------------------------------------------

wd <- r2

problem_models <-
	wd %>%
	filter(
		note != "OK" |
		!coalesce(converged, FALSE) |
		coalesce(boundary, TRUE))

if (nrow(problem_models) > 0) {
	print(select(problem_models, cohort, wave, risk, note, converged, boundary))
	stop("Review the flagged model sets before exporting S13.")
}

needed <- c(
	"cohort", "risk", "wave", "N_model", "N_cases_model",
	"N_controls_model", "PRS_SD_ratio", "R2_PGS", "R2_base",
	"R2_full", "delta_R2", "LRT_p", "converged", "boundary", "note")

stopifnot(
	all(needed %in% names(wd)),
	nrow(wd) == 60L,
	!anyDuplicated(wd[c("cohort", "risk", "wave")]),
	setequal(wd$cohort, cohort_levels),
	setequal(wd$risk, 1:5),
	setequal(wd$wave, waves),
	all(wd$converged),
	!any(wd$boundary),
	all(wd$note == "OK"),
	all(wd$N_model == wd$N_cases_model + wd$N_controls_model),
	all(wd$N_cases_model > 0),
	all(wd$N_controls_model > 0),
	all(is.finite(wd$PRS_SD_ratio)),
	all(wd$PRS_SD_ratio > 0),
	all(is.finite(wd$R2_PGS)),
	all(is.finite(wd$R2_base)),
	all(is.finite(wd$R2_full)),
	all(is.finite(wd$delta_R2)),
	all(is.finite(wd$LRT_p)),
	all(wd$LRT_p >= 0 & wd$LRT_p <= 1),
	all(abs(wd$delta_R2 - (wd$R2_full - wd$R2_base)) < 1e-7),
	all(wd$delta_R2 >= -1e-7))

# BH adjustment across all 60 exploratory nested-model comparisons.
# Keep unrounded numerical values until export formatting.
wd <-
	wd %>%
	mutate(
		delta_R2 = pmax(delta_R2, 0),
		LRT_q = p.adjust(LRT_p, method = "BH"),
		SD_ratio_pct = 100 * PRS_SD_ratio,
		cohort = factor(
			cohort, levels = cohort_levels,
			labels = c("All samples", "Males", "Females")),
		wave = factor(wave, levels = waves),
		quintile = factor(risk, levels = 1:5, labels = paste0("Q", 1:5))) %>%
	arrange(cohort, wave, risk)

# --------------------------------------------------
# Final table
# --------------------------------------------------

s13 <-
	wd %>%
	transmute(
		cohort,
		wave,
		quintile,
		N = N_model,
		cases = N_cases_model,
		SD_ratio_pct,
		R2_PGS,
		R2_base,
		R2_full,
		delta_R2,
		LRT_p,
		LRT_q)

write.csv(s13, "sup_tab13.csv", row.names = FALSE)

note <- paste0(
	"N and cases refer to the matched complete-case analytical sample. ",
	"PGS SD ratio (%) = within-model PGS SD / cohort-wide PGS SD × 100. ",
	"PGS-only model: PC-adjusted PGS; clinical baseline: age at wave + psychiatric ",
	"family history (+ sex in the full cohort); full model: clinical baseline + PGS. ",
	"All pseudo-R² values are Nagelkerke estimates expressed as percentages. ",
	"ΔR² = full minus clinical baseline (percentage points, pp), not a unique ",
	"variance component. LRT compares the nested baseline and full models. ",
	"BH q-values are exploratory across 60 comparisons, including overlapping cohorts ",
	"and repeated waves. All estimates are apparent (in-sample).")

make_table <- function(df) {

	df %>%
		select(-cohort) %>%
		flextable() %>%
		set_header_labels(
			wave = "Wave",
			quintile = "PGS quintile",
			N = "N",
			cases = "Cases",
			SD_ratio_pct = "PGS SD ratio (%)",
			R2_PGS = "PGS only",
			R2_base = "Clinical baseline",
			R2_full = "Full",
			delta_R2 = "ΔR² (pp)",
			LRT_p = "LRT p",
			LRT_q = "BH q") %>%
		add_header_row(
			values = c("", "Nagelkerke pseudo-R² (%)", ""),
			colwidths = c(5, 3, 3)) %>%
		colformat_double(j = "SD_ratio_pct", digits = 2) %>%
		colformat_double(j = c("R2_PGS", "R2_base", "R2_full"), digits = 3) %>%
		colformat_double(j = "delta_R2", digits = 4) %>%
		colformat_double(j = c("LRT_p", "LRT_q"), digits = 4) %>%
		theme_booktabs() %>%
		bold(part = "header") %>%
		align(align = "center", part = "all") %>%
		fontsize(size = 9, part = "all") %>%
		padding(padding = 2.5, part = "all") %>%
		width(j = "wave", width = 0.58) %>%
		width(j = "quintile", width = 0.72) %>%
		width(j = c("N", "cases"), width = 0.52) %>%
		width(j = "SD_ratio_pct", width = 0.92) %>%
		width(j = c("R2_PGS", "R2_base", "R2_full"), width = 0.90) %>%
		width(j = "delta_R2", width = 0.83) %>%
		width(j = c("LRT_p", "LRT_q"), width = 0.76) %>%
		set_table_properties(layout = "fixed", opts_word = list(split = FALSE))
}

# --------------------------------------------------
# Word export - one cohort per landscape page
# --------------------------------------------------
section <-
	prop_section(
		page_size = page_size(orient = "landscape", width = 8.27, height = 11.69),
		page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.65, right = 0.65))

out <- read_docx()
out <- body_set_default_section(out, section)

for (i in seq_along(levels(s13$cohort))) {

	cohort_name <- levels(s13$cohort)[i]
	if (i > 1) out <- body_add_break(out)

	out <- body_add_fpar(
		out,
		fpar(ftext(
			paste0("Supplementary Table S13. Model fit across PGS quintiles (", cohort_name, ")"),
			fp_text(bold = TRUE, font.size = 11))))

	out <- body_add_flextable(out, make_table(filter(s13, cohort == cohort_name)))
	out <- body_add_fpar(out, fpar(ftext(note, fp_text(font.size = 8))))
}

print(out, target = "sup_tab13.docx")

# Optional console check for the main Figure 5 peak.
print(filter(s13, cohort == "Males", wave == "W0", quintile == "Q5"))
