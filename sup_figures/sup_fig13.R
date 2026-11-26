# Supplementary figure - Nagelkerke pseudo-R2 across PGS quintiles
# Standalone: identical preprocessing and model fitting to S13.
# Only the expanded RDS is required; no audit CSV or other script.
# All R2 columns are already percentages; delta_R2 is in percentage points.

pacman::p_load(dplyr, tidyr, purrr, ggplot2, patchwork, envalysis)

# --------------------------------------------------
# Database and cohort preparation
# Same eligibility and fixed quintiles as corrected S16/S17
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
# Plot settings
# --------------------------------------------------

wave_colors <- c(
	W0 = "#65ADC2",
	W1 = "#233B43",
	W2 = "#E84646",
	W3 = "#9B59B6")

common_theme <-
	theme_publish(base_size = 11) +
	theme(
		panel.grid.major.y = element_line(
			color = "#c1c1c1", linetype = "dashed", linewidth = 0.2),
		strip.background = element_blank(),
		strip.text = element_text(face = "bold", size = 10),
		legend.position = "bottom",
		legend.title = element_text(face = "bold"),
		plot.tag = element_text(face = "bold", size = 14),
		panel.spacing = grid::unit(0.8, "lines"))

# --------------------------------------------------
# A. PGS-only, clinical baseline and full models
# --------------------------------------------------

model_data <-
	wd %>%
	select(cohort, wave, quintile, R2_PGS, R2_base, R2_full) %>%
	pivot_longer(
		cols = c(R2_PGS, R2_base, R2_full),
		names_to = "model",
		values_to = "R2") %>%
	mutate(
		model = factor(
			model,
			levels = c("R2_PGS", "R2_base", "R2_full"),
			labels = c("PGS only", "Clinical baseline", "Full")))

p1 <-
	ggplot(model_data, aes(quintile, R2, color = wave, group = wave)) +
	geom_line(linewidth = 0.7, alpha = 0.65) +
	geom_point(size = 2.1) +
	facet_grid(model ~ cohort, scales = "fixed") +
	scale_color_manual(values = wave_colors, name = "Wave") +
	scale_y_continuous(
		breaks = seq(0, 100, by = 5),
		expand = expansion(mult = c(0.035, 0.07))) +
	labs(
		x = "PGS quintile",
		y = "Nagelkerke pseudo-R² (%)") +
	common_theme +
	theme(legend.position = "none")

# --------------------------------------------------
# B. PGS variability and incremental pseudo-R2
# --------------------------------------------------

# The x-axis is SD within the matched model sample / cohort SD (%).
# Repeated x positions are genuine: the PGS range is fixed across waves
# when the analytical sample is unchanged. Do not jitter these values.

p2 <-
	ggplot(wd, aes(SD_ratio_pct, delta_R2, color = wave, shape = quintile)) +
	geom_hline(yintercept = 0, linewidth = 0.3, color = "#828282") +
	geom_point(size = 2.5, alpha = 0.85, stroke = 0.8) +
	facet_wrap(~cohort, nrow = 1) +
	scale_color_manual(values = wave_colors, name = "Wave") +
	scale_shape_manual(values = c(16, 17, 15, 18, 1), name = "Quintile") +
	scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
	scale_y_continuous(expand = expansion(mult = c(0.045, 0.09))) +
	guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
	labs(
		x = "Within-model PGS SD / cohort PGS SD (%)",
		y = expression(Delta~"pseudo-R"^2~"(pp)")) +
	common_theme

# --------------------------------------------------
# Final supplementary figure
# --------------------------------------------------

final <-
	(p1 / p2) +
	plot_layout(heights = c(2.8, 1.2)) +
	plot_annotation(tag_levels = "A")

# The two panels must contain 180 and 60 points, respectively.
stopifnot(
	nrow(model_data) == 180L,
	!anyNA(model_data$wave),
	!anyNA(wd$wave),
	all(is.finite(model_data$R2)),
	all(table(wd$cohort, wd$wave, wd$quintile) == 1L),
	all(table(model_data$cohort, model_data$wave, model_data$quintile, model_data$model) == 1L))

print(wd %>% count(cohort, wave, name = "N_quintiles"))
write.csv(wd, "sup_fig_R2_quintiles_data.csv", row.names = FALSE)
print(final)

# Number-neutral filename: assign the final supplementary figure number
# when the manuscript supplement order is settled.

ggsave(
	filename = "sup_fig_R2_quintiles.png",
	plot = final,
	width = 24, height = 25,
	units = "cm", dpi = 600, bg = "white")

ggsave(
	filename = "sup_fig_R2_quintiles.pdf",
	plot = final,
	width = 24, height = 25,
	units = "cm", bg = "white")
