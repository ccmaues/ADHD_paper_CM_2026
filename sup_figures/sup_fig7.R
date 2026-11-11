pacman::p_load(
	dplyr, data.table, ggplot2, ggthemr,
	envalysis, tidyr, DescTools, patchwork)

data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

# PCA correction of PRS
pc_terms <- c("PC1", "PC2", "PC3", "PC4")

all_pcs <-
	data$PCA_all_samples %>%
	inner_join(
		select(
			database,
			IID, PRS, gender,
			age_W0, age_W1, age_W2, age_W3,
			W0, W1, W2, W3),
		by = "IID") %>%
	filter(if_all(all_of(c("PRS", pc_terms)), ~ !is.na(.x)))

all_pcs$PRS_adj <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = all_pcs))

all <- all_pcs

# PGS is retained in every model, including the fully adjusted model
calc_r2 <- function(database, wave, include_sex = TRUE) {
	age <- paste0("age_", wave)
	required <- c(wave, age, "PRS", "PRS_adj", pc_terms)

	if (include_sex) {
		required <- c(required, "gender")
	}

	wave_data <-
		database %>%
		filter(if_all(all_of(required), ~ !is.na(.x)))

	if (nrow(wave_data) == 0) {
		stop(paste("No complete observations for", wave))
	}

	# Confirm binary coding. The original W0-W3 recoding above is unchanged.
	if (!all(wave_data[[wave]] %in% c(0, 1))) {
		stop(paste("Outcome is not binary (0/1) in", wave))
	}

	terms <- list(
		"PRS raw"      = "PRS",
		"PRS adjusted" = "PRS_adj",
		"PRS + PC1"    = c("PRS", "PC1"),
		"PRS + PC2"    = c("PRS", "PC2"),
		"PRS + PC3"    = c("PRS", "PC3"),
		"PRS + PC4"    = c("PRS", "PC4"),
		"PRS + PCs"    = c("PRS", pc_terms),
		"PRS + age"    = c("PRS", age)
	)

	if (include_sex) {
		terms[["PRS + sex"]] <- c("PRS", "gender")
	}

	terms[["Full model"]] <- c(
		"PRS", pc_terms, age,
		if (include_sex) "gender")

	fit_model <- function(predictors) {
		glm(
			reformulate(predictors, response = wave),
			family = "binomial",
			data = wave_data)
	}

	models <- lapply(terms, fit_model)

	r2 <- vapply(
		models,
		function(model) {
			as.numeric(PseudoR2(model, which = "Nagelkerke"))
		},
		numeric(1))

	# With the PCs included, raw and residualized PRS reparameterize
	# the same full model; their fitted probabilities should agree.
	full_adj <- fit_model(c(
		"PRS_adj", pc_terms, age,
		if (include_sex) "gender"))

	if (!isTRUE(all.equal(
		as.numeric(fitted(full_adj)),
		as.numeric(fitted(models[["Full model"]])),
		tolerance = 1e-7))) {
		warning(paste("Full-model predictions differ after PRS residualization at", wave))
	}

	# All three models below contain PGS: changes are relative to PRS-only.
	comparison <- data.frame(
		wave = wave,
		N = nrow(wave_data),
		PRS_raw_R2 = unname(r2["PRS raw"]),
		PRS_adjusted_R2 = unname(r2["PRS adjusted"]),
		PRS_plus_PCs_R2 = unname(r2["PRS + PCs"]),
		full_R2 = unname(r2["Full model"]),
		delta_adjusted_vs_raw_R2 = unname(r2["PRS adjusted"] - r2["PRS raw"]),
		delta_PRS_plus_PCs_vs_raw_R2 = unname(r2["PRS + PCs"] - r2["PRS raw"]),
		delta_full_vs_raw_R2 = unname(r2["Full model"] - r2["PRS raw"])
	)

	result <- data.frame(
		wave = wave,
		var = factor(names(r2), levels = names(terms)),
		R2 = unname(r2),
		N = nrow(wave_data)
	)

	list(result = result, comparison = comparison)
}

# W0-W3
waves <- c("W0", "W1", "W2", "W3")
models <- lapply(waves, function(wave) calc_r2(all, wave))

results <- bind_rows(lapply(models, `[[`, "result"))
comparison <- bind_rows(lapply(models, `[[`, "comparison"))

print(comparison, row.names = FALSE)

write.csv(results, "fig7sup_prs_model_R2.csv", row.names = FALSE)
write.csv(comparison, "fig7sup_prs_comparison.csv", row.names = FALSE)

# Plot: each panel contains ONE wave only, ONE red PRS-only baseline
# The .env$wave syntax avoids matching the wave column to itself.
ggthemr("grape")

results <-
	results %>%
	mutate(label = sprintf("%.2f%%", R2 * 100))

max_r2 <- max(results$R2, na.rm = TRUE)
y_max <- max(0.01, ceiling(max_r2 * 1.30 / 0.01) * 0.01)

plot_wave <- function(wave) {
	tab <- filter(results, .data$wave == .env$wave)
	baseline <- tab$R2[as.character(tab$var) == "PRS raw"]

	stopifnot(nrow(tab) > 0, length(baseline) == 1)

	ggplot(tab, aes(x = var, y = R2, fill = var)) +
		geom_col(width = 0.8) +
		geom_text(aes(label = label), vjust = -0.4, size = 2.4) +
		geom_hline(
			yintercept = baseline,
			linetype = "dashed",
			color = "red") +
		scale_fill_grey(start = 0.2, end = 0.8) +
		scale_y_continuous(
			limits = c(0, y_max),
			breaks = pretty(c(0, y_max), n = 5),
			labels = function(x) paste0(round(100 * x, 1), "%"),
			expand = expansion(mult = c(0.01, 0.03))) +
		labs(
			title = paste0(wave, " (n = ", unique(tab$N), ")"),
			x = "",
			y = "Nagelkerke pseudo-R²") +
		theme_publish() +
		theme(
			legend.position = "none",
			axis.text.x = element_text(angle = 48, hjust = 1, size = 8),
			axis.text.y = element_text(size = 8),
			axis.title = element_text(size = 9),
			plot.title = element_text(size = 10),
			panel.grid.major.x = element_blank())
}

p1 <- plot_wave("W0")
p2 <- plot_wave("W1")
p3 <- plot_wave("W2")
p4 <- plot_wave("W3")

final <-
	p1 / p2 / p3 / p4 +
	plot_annotation(tag_levels = "A")

ggsave(
	"fig7_sup.png",
	final,
	device = "png",
	width = 250,
	height = 320,
	units = "mm",
	dpi = 300,
	bg = "white")
