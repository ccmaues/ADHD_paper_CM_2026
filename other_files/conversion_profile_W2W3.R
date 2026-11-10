# Flag conversions between consecutive waves
pacman::p_load(dplyr, ggplot2, envalysis, tidyr)
data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")
database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
	mutate(
		w1_conversion = W0 == 0 & W1 == 1,
		w2_conversion = W1 == 0 & W2 == 1,
		w3_conversion = W2 == 0 & W3 == 1) %>%
	mutate(
		conversion_wave = factor(case_when(
			W0 != 0 ~ "Case at W0",
			w1_conversion ~ "W1",
			w2_conversion ~ "W2",
			w3_conversion ~ "W3",
			TRUE ~ "No conversion"),
			levels = c("Case at W0", "No conversion", "W1", "W2", "W3"))) %>%
	inner_join(., data$PCA_all_samples, by = "IID")

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

# adjusted PRS values (risk strata info)
new_PRS <-
	all_pcs %>%
	mutate(PRS = residuals(glm(
				 PRS ~ PC1 + PC2 + PC3 + PC4,
				 family = "gaussian",
				 data = .))) %>%
	select(IID, PRS)

for_plot <-
	cbind(select(database, IID, gender, conversion_wave), PRS = new_PRS$PRS) %>%
	mutate(n.risk = ntile(PRS, 10)) %>%
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	select(gender, n.risk, PC1, PC2, conversion_wave, PRS)

# Case counts by wave, with new conversions highlighted
case_data <- distinct(database, IID, .keep_all = TRUE)
case_waves <- c("W0", "W1", "W2", "W3")
total_cases <- c(
	sum(case_data$W0 == 1, na.rm = TRUE),
	sum(case_data$W1 == 1, na.rm = TRUE),
	sum(case_data$W2 == 1, na.rm = TRUE),
	sum(case_data$W3 == 1, na.rm = TRUE))
new_conversions <- c(
	0,
	sum(case_data$w1_conversion, na.rm = TRUE),
	sum(case_data$w2_conversion, na.rm = TRUE),
	sum(case_data$w3_conversion, na.rm = TRUE))

case_counts <- data.frame(
	wave = factor(rep(case_waves, times = 2), levels = case_waves),
	case_group = factor(
		rep(c("Other cases", "Conversions this wave"), each = length(case_waves)),
		levels = c("Other cases", "Conversions this wave")),
	n_cases = c(total_cases - new_conversions, new_conversions))

# Cases sample size and conversion per wave
ggplot(case_counts, aes(x = wave, y = n_cases, fill = case_group)) +
	geom_col() +
	labs(x = "Wave", y = "Number of cases", fill = "Case group") +
	scale_fill_manual(values = c(
										"Other cases" = "#c6c7c7",
										"Conversions this wave" = "#ff6600")) +
	theme_publish(base_size = 10) +
	theme(legend.position = "bottom")

# PCA of the cases
	ggplot(for_plot, aes(PC1, PC2, shape = gender, color = conversion_wave)) +
		geom_point(size = 4, alpha = 0.5) +
		labs(y = "PC2", x = "PC1", shape = "", color = "Conversion wave") +
		scale_color_manual(values = c(
											 "Case at W0" = "#7b3294",
											 "No conversion" = "#dfdfdf",
											 "W1" = "#ff6600",
											 "W2" = "#3366cc",
											 "W3" = "#109618")) +
		theme_publish(base_size = 10) +
		theme(legend.position = "bottom",
					legend.margin = margin(t = 0, r = 10, b = 0, l = 0),
					panel.grid = element_line(linewidth = 0.2))

# scatter plot testing the association of PRS against ADHD levels

# ggplot(for_plot, aes(y = PRS, color = conversion_wave)) +
# geom_histogram() +
# coord_flip()


# ADHD original labels
og_ADHD <-
	data$W3_QC %>%
	select(-p_height1) %>%
	inner_join(
		select(database, IID, site, gender, starts_with("age_W")),
		by = "IID") %>%
	inner_join(new_PRS, by = "IID")

# EDA: PRS and age by wave
age_vars <- names(og_ADHD)[grepl("^age_W", names(og_ADHD))]
eda_vars <- c("PRS", age_vars)

eda_summary <- og_ADHD %>%
	summarise(across(
		all_of(eda_vars),
		list(
			n = ~ sum(!is.na(.x)),
			missing = ~ sum(is.na(.x)),
			mean = ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE),
			sd = ~ if (sum(!is.na(.x)) < 2) NA_real_ else sd(.x, na.rm = TRUE),
			median = ~ if (all(is.na(.x))) NA_real_ else median(.x, na.rm = TRUE),
			IQR = ~ if (all(is.na(.x))) NA_real_ else IQR(.x, na.rm = TRUE),
			min = ~ if (all(is.na(.x))) NA_real_ else min(.x, na.rm = TRUE),
			max = ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)
		),
		.names = "{.col}_{.fn}"
	)) %>%
	pivot_longer(
		cols = everything(),
		names_to = c("variable", ".value"),
		names_pattern = "^(.*)_(n|missing|mean|sd|median|IQR|min|max)$"
	)
print(eda_summary)

# Tests: PRS and age associations with ADHD subtype diagnosis codes
dca_vars <- intersect(
	c("dcadhdc", "dcadhdi", "dcadhdh", "dcadhdo"),
	names(og_ADHD)
)

if (length(dca_vars) > 0) {
	dca_counts <- bind_rows(lapply(dca_vars, function(dca_var) {
		as.data.frame(table(
			DCA_code = og_ADHD[[dca_var]],
			Gender = og_ADHD$gender,
			useNA = "ifany"
		)) %>%
			mutate(DCA_type = dca_var, .before = 1)
	}))
	print(dca_counts)

	prs_test_rows <- list()
	age_test_rows <- list()

	for (dca_var in dca_vars) {
		prs_data <- og_ADHD[, c("PRS", "gender", "age_W3", dca_var)]
		names(prs_data) <- c("PRS", "gender", "age_W3", "DCA_code")
		prs_data <- prs_data[complete.cases(prs_data), ]
		prs_data$gender <- factor(prs_data$gender)
		prs_data$DCA_code <- factor(prs_data$DCA_code)

		if (nlevels(droplevels(prs_data$DCA_code)) > 1) {
			has_multiple_genders <- nlevels(droplevels(prs_data$gender)) > 1
			prs_base <- if (has_multiple_genders) {
				lm(PRS ~ gender + age_W3, data = prs_data)
			} else {
				lm(PRS ~ age_W3, data = prs_data)
			}
			prs_by_dca <- if (has_multiple_genders) {
				lm(PRS ~ DCA_code + gender + age_W3, data = prs_data)
			} else {
				lm(PRS ~ DCA_code + age_W3, data = prs_data)
			}
			association_test <- anova(prs_base, prs_by_dca)

			interaction_test <- if (has_multiple_genders) {
				anova(
					prs_by_dca,
					lm(PRS ~ DCA_code * gender + age_W3, data = prs_data)
				)
			} else {
				NULL
			}
			prs_test_rows[[dca_var]] <- data.frame(
				DCA_type = dca_var,
				n = nrow(prs_data),
				association_df = association_test$Df[2],
				association_F = association_test$F[2],
				association_p = association_test[["Pr(>F)"]][2],
				interaction_df = if (is.null(interaction_test)) NA_real_ else interaction_test$Df[2],
				interaction_F = if (is.null(interaction_test)) NA_real_ else interaction_test$F[2],
				interaction_p = if (is.null(interaction_test)) NA_real_ else interaction_test[["Pr(>F)"]][2]
			)
		} else {
			prs_test_rows[[dca_var]] <- data.frame(
				DCA_type = dca_var, n = nrow(prs_data),
				association_df = NA_real_, association_F = NA_real_, association_p = NA_real_,
				interaction_df = NA_real_, interaction_F = NA_real_, interaction_p = NA_real_
			)
		}

		for (age_var in age_vars) {
			age_data <- og_ADHD[, c(age_var, dca_var, "gender")]
			names(age_data) <- c("age", "DCA_code", "gender")
			age_data <- age_data[complete.cases(age_data), ]
			age_data$DCA_code <- factor(age_data$DCA_code)
			age_data$gender <- factor(age_data$gender)

			if (nlevels(droplevels(age_data$DCA_code)) > 1) {
				age_base <- if (nlevels(droplevels(age_data$gender)) > 1) {
					lm(age ~ gender, data = age_data)
				} else {
					lm(age ~ 1, data = age_data)
				}
				age_by_dca <- if (nlevels(droplevels(age_data$gender)) > 1) {
					lm(age ~ DCA_code + gender, data = age_data)
				} else {
					lm(age ~ DCA_code, data = age_data)
				}
				age_test <- anova(age_base, age_by_dca)
				age_test_rows[[paste(dca_var, age_var, sep = "_")]] <- data.frame(
					DCA_type = dca_var,
					age_wave = age_var,
					n = nrow(age_data),
					df = age_test$Df[2],
					F = age_test$F[2],
					p = age_test[["Pr(>F)"]][2]
				)
			} else {
				age_test_rows[[paste(dca_var, age_var, sep = "_")]] <- data.frame(
					DCA_type = dca_var, age_wave = age_var, n = nrow(age_data),
					df = NA_real_, F = NA_real_, p = NA_real_
				)
			}
		}
	}

	prs_dca_tests <- bind_rows(prs_test_rows) %>%
		mutate(
			association_p_adjusted = p.adjust(association_p, method = "BH"),
			interaction_p_adjusted = p.adjust(interaction_p, method = "BH")
		)
	age_dca_tests <- bind_rows(age_test_rows) %>%
		mutate(p_adjusted = p.adjust(p, method = "BH"))

	print(prs_dca_tests)
	print(age_dca_tests)

	# Mean and 95% CI helper for the association plots
	mean_ci <- function(x) {
		n <- sum(!is.na(x))
		mean_x <- mean(x, na.rm = TRUE)
		se_x <- sd(x, na.rm = TRUE) / sqrt(n)
		data.frame(y = mean_x, ymin = mean_x - 1.96 * se_x, ymax = mean_x + 1.96 * se_x)
	}

	# Raw PRS distributions by DCA category and gender
	dca_prs_long <- og_ADHD %>%
		select(IID, gender, PRS, all_of(dca_vars)) %>%
		pivot_longer(
			cols = all_of(dca_vars),
			names_to = "DCA_type",
			values_to = "DCA_code"
		) %>%
		mutate(
			DCA_code = factor(DCA_code),
			gender = factor(gender)
		)

	p_prs_dca <- ggplot(dca_prs_long, aes(x = DCA_code, y = PRS, color = gender)) +
		geom_point(
			position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
			alpha = 0.5,
			na.rm = TRUE
		) +
		stat_summary(
			fun.data = mean_ci,
			geom = "errorbar",
			position = position_dodge(width = 0.7),
			width = 0.15,
			na.rm = TRUE
		) +
		stat_summary(
			fun = mean,
			geom = "point",
			position = position_dodge(width = 0.7),
			shape = 18,
			size = 2.5,
			na.rm = TRUE
		) +
		facet_wrap(~ DCA_type, scales = "free_x") +
		labs(x = "DCA code", y = "PRS", color = "Gender", title = "PRS by DCA code and gender",
			subtitle = "Raw observations and means with 95% CIs; tests adjust for age and gender") +
		theme_publish(base_size = 10) +
		theme(legend.position = "bottom")
	print(p_prs_dca)

	dca_age_long <- og_ADHD %>%
		select(IID, all_of(age_vars), all_of(dca_vars)) %>%
		pivot_longer(
			cols = all_of(age_vars),
			names_to = "age_wave",
			values_to = "age"
		) %>%
		pivot_longer(
			cols = all_of(dca_vars),
			names_to = "DCA_type",
			values_to = "DCA_code"
		) %>%
		mutate(
			DCA_code = factor(DCA_code),
			age_wave = factor(age_wave, levels = age_vars)
		)

	p_age_dca <- ggplot(
		dca_age_long,
		aes(x = age_wave, y = age, color = DCA_code, group = DCA_code)
	) +
		stat_summary(fun = mean, geom = "line", linewidth = 0.7, na.rm = TRUE) +
		stat_summary(fun.data = mean_ci, geom = "errorbar", width = 0.12, na.rm = TRUE) +
		stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
		facet_wrap(~ DCA_type) +
		labs(x = "Age wave", y = "Mean age", color = "DCA code",
			title = "Age trajectory by W3 DCA code",
			subtitle = "Raw means and 95% CIs; age tests adjust for gender") +
		theme_publish(base_size = 10)
	print(p_age_dca)
}
