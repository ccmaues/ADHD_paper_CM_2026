# Description of participants available before survival analysis
# per wave

pacman::p_load(
	dplyr,
	purrr,
	tibble,
	flextable,
	officer)

data <- readRDS(
	"C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Diagnosis standardization
database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x)))

# Family history
hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(starts_with("parent_"), ~ . == 1),
			1,
			0)) %>%
	select(
		IID,
		any_hist)

wd <-
	database %>%
	inner_join(
		hist,
		by = "IID")

# --------------------------------------------------
# Function
# --------------------------------------------------

describe_wave <- function(df, wave) {

	age_var <- paste0("age_", wave)

	temp <-
		df %>%
		filter(
			!is.na(.data[[wave]]),
			!is.na(.data[[age_var]]))

	tibble(
		Wave = wave,

		N = nrow(temp),

		ADHD =
			sprintf(
				"%d (%.1f%%)",
				sum(temp[[wave]] == 1),
				100 * mean(temp[[wave]] == 1)),

		Control =
			sprintf(
				"%d (%.1f%%)",
				sum(temp[[wave]] == 0),
				100 * mean(temp[[wave]] == 0)),

		Age =
			sprintf(
				"%.1f [%.1f, %.1f]",
				median(temp[[age_var]]),
				quantile(temp[[age_var]], .25),
				quantile(temp[[age_var]], .75)),

		Female =
			sprintf(
				"%d (%.1f%%)",
				sum(temp$gender == "Female"),
				100 * mean(temp$gender == "Female")),

		Male =
			sprintf(
				"%d (%.1f%%)",
				sum(temp$gender == "Male"),
				100 * mean(temp$gender == "Male")),

		SP =
			sprintf(
				"%d (%.1f%%)",
				sum(temp$site == "SP"),
				100 * mean(temp$site == "SP")),

		RS =
			sprintf(
				"%d (%.1f%%)",
				sum(temp$site == "RS"),
				100 * mean(temp$site == "RS")),

		Family_history =
			sprintf(
				"%d (%.1f%%)",
				sum(temp$any_hist == 1),
				100 * mean(temp$any_hist == 1)))
}

# --------------------------------------------------
# Table
# --------------------------------------------------

waves <- c("W0", "W1", "W2", "W3")

final_df <-
	map_dfr(
		waves,
		~ describe_wave(wd, .x))

final <-
	final_df %>%
	flextable() %>%
	theme_booktabs() %>%
	bold(part = "header") %>%
	align(
		j = 2:ncol(final_df),
		align = "center") %>%
	add_footer_lines(
		values = c(
			"Age is reported as median [Q1, Q3].",
			"Categorical variables are reported as n (%).",
			"Participants were included in each wave summary when both diagnosis status and age were available.")) %>%
	autofit()

save_as_docx(
	"Supplementary Table S12. Participant characteristics before survival analysis by wave." = final,
	path = "sup_tab12.docx")