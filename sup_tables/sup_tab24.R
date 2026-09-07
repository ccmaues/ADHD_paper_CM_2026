# Cumulative ADHD probability by age
# Supplementary Table S24

pacman::p_load(
	dplyr,
	tidyr,
	survival,
	flextable)

source(
	"C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

# --------------------------------------------------
# Working data
# --------------------------------------------------

wd <-
	survival_data %>%
	inner_join(
		select(
			hist,
			IID,
			any_hist),
		by = "IID") %>%
	select(-IID)

# --------------------------------------------------
# Cox models
# --------------------------------------------------

cox_all <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			gender +
			site,
		data = wd)

cox_fem <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist,
		data = filter(
			wd,
			gender == "Female"))

cox_man <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			site,
		data = filter(
			wd,
			gender == "Male"))

# --------------------------------------------------
# Survival curves
# --------------------------------------------------

fit_all <- survfit(cox_all)
fit_fem <- survfit(cox_fem)
fit_man <- survfit(cox_man)

# Ages reported in original table
ages <- c(
	5,
	10,
	12,
	15,
	20,
	25)

# --------------------------------------------------
# Extract cumulative event probability
# --------------------------------------------------

get_probability_table <- function(fit) {

	sm <-
		summary(
			fit,
			times = ages,
			extend = TRUE)

	tab <-
		tibble(
			strata = as.character(
				sm$strata),

			time = sm$time,

			probability =
				(1 - sm$surv) * 100) %>%

		mutate(
			# Handles labels such as
			# "percentile=10th" or
			# "strata(percentile)=10th"
			strata =
				sub(
					"^.*=",
					"",
					strata),

			strata = recode(
				strata,
				"10th" = "Bottom 10%",
				"else" = "Else",
				"90th" = "Top 10%"),

			estimate = sprintf(
				"%.2f%%",
				probability)) %>%

		select(
			strata,
			time,
			estimate) %>%

		pivot_wider(
			names_from = strata,
			values_from = estimate) %>%

		select(
			time,
			`Bottom 10%`,
			Else,
			`Top 10%`)

	return(tab)
}

# --------------------------------------------------
# All samples
# --------------------------------------------------

tab_all <-
	get_probability_table(
		fit_all)

final_all <-
	flextable(
		tab_all) %>%
	set_header_labels(
		time = "Age (yr)",
		`Bottom 10%` = "Bottom 10%",
		Else = "Else",
		`Top 10%` = "Top 10%") %>%
	bold(
		part = "header") %>%
	align(
		align = "center",
		part = "all") %>%
	theme_booktabs() %>%
	autofit()

save_as_docx(
	"Supplementary Table S24 - Cumulative ADHD probability by age (All samples)" =
		final_all,
	path = "sup_tab24_pt1.docx")

# --------------------------------------------------
# Females
# --------------------------------------------------

tab_fem <-
	get_probability_table(
		fit_fem)

final_fem <-
	flextable(
		tab_fem) %>%
	set_header_labels(
		time = "Age (yr)",
		`Bottom 10%` = "Bottom 10%",
		Else = "Else",
		`Top 10%` = "Top 10%") %>%
	bold(
		part = "header") %>%
	align(
		align = "center",
		part = "all") %>%
	theme_booktabs() %>%
	autofit()

save_as_docx(
	"Supplementary Table S24 - Cumulative ADHD probability by age (Females)" =
		final_fem,
	path = "sup_tab24_pt2.docx")

# --------------------------------------------------
# Males
# --------------------------------------------------

tab_man <-
	get_probability_table(
		fit_man)

final_man <-
	flextable(
		tab_man) %>%
	set_header_labels(
		time = "Age (yr)",
		`Bottom 10%` = "Bottom 10%",
		Else = "Else",
		`Top 10%` = "Top 10%") %>%
	bold(
		part = "header") %>%
	align(
		align = "center",
		part = "all") %>%
	theme_booktabs() %>%
	autofit()

save_as_docx(
	"Supplementary Table S24 - Cumulative ADHD probability by age (Males)" =
		final_man,
	path = "sup_tab24_pt3.docx")