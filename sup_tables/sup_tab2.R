# descriptive table for females and males
pacman::p_load(dplyr, officer, tibble, flextable)

# Article dataset
data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Diagnosis standardization
database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	mutate(across(c(W0, W1, W2, W3), as.factor))

# Family history
hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(starts_with("parent_"), ~ . == 1),
			1,
			0)) %>%
	select(IID, any_hist)

wd <-
	inner_join(database, hist, by = "IID") %>%
	select(
		IID, site,
		W0, W1, W2, W3,
		any_hist,
		age_W0, age_W1, age_W2, age_W3,
		gender)

fem <- filter(wd, gender == "Female")
man <- filter(wd, gender == "Male")

fmt_pct <- function(x) {
	sprintf(
		"%d (%.1f%%)",
		sum(x),
		100 * mean(x))
}

fmt_med <- function(x) {
	sprintf(
		"%.1f [%.1f, %.1f]",
		median(x),
		quantile(x, .25),
		quantile(x, .75))
}

tab_s2 <-
	tribble(
		~Characteristic, ~Females, ~Males,

		"Study site", "", "",

		"São Paulo",
		fmt_pct(fem$site == "SP"),
		fmt_pct(man$site == "SP"),

		"Rio Grande do Sul",
		fmt_pct(fem$site == "RS"),
		fmt_pct(man$site == "RS"),

		"Family history", "", "",

		"Yes",
		fmt_pct(fem$any_hist == 1),
		fmt_pct(man$any_hist == 1),

		"No",
		fmt_pct(fem$any_hist == 0),
		fmt_pct(man$any_hist == 0),

		"ADHD diagnosis", "", "",

		"Wave 0",
		fmt_pct(fem$W0 == "1"),
		fmt_pct(man$W0 == "1"),

		"Wave 1",
		fmt_pct(fem$W1 == "1"),
		fmt_pct(man$W1 == "1"),

		"Wave 2",
		fmt_pct(fem$W2 == "1"),
		fmt_pct(man$W2 == "1"),

		"Wave 3",
		fmt_pct(fem$W3 == "1"),
		fmt_pct(man$W3 == "1"),

		"Age (years)", "", "",

		"Wave 0",
		fmt_med(fem$age_W0),
		fmt_med(man$age_W0),

		"Wave 1",
		fmt_med(fem$age_W1),
		fmt_med(man$age_W1),

		"Wave 2",
		fmt_med(fem$age_W2),
		fmt_med(man$age_W2),

		"Wave 3",
		fmt_med(fem$age_W3),
		fmt_med(man$age_W3))

section_rows <-
	which(
		tab_s2$Females == "" &
		tab_s2$Males == "")

ft_s2 <-
	tab_s2 %>%
	flextable() %>%
	theme_booktabs() %>%
	bold(part = "header") %>%
	bold(i = section_rows) %>%
	bg(i = section_rows, bg = "#F2F2F2") %>%
	hline(
		i = section_rows,
		border = officer::fp_border(width = 1)) %>%
	align(
		j = c("Females", "Males"),
		align = "center") %>%
	add_footer_lines(
		values = c(
			"Continuous variables are reported as median [Q1, Q3].",
			"Categorical variables are reported as n (%).")) %>%
	autofit()

save_as_docx(
	"Supplementary Table S2. Cohort characteristics per gender." = ft_s2,
	path = "sup_tab2.docx")