# descriptive table
pacman::p_load(dplyr, officer, tibble, flextable)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
	mutate(across(c(W0, W1, W2, W3), as.factor))

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

wd <-
	inner_join(database, hist, by = "IID") %>%
	select(
		IID, site, W0, W1, W2, W3, any_hist,
		age_W0, age_W1, age_W2, age_W3, gender)

female_txt <- sprintf("%d (%.1f%%)", sum(wd$gender == "Female"), 100 * mean(wd$gender == "Female"))
sp_txt <- sprintf("%d (%.1f%%)", sum(wd$site == "SP"), 100 * mean(wd$site == "SP"))
fh_yes_txt <- sprintf("%d (%.1f%%)", sum(wd$any_hist == 1), 100 * mean(wd$any_hist == 1))
w0_txt <- sprintf("%d (%.1f%%)", sum(wd$W0 == "1"), 100 * mean(wd$W0 == "1"))
w1_txt <- sprintf("%d (%.1f%%)", sum(wd$W1 == "1"), 100 * mean(wd$W1 == "1"))
w2_txt <- sprintf("%d (%.1f%%)", sum(wd$W2 == "1"), 100 * mean(wd$W2 == "1"))
w3_txt <- sprintf("%d (%.1f%%)", sum(wd$W3 == "1"), 100 * mean(wd$W3 == "1"))
age0_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W0), quantile(wd$age_W0, .25), quantile(wd$age_W0, .75))
age1_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W1), quantile(wd$age_W1, .25), quantile(wd$age_W1, .75))
age2_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W2), quantile(wd$age_W2, .25), quantile(wd$age_W2, .75))
age3_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W3), quantile(wd$age_W3, .25), quantile(wd$age_W3, .75))


tab_s1 <- tribble(
	~Characteristic, ~Overall,
	"Demographics", "",
	"Female", female_txt,
	"São Paulo", sp_txt,
	"Family history", "",
	"ADHD diagnosis", "",
	"Wave 0", w0_txt,
	"Wave 1", w1_txt,
	"Wave 2", w2_txt,
	"Wave 3", w3_txt,
	"Age (years), median (IQR)", "",
	"Wave 0", age0_txt,
	"Wave 1", age1_txt,
	"Wave 2", age2_txt,
	"Wave 3", age3_txt)

section_rows <- which(tab_s1$Overall == "")

ft_s1 <-
  tab_s1 %>%
  flextable() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = section_rows, bold = TRUE) %>%
  bg(i = section_rows, bg = "#F2F2F2") %>%
  hline(
    i = section_rows,
    border = fp_border(width = 1)) %>%
  align(j = "Overall", align = "center") %>%
  autofit()

save_as_docx(
  "Table 1. Cohort characteristics." = ft_s1,
  path = "tab1.docx")