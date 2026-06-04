# descriptive table
pacman::p_load(dplyr, officer, tibble, flextable)

# Article dataset
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Diagnosis standardtization
database <-
	data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	mutate(across(c(W0, W1, W2), as.factor))

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

wd <-
  inner_join(database, hist, by = "IID") %>%
  select(IID, site, W0, W1, W2, any_hist, age_W0, age_W1, age_W2, gender, site)

female_txt <- sprintf("%d (%.1f%%)", sum(wd$gender == "Female"), 100 * mean(wd$gender == "Female"))
rs_txt <- sprintf("%d (%.1f%%)", sum(wd$site == "RS"), 100 * mean(wd$site == "RS"))
fh_yes_txt <- sprintf("%d (%.1f%%)", sum(wd$any_hist == 1), 100 * mean(wd$any_hist == 1))
w0_txt <- sprintf("%d (%.1f%%)", sum(wd$W0 == "1"), 100 * mean(wd$W0 == "1"))
w1_txt <- sprintf("%d (%.1f%%)", sum(wd$W1 == "1"), 100 * mean(wd$W1 == "1"))
w2_txt <- sprintf("%d (%.1f%%)", sum(wd$W2 == "1"), 100 * mean(wd$W2 == "1"))
age0_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W0), quantile(wd$age_W0, .25), quantile(wd$age_W0, .75))
age1_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W1), quantile(wd$age_W1, .25), quantile(wd$age_W1, .75))
age2_txt <- sprintf("%.1f (%.1f–%.1f)", median(wd$age_W2), quantile(wd$age_W2, .25), quantile(wd$age_W2, .75))

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
  "Age (years), median (IQR)", "",
  "Wave 0", age0_txt,
  "Wave 1", age1_txt,
  "Wave 2", age2_txt)

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