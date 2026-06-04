# descriptive table for females and males
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

fem <- filter(wd, gender == "Female")
man <- filter(wd, gender == "Male")

fmt_pct <- function(x) {sprintf("%d (%.1f%%)", sum(x), 100 * mean(x))}
fmt_med <- function(x) {sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, .25), quantile(x, .75))}

tab_s1 <- tribble(
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

  "Age (years)", "", "",
  "Wave 0",
  fmt_med(fem$age_W0),
  fmt_med(man$age_W0),

  "Wave 1",
  fmt_med(fem$age_W1),
  fmt_med(man$age_W1),

  "Wave 2",
  fmt_med(fem$age_W2),
  fmt_med(man$age_W2))

section_rows <- which(
  tab_s1$Females == "" &
    tab_s1$Males == "")

ft_s1 <-
  tab_s1 %>%
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
  "Supplementary Table S2. Cohort characteristics per gender." = ft_s1,
  path = "sup_tab2.docx")
