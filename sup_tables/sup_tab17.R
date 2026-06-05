# figure 5 - part 1 - table
pacman::p_load(dplyr, flextable)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")
database <-
  data$proband_data %>%
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS = new_PRS)

# Stratify by sex
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

calc_prev <- function(data, n, column_name, wave) {
  df <-
    select(data, all_of(column_name), wave) %>%
    mutate(quantile_n = ntile(!!sym(column_name), n)) %>%
    select(quantile_n, wave)
  p <-
    as.data.frame(cbind(1:n, table(df$quantile_n, df[[wave]]))) %>%
    rename("ntile" = 1, "Control" = 2, "Case" = 3) %>%
    mutate(
      prev = Case / (Control + Case),
      ntile = factor(ntile)) %>%
    select(-Control, -Case)
  return(p)}

## for all datas
p1 <- calc_prev(database, 5, "PRS", "W0") %>% rename(W0 = 2)
p2 <- calc_prev(database, 5, "PRS", "W1") %>% rename(W1 = 2)
p3 <- calc_prev(database, 5, "PRS", "W2") %>% rename(W2 = 2)

## Prevalence calculation (Female)
p1_fem <- calc_prev(females, 5, "PRS", "W0") %>% rename(W0 = 2)
p2_fem <- calc_prev(females, 5, "PRS", "W1") %>% rename(W1 = 2)
p3_fem <- calc_prev(females, 5, "PRS", "W2") %>% rename(W2 = 2)

## Prevalence calculation (Male)
p1_male <- calc_prev(males, 5, "PRS", "W0") %>% rename(W0 = 2)
p2_male <- calc_prev(males, 5, "PRS", "W1") %>% rename(W1 = 2)
p3_male <- calc_prev(males, 5, "PRS", "W2") %>% rename(W2 = 2)

# Delta plot preparation
all_plot_delta <-
	inner_join(p1, p3, by = "ntile") %>%
  mutate(deltaW0W2 = W2 - W0) %>%
  select(1, 4) %>%
  rename(delta = 2)

female_plot_delta <-
	inner_join(p1_fem, p3_fem, by = "ntile") %>%
  mutate(deltaW0W2 = W2 - W0) %>%
  select(1, 4) %>%
  rename(delta = 2)

male_plot_delta <-
	inner_join(p1_male, p3_male, by = "ntile") %>%
  mutate(deltaW0W2 = W2 - W0) %>%
  select(1, 4) %>%
  rename(delta = 2)

# Plot data preparation
for_plot_overall <-
	plyr::join_all(list(p1, p2, p3), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

for_plot_female <-
	plyr::join_all(list(p1_fem, p2_fem, p3_fem), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

for_plot_male <-
	plyr::join_all(list(p1_male, p2_male, p3_male), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

database_long <-
  select(database, W0, W1, W2, PRS) %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

females_long <-
  select(database, gender, W0, W1, W2, PRS) %>%
  filter(gender == "Female") %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

males_long <-
  select(database, gender, W0, W1, W2, PRS) %>%
  filter(gender == "Male") %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

final <-
  bind_rows(
    mutate(plyr::join_all(list(p1, p2, p3), by = "ntile"), cohort = "All"),
    mutate(plyr::join_all(list(p1_male, p2_male, p3_male), by = "ntile"), cohort = "Male"),
    mutate(plyr::join_all(list(p1_fem, p2_fem, p3_fem), by = "ntile"), cohort = "Female")) %>%
  mutate(across(starts_with("W"), ~ round(.x * 100, 2)), ntile = factor(ntile, levels = 1:5, labels = c("1st", "2nd", "3rd", "4th", "5th"))) %>%
  relocate(cohort, ntile) %>%
	mutate(delta = round((W2 - W0), 2)) %>%
	rename(`Δ W0–W2` = delta) %>%
  flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::align(part = "all", align = "center") %>%
	flextable::add_header_row(
  values = c("", "", "Prevalence (%)", ""),
  colwidths = c(1, 1, 3, 1))

flextable::save_as_docx(
  "Supplementary Table S17" = final,
  path = "sup_tab17.docx")
