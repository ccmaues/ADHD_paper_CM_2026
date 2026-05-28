# suplementary descriptive table for females
pacman::p_load(dplyr, data.table, tidyr, flextable, gtsummary)

# time-dependent variables (age, diagnosis per gender)
database <-
    readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$proband_data %>%
    filter(gender == "Female")

tab_diagnosis <-
  select(database, IID, W0, W1, W2) %>%
  pivot_longer(
    cols = c(W0, W1, W2),
    names_to = "wave",
    values_to = "diagnosis")

tab_age <-
  select(database, IID, age_W0, age_W1, age_W2) %>%
  pivot_longer(
    cols = c(age_W0, age_W1, age_W2),
    names_to = "wave",
    values_to = "age") %>%
  mutate(wave = gsub("age_", "", wave))

tab_wd <-
  inner_join(tab_age, tab_diagnosis, by = c("IID", "wave")) %>%
  mutate(
    diagnosis = ifelse(diagnosis == 0, "Control", "Case"),
    wave = case_when(
      wave == "W0" ~ "Wave 0",
      wave == "W1" ~ "Wave 1",
      wave == "W2" ~ "Wave 2"),
    across(c(diagnosis, wave), as.factor)) %>%
    select(-IID)

final <-
  tab_wd %>%
  tbl_summary(
    by = wave,
    include = c(age, diagnosis),
    statistic = list(
      age ~ "{mean} ({sd})",
      diagnosis ~ "{n} / {N} ({p}%)")) %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S1" = final,
  path = "sup_tab1.docx")
    