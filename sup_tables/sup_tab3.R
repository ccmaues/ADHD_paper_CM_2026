# suplementary descriptive table for females and males
pacman::p_load(dplyr, tidyr, flextable, gtsummary)

# then others (site, family history per gender)
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  inner_join(data$proband_data, data$family_history, by = "IID") %>%
	select(gender, site, PRS, starts_with("parent_")) %>%
	# filter(gender == "Female") %>%
	mutate(any_hist = if_else(
		if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(!starts_with("parent_"))


final <-
	database %>%
  tbl_summary(
    by = gender,
    include = c(site, PRS, any_hist),
    statistic = list(
      PRS ~ "{mean} ({sd})",
      any_hist ~ "{n} / {N} ({p}%)")) %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S3" = final,
  path = "sup_tab3.docx")
    