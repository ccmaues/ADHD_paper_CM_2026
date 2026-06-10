# survival model
pacman::p_load(envalysis, purrr, tibble, dplyr)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

## Cox-regression
# NUll
cox1 <- coxph(Surv(time, status) ~ strata(percentile), data = filter(wd, gender == "Male"))
cox2 <- coxph(Surv(time, status) ~ strata(percentile) + site, data = filter(wd, gender == "Male"))
cox3 <- coxph(Surv(time, status) ~ strata(percentile) + any_hist, data = filter(wd, gender == "Male"))
cox4 <- coxph(Surv(time, status) ~ strata(percentile) + site + any_hist, data = filter(wd, gender == "Male"))

models <- list(
  "Stratified baseline" = cox1,
  "+ Site" = cox2,
  "+ Family history" = cox3,  
  "+ Site + Family history" = cox4)

comparison_table <-
  imap_dfr(models, ~{
    s <- summary(.x)
    tibble(
      Model = .y,
      Parameters = length(coef(.x)),
      logLik = as.numeric(logLik(.x)),
      AIC = AIC(.x),
      Concordance = s$concordance[1])}) %>%
  mutate(delta_AIC = AIC - min(AIC))

lrt_site <- anova(cox1, cox2, test = "LRT")
lrt_hist <- anova(cox1, cox3, test = "LRT")
lrt_site_hist <- anova(cox1, cox4, test = "LRT")

lrt_table <- tibble(
  Model = c(
    "Stratified baseline",
    "+ Site",
    "+ Family history",
    "+ Site + Family history"),
  LRT_ChiSq = c(
    NA,
    lrt_site$Chisq[2],
    lrt_hist$Chisq[2],
    lrt_site_hist$Chisq[2]),
  LRT_p = c(
    NA,
    lrt_site$`Pr(>|Chi|)`[2],
    lrt_hist$`Pr(>|Chi|)`[2],
    lrt_site_hist$`Pr(>|Chi|)`[2]))

final <-
  comparison_table %>%
  left_join(lrt_table, by = "Model") %>%
  relocate(Model) %>%
  mutate(across(c(logLik, AIC, delta_AIC, Concordance, LRT_ChiSq), ~ round(.x, 3)),
    LRT_p = scales::pvalue(LRT_p)) %>%
    flextable::flextable() %>%
    flextable::bold(part = "header") %>%
    flextable::align(part = "all", align = "center") %>%
    flextable::theme_booktabs() %>%
    flextable::autofit()

flextable::save_as_docx(
  "Supplementary Table S19_pt1_males" = final,
  path = "sup_tab19_part1.docx")

cox1 <- coxph(Surv(time, status) ~ strata(percentile), data = filter(wd, gender == "Female"))
cox2 <- coxph(Surv(time, status) ~ strata(percentile) + site, data = filter(wd, gender == "Female"))
cox3 <- coxph(Surv(time, status) ~ strata(percentile) + any_hist, data = filter(wd, gender == "Female"))
cox4 <- coxph(Surv(time, status) ~ strata(percentile) + site + any_hist, data = filter(wd, gender == "Female"))

models <- list(
  "Stratified baseline" = cox1,
  "+ Site" = cox2,
  "+ Family history" = cox3,  
  "+ Site + Family history" = cox4)

comparison_table <-
  imap_dfr(models, ~{
    s <- summary(.x)
    tibble(
      Model = .y,
      Parameters = length(coef(.x)),
      logLik = as.numeric(logLik(.x)),
      AIC = AIC(.x),
      Concordance = s$concordance[1])}) %>%
  mutate(delta_AIC = AIC - min(AIC))

lrt_site <- anova(cox1, cox2, test = "LRT")
lrt_hist <- anova(cox1, cox3, test = "LRT")
lrt_site_hist <- anova(cox1, cox4, test = "LRT")

lrt_table <- tibble(
  Model = c(
    "Stratified baseline",
    "+ Site",
    "+ Family history",
    "+ Site + Family history"),
  LRT_ChiSq = c(
    NA,
    lrt_site$Chisq[2],
    lrt_hist$Chisq[2],
    lrt_site_hist$Chisq[2]),
  LRT_p = c(
    NA,
    lrt_site$`Pr(>|Chi|)`[2],
    lrt_hist$`Pr(>|Chi|)`[2],
    lrt_site_hist$`Pr(>|Chi|)`[2]))

final <-
  comparison_table %>%
  left_join(lrt_table, by = "Model") %>%
  relocate(Model) %>%
  mutate(across(c(logLik, AIC, delta_AIC, Concordance, LRT_ChiSq), ~ round(.x, 3)),
    LRT_p = scales::pvalue(LRT_p)) %>%
    flextable::flextable() %>%
    flextable::bold(part = "header") %>%
    flextable::align(part = "all", align = "center") %>%
    flextable::theme_booktabs() %>%
    flextable::autofit()

flextable::save_as_docx(
  "Supplementary Table S19_pt2_females" = final,
  path = "sup_tab19_part2.docx")