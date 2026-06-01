# survival model
pacman::p_load(envalysis, purrr, tibble, dplyr)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")
## Cox-regression
# NUll
cox1 <- coxph(Surv(time, status) ~ strata(percentile), data = survival_data)
cox2 <- coxph(Surv(time, status) ~ strata(percentile) + site, data = survival_data)
cox3 <- coxph(Surv(time, status) ~ strata(percentile) + gender, data = survival_data)
cox4 <- coxph(Surv(time, status) ~ strata(percentile) + gender + site, data = survival_data)

models <- list(
  "Stratified baseline" = cox1,
  "+ Site" = cox2,
  "+ Gender" = cox3,
  "+ Site + Gender" = cox4)

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

lrt_site   <- anova(cox1, cox2, test = "LRT")
lrt_gender <- anova(cox1, cox3, test = "LRT")
lrt_both   <- anova(cox1, cox4, test = "LRT")

lrt_table <- tibble(
  Model = c(
    "Stratified baseline",
    "+ Site",
    "+ Gender",
    "+ Site + Gender"),
  LRT_ChiSq = c(
    NA,
    lrt_site$Chisq[2],
    lrt_gender$Chisq[2],
    lrt_both$Chisq[2]),
  LRT_p = c(
    NA,
    lrt_site$`Pr(>|Chi|)`[2],
    lrt_gender$`Pr(>|Chi|)`[2],
    lrt_both$`Pr(>|Chi|)`[2]))

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
  "Supplementary Table S12" = final,
  path = "sup_tab12.docx")
