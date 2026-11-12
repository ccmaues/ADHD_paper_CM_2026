# Supplementary Table S26 - Longitudinal PGS x wave interaction (GEE)
pacman::p_load(dplyr, tidyr, geepack, flextable, officer)

options(scipen = 999)

data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")
#---------------------------------
# Proband and family history data
#---------------------------------
hist <-
  data$family_history %>%
  select(IID, any_hist)

database <-
  data$proband_data %>%
  mutate(
    across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
  select(IID, gender, site, W0, W1, W2, W3, PRS) %>%
  inner_join(hist, by = "IID")

#---------------------------------
# PCA correction and risk groups
#---------------------------------
make_cohort <- function(cohort) {

  df <- switch(
    cohort,
    "All" = database,
    "Females" = filter(database, gender == "Female"),
    "Males" = filter(database, gender == "Male"))

  pcs <- if (cohort == "All") {
    data$PCA_all_samples
        } else {
    data$PCA_by_sex}

  df <-
    df %>%
    inner_join(select(pcs, IID, PC1, PC2, PC3, PC4), by = "IID") %>%
    drop_na(IID, gender, site, any_hist,
            W0, W1, W2, W3, PRS,
            PC1, PC2, PC3, PC4)
  prs_fit <-
    glm(PRS ~ PC1 + PC2 + PC3 + PC4,
        family = gaussian,
        data = df)

  df <-
    df %>%
    mutate(PRS = as.numeric(residuals(prs_fit)),
           PRS_z = as.numeric(scale(PRS)),
           decile = ntile(PRS, 10),
           type1 = case_when(decile == 1 ~ "low",
                             decile == 10 ~ "high",
                             TRUE ~ "else"),
      type1 = factor(type1, levels = c("low", "else", "high")))
  df}

#---------------------------------
# Long-format data and cluster IDs
#---------------------------------

make_long <- function(df) {

  df_long <-
    df %>%
    select(IID, gender, site, any_hist,
           W0, W1, W2, W3,
           PRS_z, decile, type1) %>%
    mutate(cluster_id = as.integer(factor(IID))) %>%
    pivot_longer(cols = c(W0, W1, W2, W3),
                 names_to = "wave",
                 values_to = "diagnosis") %>%
    mutate(wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
           wave_id = as.integer(wave),
           diagnosis = as.integer(diagnosis)) %>%
    arrange(cluster_id, wave_id)

  df_long}

#---------------------------------
# Robust joint Wald interaction test
#---------------------------------

test_interaction <- function(model) {
  # The models contain exactly one interaction term: PGS x wave.
  idx <- which(grepl(":wave", names(coef(model)), fixed = TRUE))
  stopifnot(length(idx) > 0)

  beta <- coef(model)[idx]
  V <- vcov(model)[idx, idx, drop = FALSE]

  if (qr(V)$rank != length(idx)) {
    stop("Interaction covariance matrix is singular; inspect this model.")}

  chisq <- as.numeric(t(beta) %*% solve(V, beta))

  tibble(ChiSq = chisq,
         DF = length(beta),
         P_value = pchisq(chisq,
                          df = length(beta),
                          lower.tail = FALSE))}

#---------------------------------
# Fit one GEE interaction model
#---------------------------------

fit_gee <- function(df_long, cohort, analysis) {
  # For high vs low: retain only the two extreme deciles.
  model_data <- if (analysis == "Top vs bottom decile") {
    df_long %>%
      filter(type1 %in% c("low", "high")) %>%
      mutate(type1 = droplevels(type1)) %>%
      arrange(cluster_id, wave_id)
        } else {
    df_long}

  predictor <- if (analysis == "Continuous PGS (per SD)") {
    "PRS_z"
        } else {
    "type1"}

  # Gender is included only for the pooled sample.
  covars <- c(
    if (cohort == "All") "gender",
    "any_hist", "site")

  formula_gee <- reformulate(
    c(covars, paste0(predictor, " * wave")),
    response = "diagnosis")

  model <-
    geepack::geeglm(
      formula_gee,
      id = cluster_id,
      waves = wave_id,
      data = model_data,
      family = binomial,
      corstr = "exchangeable",
      std.err = "san.se",
      control = geepack::geese.control(
        maxit = 40,
        trace = FALSE))

  n_participants <- n_distinct(model_data$IID)

  test <- test_interaction(model)

  result <-
    test %>%
    mutate(
      Cohort = cohort,
      Analysis = analysis,
      N = n_participants,
      N_observations = nrow(model_data),
      Alpha = as.numeric(model$geese$alpha),
      .before = 1)

  list(result = result, model = model)}

#---------------------------------
# Nine analyses: 3 cohorts x 3 tests
#---------------------------------
cohorts <- c("All", "Females", "Males")
analyses <- c(
  "Three PGS groups",
  "Top vs bottom decile",
  "Continuous PGS (per SD)")

gee_models <- list()
result_list <- list()
group_list <- list()

for (cohort in cohorts) {

  df_long <- make_long(make_cohort(cohort))

  # Verify that deciles have the expected group sizes.
  group_list[[cohort]] <-
    df_long %>%
    distinct(IID, type1) %>%
    count(type1, name = "n_participants") %>%
    mutate(Cohort = cohort, .before = 1)

  gee_models[[cohort]] <- list()

  for (analysis in analyses) {
    fit <- fit_gee(df_long, cohort, analysis)
    result_list[[paste(cohort, analysis, sep = "_")]] <- fit$result
    gee_models[[cohort]][[analysis]] <- fit$model}}

#---------------------------------
# Supplementary Table S26
#---------------------------------
tab_s26 <-
  bind_rows(result_list) %>%
  mutate(
    Cohort = factor(Cohort, levels = cohorts),
    Analysis = factor(Analysis, levels = analyses)) %>%
  arrange(Cohort, Analysis) %>%
  mutate(
    # Exploratory BH correction across all nine interaction tests.
    BH_q = p.adjust(P_value, method = "BH"))

cat("\nParticipant counts by risk group:\n")
print(bind_rows(group_list), n = Inf)

cat("\nS26 interaction results (unrounded):\n")
print(tab_s26, n = Inf)

write.csv(
  tab_s26,
  file = "sup_tab26_raw.csv",
  row.names = FALSE)

#---------------------------------
# Word table
#---------------------------------
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))}

tab_s26_display <-
  tab_s26 %>%
  transmute(
    Cohort = as.character(Cohort),
    Comparison = as.character(Analysis),
    N = N,
    `Person-waves` = N_observations,
    `Wald Chi2` = sprintf("%.2f", ChiSq),
    df = DF,
    `p-value` = format_p(P_value),
    `BH q-value` = format_p(BH_q))

ft <-
  flextable::flextable(tab_s26_display) %>%
  flextable::theme_booktabs() %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Cohort", "Comparison"),
                   align = "left", part = "body") %>%
  flextable::add_footer_lines(values = paste(
    "Binomial-logit GEE with exchangeable working correlation and robust",
    "sandwich Wald tests. All: adjusted for sex, study site, and psychiatric",
    "family history; sex-specific models: adjusted for study site and",
    "family history. Low/high = bottom/top within-cohort PGS decile;",
    "middle group = deciles 2-9. Continuous PGS standardized within cohort.",
    "PGS residualized on PC1-PC4 (pooled PCs for All, sex-specific PCs",
    "for sex-stratified models). BH q-values are exploratory across nine tests.",
    "Outcome is cumulative recorded ADHD diagnosis, not biological onset.")) %>%
  flextable::autofit()

section <- officer::prop_section(
  page_size = officer::page_size(orient = "landscape"),
  page_margins = officer::page_mar(
    top = 0.6, bottom = 0.6, left = 0.55, right = 0.55))

flextable::save_as_docx(
  "Supplementary Table S26. Longitudinal PGS-by-wave interactions" = ft,
  path = "sup_tab26.docx",
  pr_section = section)

# Optional model inspection:
summary(gee_models[["All"]][["Three PGS groups"]])
summary(gee_models[["All"]][["Top vs bottom decile"]])
summary(gee_models[["All"]][["Continuous PGS (per SD)"]])
