pacman::p_load(dplyr, tidyr, geepack, flextable, officer)
options(scipen = 999)

data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")
#---------------------------------
# Dataset and covariates
#---------------------------------
hist <-
  data$family_history %>%
  select(IID, any_hist)

pcs <-
  data$PCA_all_samples %>%
  select(IID, PC1, PC2, PC3, PC4)

database <-
  data$proband_data %>%
  select(IID, gender, site, W0, W1, W2, W3, PRS) %>%
  mutate(
    across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
  inner_join(hist, by = "IID") %>%
  inner_join(pcs, by = "IID") %>%
  drop_na(IID, gender, site, any_hist,
          W0, W1, W2, W3, PRS,
          PC1, PC2, PC3, PC4) %>%
  mutate(
    gender = factor(gender, levels = c("Female", "Male")),
    site = factor(site))

prs_fit <-
  glm(
    PRS ~ PC1 + PC2 + PC3 + PC4,
    family = gaussian,
    data = database)

database <-
  database %>%
  mutate(
    PRS_adjusted = as.numeric(residuals(prs_fit)),
    PRS_z = as.numeric(scale(PRS_adjusted)),
    decile = ntile(PRS_adjusted, 10),
    type1 = case_when(
      decile == 1 ~ "low",
      decile == 10 ~ "high",
      TRUE ~ "else"),
    type1 = factor(type1, levels = c("low", "else", "high")),
    cluster_id = as.integer(factor(IID)))

#---------------------------------
# Repeated observations
#---------------------------------
database_long <-
  database %>%
  select(IID, cluster_id, gender, site, any_hist,
         W0, W1, W2, W3, PRS_z, type1) %>%
  pivot_longer(
    cols = c(W0, W1, W2, W3),
    names_to = "wave",
    values_to = "diagnosis") %>%
  mutate(
    wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
    wave_id = as.integer(wave),
    diagnosis = as.integer(diagnosis)) %>%
  arrange(cluster_id, wave_id)

cluster_check <-
  database_long %>%
  count(cluster_id, name = "n_waves")

cat("\nSex by PGS group counts (distinct participants):\n")
print(
  database %>%
    count(gender, type1, name = "N"),
  n = Inf)

case_counts_s27 <-
  database_long %>%
  group_by(gender, type1, wave) %>%
  summarise(
    N = n(),
    cases = sum(diagnosis == 1),
    controls = sum(diagnosis == 0),
    cumulative_pct = 100 * mean(diagnosis == 1),
    .groups = "drop")

test_threeway <- function(model, predictor, expected_df) {

  terms <- names(coef(model))
  parts <- strsplit(terms, ":", fixed = TRUE)

  idx <- which(vapply(parts, function(x) {
    length(x) == 3 &&
      any(startsWith(x, predictor)) &&
      any(startsWith(x, "gender")) &&
      any(startsWith(x, "wave"))},
      logical(1)))

  beta <- coef(model)[idx]
  V <- vcov(model)[idx, idx, drop = FALSE]

  chisq <- as.numeric(crossprod(beta, solve(V, beta)))

  tibble(
    Wald_ChiSq = chisq,
    DF = length(idx),
    P_value = pchisq(chisq, df = length(idx), lower.tail = FALSE))}

analyses <- c(
  "Three PGS groups",
  "Top vs bottom decile",
  "Continuous PGS (per SD)")

gee_models_s27 <- list()
results_s27 <- list()

for (analysis in analyses) {

  model_data <- if (analysis == "Top vs bottom decile") {
    database_long %>%
      filter(type1 %in% c("low", "high")) %>%
      mutate(type1 = droplevels(type1)) %>%
      arrange(cluster_id, wave_id)
        } else {
    database_long}

  predictor <- if (analysis == "Continuous PGS (per SD)") {
    "PRS_z"
        } else {
    "type1"}

  formula_gee <- reformulate(
    c("site", "any_hist", paste0(predictor, " * gender * wave")),
    response = "diagnosis")

  model <-
    geepack::geeglm(formula = formula_gee,
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

  expected_df <- if (analysis == "Three PGS groups") 6L else 3L
  test <- test_threeway(model, predictor, expected_df)

  results_s27[[analysis]] <-
    test %>%
    mutate(Comparison = analysis,
           N = n_participants,
           Person_waves = nrow(model_data),
           Alpha = as.numeric(model$geese$alpha),
           .before = 1)

  gee_models_s27[[analysis]] <- model}

tab_s27 <-
  bind_rows(results_s27) %>%
  mutate(Comparison = factor(Comparison, levels = analyses)) %>%
  arrange(Comparison) %>%
  mutate(BH_q = p.adjust(P_value, method = "BH"))

cat("\nSupplementary Table S27 (unrounded):\n")
print(tab_s27, n = Inf)

write.csv(
  tab_s27,
  file = "sup_tab27_raw.csv",
  row.names = FALSE)

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))}

tab_s27_display <-
  tab_s27 %>%
  transmute(
    Comparison = as.character(Comparison),
    N = N,
    `Person-waves` = Person_waves,
    `Wald Chi2` = sprintf("%.2f", Wald_ChiSq),
    df = DF,
    `p-value` = format_p(P_value),
    `BH q-value` = format_p(BH_q))

ft <-
  flextable::flextable(tab_s27_display) %>%
  flextable::theme_booktabs() %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = "Comparison", align = "left", part = "body") %>%
  flextable::add_footer_lines(values = paste(
    "Pooled binomial-logit GEEs with exchangeable working correlation and",
    "robust sandwich Wald tests. Each model includes site and psychiatric",
    "family history; all lower-order PGS, gender and wave terms are retained.",
    "PGS residualized on pooled PC1-PC4 and standardized across the full cohort.",
    "The low/high groups are the bottom/top pooled PGS deciles; the middle",
    "group contains deciles 2-9. The joint test targets the PGS x gender x",
    "wave interaction only. BH q-values are exploratory across three tests.",
    "Outcome is cumulative recorded ADHD diagnosis, not biological onset.")) %>%
  flextable::autofit()

section <-
  officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 0.6, bottom = 0.6,
                                     left = 0.55, right = 0.55))

flextable::save_as_docx(
  "Supplementary Table S27. PGS-by-gender-by-wave interaction" = ft,
  path = "sup_tab27.docx",
  pr_section = section)
