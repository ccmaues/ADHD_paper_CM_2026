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

cat("\nCase counts by sex / group / wave:\n")
print(case_counts_s27, n = Inf)

if (any(case_counts_s27$cases == 0 | case_counts_s27$controls == 0)) {
  warning("Some risk-group x sex x wave cells have zero cases or controls; inspect model stability.")}

test_threeway <- function(model, predictor, expected_df) {

  terms <- names(coef(model))
  parts <- strsplit(terms, ":", fixed = TRUE)

  idx <- which(vapply(parts, function(x) {
    length(x) == 3 &&
      any(startsWith(x, predictor)) &&
      any(startsWith(x, "gender")) &&
      any(startsWith(x, "wave"))},
      logical(1)))

  if (length(idx) != expected_df) {
    stop("Expected ", expected_df, " three-way interaction terms; found ",
         length(idx), ". Check factor levels and model coefficients.")}

  beta <- coef(model)[idx]
  V <- vcov(model)[idx, idx, drop = FALSE]

  if (any(!is.finite(beta)) || any(!is.finite(V)) ||
      qr(V)$rank != length(idx)) {
    stop("Three-way interaction coefficients/covariance are invalid or singular.")}

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


model <- gee_models_s27[["Continuous PGS (per SD)"]]

beta <- coef(model)
V <- vcov(model)

n_by_gender <-
  database_long %>%
  distinct(IID, gender) %>%
  count(gender, name = "N") %>%
  mutate(gender = as.character(gender))

#---------------------------------
# Wave- and sex-specific linear contrasts
#---------------------------------
get_wave_or <- function(gender, wave) {

  L <- setNames(rep(0, length(beta)), names(beta))
  required <- "PRS_z"

  if (gender == "Male") {
    required <- c(required, "PRS_z:genderMale")
  }
  if (wave != "W0") {
    required <- c(required, paste0("PRS_z:wave", wave))
    if (gender == "Male") {
      required <- c(required, paste0("PRS_z:genderMale:wave", wave))
    }
  }

  missing_terms <- setdiff(required, names(beta))
  if (length(missing_terms)) {
    stop("Coefficient(s) not present in pooled GEE: ",
         paste(missing_terms, collapse = ", "))
  }
  L[required] <- 1

  log_OR <- sum(L * beta)
  SE <- sqrt(as.numeric(t(L) %*% V %*% L))
  if (!is.finite(SE) || SE <= 0) stop("Invalid SE for ", gender, " ", wave)

  tibble(
    gender = gender,
    wave = wave,
    log_OR = log_OR,
    SE_log_OR = SE,
    OR = exp(log_OR),
    CI_low = exp(log_OR - qnorm(0.975) * SE),
    CI_high = exp(log_OR + qnorm(0.975) * SE),
    p_value = 2 * pnorm(-abs(log_OR / SE)))
}

specs <- tidyr::expand_grid(
  gender = c("Female", "Male"),
  wave = c("W0", "W1", "W2", "W3"))

tab_s28 <-
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    get_wave_or(specs$gender[i], specs$wave[i])
  })) %>%
  left_join(n_by_gender, by = "gender") %>%
  mutate(
    # Exploratory adjustment of the EIGHT wave-specific post-hoc tests.
    BH_q = p.adjust(p_value, method = "BH")) %>%
  select(gender, wave, N, OR, CI_low, CI_high,
         p_value, BH_q, log_OR, SE_log_OR)

fmt_p <- function(p) ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))

tab_s28_display <-
  tab_s28 %>%
  transmute(
    Sex = gender,
    Wave = wave,
    N = format(N, big.mark = ",", scientific = FALSE),
    OR = sprintf("%.2f", OR),
    `95% CI` = sprintf("%.3f–%.3f", CI_low, CI_high),
    `p-value` = fmt_p(p_value),
    `BH q-value` = fmt_p(BH_q))

notes <- paste(
  "Odds ratios (ORs) per one standard deviation increase in the pooled",
  "PC1–PC4-residualized PGS. Derived from a single binomial-logit",
  "generalized estimating equation model with PGS × sex × wave,",
  "adjusted for study site and psychiatric family history.",
  "Exchangeable working correlation; cluster-robust sandwich 95% CIs.",
  "BH q-values are exploratory adjustments across these eight",
  "post-hoc wave-specific tests. Outcome: cumulative recorded ADHD diagnosis.")

ft <-
  flextable(tab_s28_display) %>%
  theme_booktabs() %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(j = c("Wave", "N", "OR", "95% CI", "p-value", "BH q-value"),
        align = "center", part = "all") %>%
  width(j = "Sex", width = 0.75) %>%
  width(j = "Wave", width = 0.62) %>%
  width(j = "N", width = 0.64) %>%
  width(j = "OR", width = 0.60) %>%
  width(j = "95% CI", width = 1.34) %>%
  width(j = "p-value", width = 0.87) %>%
  width(j = "BH q-value", width = 1.05) %>%
  set_table_properties(layout = "fixed") %>%
  add_footer_lines(notes) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer")

out_csv <- "sup_tab28_raw.csv"
out_docx <- "sup_tab28.docx"

write.csv(tab_s28, out_csv, row.names = FALSE)

report <- officer::read_docx() %>%
  officer::body_add_par("Supplementary Table S28", style = "Normal") %>%
  officer::body_add_par(
    "Wave-specific PGS associations with cumulative ADHD diagnosis, derived from pooled GEE",
    style = "Normal")
report <- flextable::body_add_flextable(report, ft)
print(report, target = out_docx)

cat("\nS28 generated:\n", out_csv, "\n", out_docx, "\n")
