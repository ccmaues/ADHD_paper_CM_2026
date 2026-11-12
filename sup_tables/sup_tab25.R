# PRS incrementaal contribution
pacman::p_load(dplyr, tidyr, nsROC, PRROC, DescTools, purrr, flextable, officer)

data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")

# Family history
hist <-
	data$family_history %>%
	select(IID, any_hist)

# Proband data
database <-
  data$proband_data %>% # latest version
  inner_join(., hist, by = "IID") %>%
 	mutate(across(c(W0, W1, W2), ~ifelse(.x == 2, 1, 0))) # W3 is already as 1/0

datasets <- list(
  all = database %>%
    inner_join(data$PCA_all_samples, by = "IID"),
  females = database %>%
    inner_join(data$PCA_by_sex, by = "IID") %>%
    filter(gender == "Female"),
  males = database %>%
    inner_join(data$PCA_by_sex, by = "IID") %>%
    filter(gender == "Male"))

# Functions ----------------------------
# PRS correction
correct_prs <- function(df) {
  new_prs <- residuals(glm(
      PRS ~ PC1 + PC2 + PC3 + PC4,
      family = "gaussian",
      data = df))
  cbind(select(df, -PRS), PRS = new_prs)}

# evaluation
get_metrics <- function(wave, df, sex_adjust = TRUE) {

  age_var <- paste0("age_", wave)
  formula <- if (sex_adjust) {
    as.formula(paste0(wave, " ~ PRS + any_hist + gender + ", age_var))
      } else {
    as.formula(paste0(wave, " ~ PRS + any_hist + ", age_var))}

  model <- glm(formula, family = "binomial", data = df)
  observed <- model.response(model.frame(model))
  predicted <- as.numeric(predict(model, type = "response"))

  r2 <- as.numeric(PseudoR2(model, which = "Nagelkerke"))

  auroc <- as.numeric(PRROC::roc.curve(
    scores.class0 = predicted,
    weights.class0 = observed,
    curve = FALSE)$auc)

  aucpr <- as.numeric(PRROC::pr.curve(
    scores.class0 = predicted,
    weights.class0 = observed,
    curve = FALSE)$auc.integral)

  tibble(wave = wave,
         predictor = c("R2", "AUROC", "AUCPR"),
         value = c(r2, auroc, aucpr))}
# -----------------------

# PGS contribution
get_pgs_contribution <- function(wave, df, sex_adjust = TRUE) {

  age_var <- paste0("age_", wave)

  covars <- c("any_hist", if (sex_adjust) "gender", age_var)

  df_model <- df %>%
    select(all_of(c(wave, "PRS", covars))) %>%
    drop_na()

  # Model without PGS
  formula_base <-reformulate(covars, response = wave)

  model_base <- glm(
    formula_base,
    family = "binomial",
    data = df_model)

  # Model with PGS
  formula_full <- reformulate(c("PRS", covars), response = wave)

  model_full <- glm(
    formula_full,
    family = "binomial",
    data = df_model)

  # Performance metrics
  calc_metrics <- function(model) {

    observed <- model.response(model.frame(model))
    predicted <- predict(model, type = "response")

    c(R2 = as.numeric(PseudoR2(model, which = "Nagelkerke")),

      AUROC = as.numeric(
        roc.curve(
          scores.class0 = predicted,
          weights.class0 = observed,
          curve = FALSE)$auc),

      AUCPR = as.numeric(
        pr.curve(
          scores.class0 = predicted,
          weights.class0 = observed,
          curve = FALSE)$auc.integral))}

  base <- calc_metrics(model_base)
  full <- calc_metrics(model_full)

  lrt <- anova(model_base, model_full, test = "Chisq")

  # Final output
  tibble(wave = wave, n = nrow(df_model),
    R2_base = base["R2"],
    R2_full = full["R2"],
    delta_R2 = full["R2"] - base["R2"],

    AUROC_base = base["AUROC"],
    AUROC_full = full["AUROC"],
    delta_AUROC = full["AUROC"] - base["AUROC"],

    AUCPR_base = base["AUCPR"],
    AUCPR_full = full["AUCPR"],
    delta_AUCPR = full["AUCPR"] - base["AUCPR"],

    p_PGS = lrt[["Pr(>Chi)"]][2])}

tabs <- lapply(datasets, correct_prs)

waves <- c("W0", "W1", "W2", "W3")

# Analysis
pgs_contribution <-
  imap_dfr(tabs,
    ~ map_dfr(waves, get_pgs_contribution, df = .x, sex_adjust = (.y == "all")) %>%
    mutate(subset = .y, .before = 1))

format_metric <- function(x) sprintf("%.5f", x)
format_delta <- function(x) ifelse(abs(x) < 0.0001,
                                   sprintf("%.7f", x),
                                   sprintf("%.5f", x))
format_p <- function(x) ifelse(x < 0.001,
                               sprintf("%.6f", x),
                               ifelse(x < 0.01,
                                      sprintf("%.5f", x),
                                      sprintf("%.4f", x)))

table_s25 <-
  pgs_contribution %>%
  mutate(
    subset = factor(subset, levels = c("all", "females", "males")),
    wave = factor(wave, levels = c("W0", "W1", "W2", "W3"))) %>%
  arrange(subset, wave) %>%
  mutate(
    Cohort = case_when(
      subset == "all" ~ paste0("All (n = ", n, ")"),
      subset == "females" ~ paste0("Females (n = ", n, ")"),
      subset == "males" ~ paste0("Males (n = ", n, ")")),
    R2_base = format_metric(R2_base),
    R2_full = format_metric(R2_full),
    delta_R2 = format_delta(delta_R2),
    AUROC_base = format_metric(AUROC_base),
    AUROC_full = format_metric(AUROC_full),
    delta_AUROC = format_delta(delta_AUROC),
    AUCPR_base = format_metric(AUCPR_base),
    AUCPR_full = format_metric(AUCPR_full),
    delta_AUCPR = format_delta(delta_AUCPR),
    p_PGS = format_p(p_PGS)) %>%
  select(Cohort, wave, R2_base, R2_full, delta_R2,
         AUROC_base, AUROC_full, delta_AUROC,
         AUCPR_base, AUCPR_full, delta_AUCPR, p_PGS)

ft <-
  flextable(table_s25) %>%
  set_header_labels(
    Cohort = "Cohort (n)", wave = "Wave",
    R2_base = "Base", R2_full = "+ PGS", delta_R2 = "Delta",
    AUROC_base = "Base", AUROC_full = "+ PGS", delta_AUROC = "Delta",
    AUCPR_base = "Base", AUCPR_full = "+ PGS", delta_AUCPR = "Delta",
    p_PGS = "LRT p") %>%
  add_header_row(
    values = c("", "", "Nagelkerke pseudo-R2", "AUROC", "AUCPR", ""),
    colwidths = c(1, 1, 3, 3, 3, 1),
    top = TRUE) %>%
  theme_booktabs() %>%
  fontsize(size = 8, part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = "Cohort", align = "left", part = "body") %>%
  padding(padding = 3, part = "all") %>%
  set_caption(
    caption = "Supplementary Table S25. Incremental predictive performance after addition of ADHD polygenic scores") %>%
  add_footer_lines(values = paste0(
    "Base = age at wave + psychiatric family history (+ sex for the full cohort); ",
    "+ PGS = same covariates plus PC-adjusted PGS. ",
    "Delta = (+ PGS) - Base. LRT p = unadjusted likelihood-ratio-test p-value. ",
    "All performance metrics are apparent (in-sample), not validated.")) %>%
  autofit()

section <- officer::prop_section(
  page_size = officer::page_size(orient = "landscape"),
  page_margins = officer::page_mar(top = .6, bottom = .6, left = .55, right = .55))

flextable::save_as_docx(
  ft,
  path = "sup_tab25.docx",
  pr_section = section)

write.csv(pgs_contribution,
          "sup_tab25.csv",
          row.names = FALSE)
