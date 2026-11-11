# Table with the description of people before
# the cut for survival analysis in each wave
# add the info to the figure 1 too (just like how marcos asked)

pacman::p_load(dplyr, tidyr, purrr, tibble)

if (!exists("external_dir")) {
  stop("Define external_dir first: folder containing Santoro_278BHRC_2026_04_23.rds")
}

output_dir <- getwd()
fix_idate <- function(df) {
  malformed <- character()

  for (v in names(df)) {
    x <- df[[v]]

    if (inherits(x, "IDate") && typeof(x) != "integer") {
      malformed <- c(malformed, v)
      df[[v]] <- as.Date(
        as.numeric(unclass(x)), origin = "1970-01-01")
    }
  }

  if (length(malformed) > 0) {
    message("Repaired malformed IDate columns: ",
            paste(malformed, collapse = ", "))
  }

  df
}

raw_bhrc <- readRDS(
  file.path(external_dir, "Santoro_278BHRC_2026_04_23.rds")) %>%
  fix_idate()

wave_names <- c(
  wave0_arm_1 = "W0",
  wave1_arm_1 = "W1",
  wave2_arm_1 = "W2",
  wave3_arm_1 = "W3")

wave_levels <- c("W0", "W1", "W2", "W3")

# -------------------------------------------------------------------------
# Original BHRC / Figure 1 reference counts (NOT reconstructed from RDS)
# -------------------------------------------------------------------------
figure1_reference <- tibble(
  wave = wave_levels,
  N_original_reference = c(2512L, 2010L, 1944L, NA_integer_),
  reference_source = c(
    "BHRC GitLab README: baseline protocol (also states cohort N=2511)",
    "Original Figure 1A; W1=2010 independently documented in BHRC literature",
    "Original Figure 1A; W2=1944 not verified in public GitLab README",
    "Current 2026 RDS: observed evaluation evidence"
  )
)

# The named instruments/variables were checked in the 2026-04-23 source.
required_vars <- c(
  "ext_genid", "redcap_event_name", "dawbaage", "dawbagender",
  "sdqhk", "dcanyhk", "apm_completed", "apm_date1",
  "p_psico_date1", "age", "gender", "site")

missing_vars <- setdiff(required_vars, names(raw_bhrc))
if (length(missing_vars) > 0) {
  stop("Missing raw variables: ", paste(missing_vars, collapse = ", "))
}

has_value <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

fmt_n_pct <- function(n, total) {
  if (is.na(total) || total == 0) return(NA_character_)
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

fmt_age <- function(median, q1, q3) {
  if (is.na(median)) return("NA")
  sprintf("%.1f [%.1f, %.1f]", median, q1, q3)
}

wave_db <-
  raw_bhrc %>%
  filter(redcap_event_name %in% names(wave_names), !is.na(ext_genid)) %>%
  mutate(wave = unname(wave_names[as.character(redcap_event_name)]))

# Check for duplicate participant/event records before deduplication.
wave_duplicates <-
  wave_db %>%
  count(wave, ext_genid) %>%
  filter(n > 1)

if (nrow(wave_duplicates) > 0) {
  stop("Duplicate ext_genid within wave. Resolve these records before summarizing.")
}

# -------------------------------------------------------------------------
# 1. Evaluation evidence: actual assessment data, not merely REDCap rows
# -------------------------------------------------------------------------
wave_db <-
  wave_db %>%
  mutate(
    has_DAWBA = has_value(dawbaage) |
      has_value(dawbagender) |
      has_value(sdqhk),
    has_APM = has_value(apm_date1) &
      !is.na(apm_completed) &
      as.character(apm_completed) == "1",
    has_psico_date = has_value(p_psico_date1),
    has_evaluation = has_DAWBA | has_APM | has_psico_date,

    # 0 = control, 1 = uncertain, 2 = ADHD case
    ADHD_raw = as.character(dcanyhk),

    # W3 has DAWBA age although the general 'age' field is missing.
    # Do NOT estimate ages from old waves here.
    age_observed = coalesce(as_numeric(dawbaage), as_numeric(age)),

    # Original REDCap categories: 2 = female, 1 = male
    sex_raw = case_when(
      as.character(gender) %in% c("2", "Female") ~ "Female",
      as.character(gender) %in% c("1", "Male") ~ "Male",
      TRUE ~ NA_character_),

    # Original REDCap categories: 1 = RS, 2 = SP
    site_raw = case_when(
      as.character(site) %in% c("1", "RS") ~ "RS",
      as.character(site) %in% c("2", "SP") ~ "SP",
      TRUE ~ NA_character_))

wave_db$FH <- NA_character_

if (exists("mi") && exists("family_hist")) {
  id_crosswalk <- mi %>%
    transmute(ext_genid = as.character(ext_genid), IID_match = as.character(IID)) %>%
    distinct()

  if (anyDuplicated(id_crosswalk$ext_genid)) {
    stop("The ID dictionary maps some ext_genid values to multiple IIDs")
  }

  fh_lookup <- family_hist %>%
    transmute(IID_match = as.character(IID), FH_linked = as.character(any_hist)) %>%
    distinct(IID_match, .keep_all = TRUE)

  wave_db <- wave_db %>%
    left_join(id_crosswalk, by = "ext_genid") %>%
    left_join(fh_lookup, by = "IID_match") %>%
    mutate(FH = FH_linked) %>%
    select(-FH_linked)
}

# -------------------------------------------------------------------------
# 2. Counts per wave, before phenotype reconstruction or PCA filtering
# -------------------------------------------------------------------------
evaluation_counts <-
  wave_db %>%
  group_by(wave) %>%
  summarise(
    N_event_records = n(),
    N_DAWBA = sum(has_DAWBA),
    N_APM_only = sum(!has_DAWBA & has_APM),
    N_psico_only = sum(!has_DAWBA & !has_APM & has_psico_date),
    N_evaluated = sum(has_evaluation),
    .groups = "drop") %>%
  arrange(match(wave, wave_levels))

# One participant can have uncertain/missing ADHD and still be evaluated.
s14_summary <-
  wave_db %>%
  filter(has_evaluation) %>%
  group_by(wave) %>%
  summarise(
    N_evaluated = n(),
    N_case = sum(ADHD_raw == "2", na.rm = TRUE),
    N_control = sum(ADHD_raw == "0", na.rm = TRUE),
    N_uncertain = sum(ADHD_raw == "1", na.rm = TRUE),
    N_ADHD_missing_or_other = sum(
      is.na(ADHD_raw) | !ADHD_raw %in% c("0", "1", "2")),
    N_age = sum(!is.na(age_observed)),
    age_median = if (all(is.na(age_observed))) NA_real_
      else median(age_observed, na.rm = TRUE),
    age_q1 = if (all(is.na(age_observed))) NA_real_
      else unname(quantile(age_observed, .25, na.rm = TRUE)),
    age_q3 = if (all(is.na(age_observed))) NA_real_
      else unname(quantile(age_observed, .75, na.rm = TRUE)),
    N_female = sum(sex_raw == "Female", na.rm = TRUE),
    N_male = sum(sex_raw == "Male", na.rm = TRUE),
    N_gender_missing = sum(is.na(sex_raw)),
    N_SP = sum(site_raw == "SP", na.rm = TRUE),
    N_RS = sum(site_raw == "RS", na.rm = TRUE),
    N_site_missing = sum(is.na(site_raw)),
    N_FH_available = sum(FH %in% c("0", "1")),
    N_FH_yes = sum(FH == "1", na.rm = TRUE),
    N_observed_ADHD_and_age = sum(
      ADHD_raw %in% c("0", "2") & !is.na(age_observed)),
    .groups = "drop") %>%
  arrange(match(wave, wave_levels))

if (any(with(s14_summary,
             N_case + N_control + N_uncertain +
             N_ADHD_missing_or_other != N_evaluated))) {
  stop("ADHD diagnostic categories do not sum to the evaluated total")
}

s14_formatted <-
  s14_summary %>%
  transmute(
    wave,
    `Evaluated participants (2026 RDS), N` = as.character(N_evaluated),
    `ADHD cases, n (%)` = mapply(fmt_n_pct, N_case, N_evaluated),
    `Controls, n (%)` = mapply(fmt_n_pct, N_control, N_evaluated),
    `Uncertain ADHD status, n` = as.character(N_uncertain),
    `Missing/other ADHD status, n` = as.character(N_ADHD_missing_or_other),
    `Age available, n` = as.character(N_age),
    `Age (years), median [Q1, Q3]` = mapply(
      fmt_age, age_median, age_q1, age_q3),
    `Female, n (%)` = mapply(fmt_n_pct, N_female, N_evaluated),
    `Male, n (%)` = mapply(fmt_n_pct, N_male, N_evaluated),
    `Gender missing, n` = as.character(N_gender_missing),
    `São Paulo, n (%)` = mapply(fmt_n_pct, N_SP, N_evaluated),
    `Rio Grande do Sul, n (%)` = mapply(fmt_n_pct, N_RS, N_evaluated),
    `Site missing, n` = as.character(N_site_missing),
    `Family history available, n` = as.character(N_FH_available),
    `Family history positive, n / available` = paste0(
      N_FH_yes, " / ", N_FH_available),
    `Observed binary ADHD + age available, n` = as.character(
      N_observed_ADHD_and_age)) %>%
  pivot_longer(-wave, names_to = "Characteristic", values_to = "value") %>%
  mutate(Characteristic = factor(Characteristic, levels = unique(Characteristic))) %>%
  pivot_wider(names_from = wave, values_from = value) %>%
  arrange(Characteristic) %>%
  mutate(Characteristic = as.character(Characteristic))

if (!any(s14_summary$N_FH_available > 0)) {
  s14_formatted <- s14_formatted %>%
    filter(!grepl("Family history", Characteristic))
}

figure1_counts <-
  evaluation_counts %>%
  left_join(figure1_reference, by = "wave") %>%
  mutate(
    N_figure1_reference = coalesce(N_original_reference, N_evaluated),
    difference_evaluated_vs_reference = N_evaluated - N_figure1_reference
  ) %>%
  arrange(match(wave, wave_levels))

figure1_row <-
  figure1_counts %>%
  transmute(
    wave,
    Characteristic = "Original cohort / Figure 1A reference, N",
    value = as.character(N_figure1_reference)) %>%
  pivot_wider(names_from = wave, values_from = value) %>%
  select(Characteristic, all_of(wave_levels))

s14_formatted <-
  bind_rows(
    figure1_row,
    s14_formatted %>% select(Characteristic, all_of(wave_levels)))

source_note <- paste0(
  "Original cohort / Figure 1A reference: W0=2,512 (BHRC GitLab README " ,
  "baseline protocol, which also calls the ongoing cohort N=2,511); " ,
  "W1=2,010 (prior Figure 1A, also reported in BHRC follow-up research); " ,
  "W2=1,944 (prior Figure 1A only; NOT verified in the public GitLab README); " ,
  "W3=", figure1_counts$N_figure1_reference[figure1_counts$wave == "W3"],
  " (2026 RDS assessment evidence; not a GitLab count). " ,
  "BHRC repository: https://gitlab.com/bhrc/bhrc . " ,
  "The remaining rows describe evaluated participants in the 2026 raw RDS " ,
  "before phenotype reconstruction and survival-analysis selection; " ,
  "the first and second rows have different inclusion definitions."
)

# -------------------------------------------------------------------------
# 3. Figure 1 reference source metadata is included in figure1_counts above
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 4. Optional Figure 1 filtering flow if original data-processing objects exist
# -------------------------------------------------------------------------
needed <- c("dob_info", "mi", "data", "imputed_pheno", "all_pca", "proband_data")

if (all(vapply(needed, exists, logical(1), inherits = TRUE))) {
  participants <- wave_db %>%
    filter(has_evaluation) %>%
    distinct(wave, ext_genid)

  stage_dob_id <- participants %>%
    semi_join(dob_info %>% distinct(ext_genid), by = "ext_genid")

  id_lookup <- mi %>%
    transmute(ext_genid = as.character(ext_genid), IID = as.character(IID)) %>%
    distinct()

  stage_IID <- stage_dob_id %>%
    inner_join(id_lookup, by = "ext_genid")

  stage_PRS <- stage_IID %>%
    semi_join(data %>% transmute(IID = as.character(IID)) %>% distinct(),
              by = "IID")

  stage_phenotype <- stage_PRS %>%
    semi_join(imputed_pheno %>% transmute(IID = as.character(IID)) %>% distinct(),
              by = "IID")

  stage_PCA <- stage_phenotype %>%
    semi_join(all_pca %>% transmute(IID = as.character(IID)) %>% distinct(),
              by = "IID")

  stage_final <- stage_PCA %>%
    semi_join(proband_data %>% transmute(IID = as.character(IID)) %>% distinct(),
              by = "IID")

  stages <- list(
    "01 Evaluated in wave" = participants,
    "02 DOB-source/ID overlap" = stage_dob_id,
    "03 Linked genetic IID" = stage_IID,
    "04 Included in PRS-merged database" = stage_PRS,
    "05 Retained by phenotype reconstruction" = stage_phenotype,
    "06 Available in PCA" = stage_PCA,
    "07 Included in final longitudinal cohort" = stage_final)

  figure1_filter_flow <-
    purrr::imap_dfr(stages, function(ids, stage) {
      tibble(wave = wave_levels) %>%
        left_join(ids %>% count(wave, name = "N_remaining"), by = "wave") %>%
        mutate(stage = stage,
               N_remaining = replace_na(N_remaining, 0L))
    }) %>%
    mutate(
      stage = factor(stage, levels = names(stages)),
      wave = factor(wave, levels = wave_levels)) %>%
    arrange(wave, stage) %>%
    group_by(wave) %>%
    mutate(N_excluded_since_previous = lag(N_remaining) - N_remaining) %>%
    ungroup() %>%
    mutate(wave = as.character(wave), stage = as.character(stage))

  write.csv(figure1_filter_flow,
    file.path(output_dir, "Figure1_filtering_by_wave.csv"),
    row.names = FALSE)

  cat("\n=== Figure 1: filtering by wave ===\n")
  print(figure1_filter_flow, n = Inf, width = Inf)
} else {
  cat("\nNote: filtering flow not exported. Run the original preprocessing script first.\n")
  cat("Missing objects: ",
      paste(needed[!vapply(needed, exists, logical(1), inherits = TRUE)],
            collapse = ", "), "\n")
}

# -------------------------------------------------------------------------
# 5. Export
# -------------------------------------------------------------------------
write.csv(s14_summary,
  file.path(output_dir, "sup_tab14_numeric.csv"), row.names = FALSE)
write.csv(s14_formatted,
  file.path(output_dir, "sup_tab14_formatted.csv"), row.names = FALSE)
write.csv(figure1_counts,
  file.path(output_dir, "Figure1_wave_counts.csv"), row.names = FALSE)
writeLines(source_note, con = file.path(output_dir, "sup_tab14_source_note.txt"))

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(s14_formatted) %>%
    flextable::theme_booktabs() %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::width(j = 1, width = 3.75) %>%
    flextable::width(j = 2:5, width = 1.45) %>%
    flextable::add_footer_lines(values = source_note) %>%
    flextable::fontsize(size = 7, part = "footer") %>%
    flextable::set_table_properties(layout = "fixed")

  # A landscape page prevents the characteristic names from being clipped.
  section <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(
      top = 0.5, bottom = 0.5, left = 0.5, right = 0.5))

  flextable::save_as_docx(
    "Supplementary Table S14. Observed participant characteristics before survival-analysis selection" = ft,
    path = file.path(output_dir, "sup_tab14.docx"),
    pr_section = section)
} else {
  cat("\nInstall flextable and officer for automatic Word export: ",
      "install.packages(c('flextable', 'officer'))\n")
}

cat("\n=== Supplementary Table S14 ===\n")
print(s14_formatted, n = Inf, width = Inf)
cat("\n=== Figure 1: observed evaluation counts ===\n")
print(figure1_counts, n = Inf, width = Inf)
cat("\nNOTE: Original BHRC/Figure 1A reference counts are labelled with sources.\n")
cat("W0 from GitLab baseline protocol; W1 corroborated by BHRC literature;\n")
cat("W2 = 1,944 is from the manuscript's previous Figure 1A, not verified in GitLab;\n")
cat("W3 is calculated from the 2026 RDS, not an official GitLab count.\n")