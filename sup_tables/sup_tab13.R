# phenotype imputation
# removal and stats
pacman::p_load(dplyr, tidyr, purrr, tibble, ggplot2)

# -------------------------------------------------------------------------
# 1. Paths and required objects
# -------------------------------------------------------------------------
if (!exists("external_dir")) {
  stop("The object 'external_dir' is missing. Run the main preprocessing script first.")
}

required_objects <- c(
  "mi", "dob_info", "bhrc_raw", "demographics", "longitudinal_wide",
  "prs", "data", "imputed_pheno", "new_age", "fam", "pca_files_IDs",
  "all_pca", "proband_data")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1),
                                          inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop("Missing preprocessing objects: ", paste(missing_objects, collapse = ", "))
}

output_dir <- file.path(getwd(), "S13_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source_file <- file.path(external_dir, "Santoro_278BHRC_2026_04_23.rds")
if (!file.exists(source_file)) {
  stop("Cannot find the source RDS: ", source_file)
}

waves <- c("W0", "W1", "W2", "W3")

source_bhrc <- readRDS(source_file)
raw_ext <- as.character(source_bhrc[["ext_genid"]])
raw_events <- as.character(source_bhrc[["redcap_event_name"]])
rm(source_bhrc)

raw_w0_ids <- unique(raw_ext[raw_events == "wave0_arm_1" &
                              !is.na(raw_ext) & nzchar(raw_ext)])


ids <- function(df, column = "IID") {
  if (!is.data.frame(df) || !column %in% names(df)) {
    stop("Expected a data.frame with column: ", column)
  }
  x <- as.character(df[[column]])
  unique(x[!is.na(x) & nzchar(x)])
}

assert_unique_ids <- function(df, object_name, column = "IID") {
  x <- as.character(df[[column]])
  x <- x[!is.na(x)]
  if (anyDuplicated(x)) {
    stop(object_name, " contains duplicate ", column,
         " values. Resolve these before interpreting the audit.")
  }
}

for (nm in c("data", "imputed_pheno", "new_age", "proband_data")) {
  assert_unique_ids(get(nm), nm)
}

id_map <- mi %>%
  transmute(ext_genid = as.character(ext_genid),
            IID = as.character(IID)) %>%
  filter(!is.na(ext_genid), !is.na(IID)) %>%
  distinct()
if (anyDuplicated(id_map$ext_genid) || anyDuplicated(id_map$IID)) {
  stop("mi is not a one-to-one ext_genid <-> IID map; examine its duplicates.")
}

# -------------------------------------------------------------------------
# 2. Sequential participant-level selection audit
# -------------------------------------------------------------------------
map_ext_to_iid <- function(x) {
  z <- id_map$IID[match(x, id_map$ext_genid)]
  unique(z[!is.na(z) & nzchar(z)])
}

raw_iid_matched <- map_ext_to_iid(raw_w0_ids)
dob_iid <- map_ext_to_iid(intersect(
  raw_w0_ids, ids(dob_info, "ext_genid")))

flow_id_sets <- list(
  `Raw W0 source records` = raw_w0_ids,
  `Linked to genetic ID dictionary` = intersect(raw_w0_ids,
                                                 ids(id_map, "ext_genid")),
  `Birth date and ID linkage` = intersect(raw_iid_matched, dob_iid),
  `Longitudinal phenotype + demographics` = intersect(
    dob_iid, intersect(ids(demographics), ids(longitudinal_wide))),
  `PRS available / database assembled` = ids(data),
  `Eligible for phenotype reconstruction` = ids(imputed_pheno),
  `Present in genetic FAM file` = ids(pca_files_IDs),
  `PCA eigenvectors available` = intersect(ids(pca_files_IDs), ids(all_pca)),
  `Final analysis sample` = ids(proband_data)
)

flow_check_sets <- list(
  flow_id_sets[[1]],
  flow_id_sets[[2]],
  dob_iid,
  intersect(dob_iid, flow_id_sets[[4]]),
  intersect(intersect(dob_iid, flow_id_sets[[4]]), ids(data)),
  intersect(intersect(intersect(dob_iid, flow_id_sets[[4]]), ids(data)),
            ids(imputed_pheno)),
  NULL,
  NULL,
  NULL
)

flow_check_sets[[7]] <- intersect(flow_check_sets[[6]], ids(pca_files_IDs))
flow_check_sets[[8]] <- intersect(flow_check_sets[[7]], ids(all_pca))
flow_check_sets[[9]] <- intersect(flow_check_sets[[8]], ids(proband_data))

for (i in 4:length(flow_check_sets)) {
  if (!all(flow_check_sets[[i]] %in% flow_check_sets[[i - 1]])) {
    stop("Non-nested sample sets at flow step ", i)
  }
}
if (length(flow_check_sets[[3]]) > length(flow_check_sets[[2]])) {
  stop("DOB/IID linkage unexpectedly increases participant count")
}

flow_labels <- names(flow_id_sets)
flow_n <- vapply(flow_check_sets, length, integer(1))
flow <- tibble(
  step_order = seq_along(flow_n),
  step = flow_labels,
  N_before = c(NA_integer_, head(flow_n, -1)),
  N_excluded = c(NA_integer_, head(flow_n, -1) - tail(flow_n, -1)),
  N_remaining = flow_n
)

checks <- tibble(
  object = c("data", "imputed_pheno", "pca_files_IDs", "proband_data"),
  N_object = c(length(ids(data)), length(ids(imputed_pheno)),
               length(ids(pca_files_IDs)), length(ids(proband_data))),
  N_in_sequential_audit = flow_n[c(5, 6, 7, 9)]
) %>% mutate(matches = N_object == N_in_sequential_audit)
if (any(!checks$matches)) {
  warning("Some existing objects contain IDs outside the traced baseline cohort; ",
          "inspect S13_object_checks.csv before using Figure 1C.")
}

genetic_qc_reference <- tibble(
  metric = c("Genotyped samples before QC (manuscript)",
             "Genotyped samples after QC (manuscript)",
             "Genotyped samples excluded by QC (manuscript)"),
  N = c(2511L, 2190L, 321L),
  source = "Current manuscript Methods, Section 2.1; not audited by IID here"
)

# -------------------------------------------------------------------------
# 3. Original longitudinal ADHD status and eligibility groups
# -------------------------------------------------------------------------
pre <- data %>%
  transmute(IID = as.character(IID),
            across(all_of(waves), ~ as.character(.x)))

status_mat <- as.matrix(pre[, waves])
first_idx <- apply(status_mat, 1, function(x) {
  p <- which(x == "2")
  if (length(p)) p[1] else NA_integer_
})
all_control <- apply(status_mat, 1, function(x) all(!is.na(x) & x == "0"))

eligibility <- tibble(
  IID = pre$IID,
  first_case_idx = first_idx,
  participant_group = ifelse(
    !is.na(first_idx), paste0("First ADHD case at ", waves[first_idx]),
    ifelse(all_control, "Control at all four waves",
           "Excluded: not diagnosed and not all-wave control")),
  selected_by_rule = !is.na(first_idx) | all_control,
  present_in_imputed_pheno = pre$IID %in% ids(imputed_pheno),
  present_in_final = pre$IID %in% ids(proband_data)
)

if (!identical(eligibility$selected_by_rule,
               eligibility$present_in_imputed_pheno)) {
  stop("Actual imputed_pheno participants do not match the stated phenotype rules. ",
       "Review the main preprocessing script before reporting S13.")
}

group_levels <- c(paste0("First ADHD case at ", waves),
                  "Control at all four waves",
                  "Excluded: not diagnosed and not all-wave control")

group_counts <- eligibility %>%
  mutate(participant_group = factor(participant_group, levels = group_levels)) %>%
  group_by(participant_group, .drop = FALSE) %>%
  summarise(
    N_before_reconstruction = n(),
    N_retained_after_reconstruction = sum(present_in_imputed_pheno),
    N_in_final_analytic_sample = sum(present_in_final),
    .groups = "drop"
  ) %>%
  mutate(participant_group = as.character(participant_group))

stopifnot(sum(group_counts$N_before_reconstruction) == nrow(pre))
stopifnot(sum(group_counts$N_retained_after_reconstruction) == nrow(imputed_pheno))
stopifnot(sum(group_counts$N_in_final_analytic_sample) == nrow(proband_data))

# -------------------------------------------------------------------------
# 4. Per-wave phenotype changes: compare EXACT SAME retained individuals
# -------------------------------------------------------------------------

post <- imputed_pheno %>%
  transmute(IID = as.character(IID),
            across(all_of(waves), ~ as.character(.x)))

case_first_post <- eligibility$first_case_idx[match(post$IID, eligibility$IID)]
for (j in seq_along(waves)) {
  expected <- ifelse(!is.na(case_first_post) & case_first_post <= j,
                     "2", "0")
  actual <- post[[waves[j]]]
  if (!identical(as.character(expected), as.character(actual))) {
    stop("Reconstructed statuses differ from the original rule at ", waves[j],
         ". Reinspect the main imputation code.")
  }
}

change_wide <- pre %>%
  semi_join(post, by = "IID") %>%
  left_join(post, by = "IID", suffix = c("_observed", "_reconstructed")) %>%
  left_join(eligibility %>% select(IID, first_case_idx), by = "IID")

pheno_long <- purrr::map_dfr(seq_along(waves), function(j) {
  w <- waves[j]
  observed <- change_wide[[paste0(w, "_observed")]]
  reconstructed <- change_wide[[paste0(w, "_reconstructed")]]
  tibble(
    IID = change_wide$IID,
    wave = w,
    wave_index = j,
    observed = observed,
    reconstructed = reconstructed,
    first_case_idx = change_wide$first_case_idx,
    original_class = case_when(
      observed == "2" ~ "Case (observed)",
      observed == "0" ~ "Control (observed)",
      observed == "1" ~ "Uncertain (observed)",
      TRUE ~ "Missing/other (observed)"),
    reconstructed_class = case_when(
      reconstructed == "2" ~ "Case (reconstructed)",
      reconstructed == "0" ~ "Control (reconstructed)",
      TRUE ~ "Unexpected status"),
    # A change from NA is counted as a reconstruction, not an observation.
    status_changed = is.na(observed) | observed != reconstructed,
    carried_forward_case = !is.na(first_case_idx) &
      wave_index > first_case_idx & reconstructed == "2" &
      (is.na(observed) | observed != "2")
  )
})

if (any(!pheno_long$reconstructed %in% c("0", "2"))) {
  stop("Post-reconstruction ADHD contains an unexpected value")
}

phenotype_by_wave <- pheno_long %>%
  group_by(wave) %>%
  summarise(
    N_retained = n_distinct(IID),
    N_observed_case = sum(observed == "2", na.rm = TRUE),
    N_observed_control = sum(observed == "0", na.rm = TRUE),
    N_observed_uncertain = sum(observed == "1", na.rm = TRUE),
    N_observed_missing_other = sum(is.na(observed) | !observed %in% c("0", "1", "2")),
    N_reconstructed_case = sum(reconstructed == "2"),
    N_reconstructed_control = sum(reconstructed == "0"),
    N_status_changed = sum(status_changed),
    N_case_carried_forward = sum(carried_forward_case),
    N_observed_control_to_case = sum(observed == "0" & reconstructed == "2", na.rm = TRUE),
    N_uncertain_missing_to_case = sum(
      (is.na(observed) | !observed %in% c("0", "2")) &
      reconstructed == "2"),
    N_uncertain_missing_to_control = sum(
      (is.na(observed) | !observed %in% c("0", "2")) &
      reconstructed == "0"),
    .groups = "drop"
  ) %>%
  arrange(match(wave, waves))

if (any(phenotype_by_wave$N_reconstructed_case +
        phenotype_by_wave$N_reconstructed_control != phenotype_by_wave$N_retained)) {
  stop("Post-reconstruction case/control counts do not equal retained N")
}
if (any(phenotype_by_wave$N_status_changed !=
        (phenotype_by_wave$N_observed_control_to_case +
         phenotype_by_wave$N_uncertain_missing_to_case +
         phenotype_by_wave$N_uncertain_missing_to_control))) {
  stop("Diagnosis transition counts do not reconcile with changed statuses")
}

# Optional diagnosis counts in the FINAL filtered cohort, which uses W3=1.
final_case_counts <- purrr::map_dfr(waves, function(w) {
  code <- as.character(proband_data[[w]])
  case_code <- if (w == "W3") "1" else "2"
  tibble(wave = w,
         N_final = length(code),
         N_cases_final = sum(code == case_code, na.rm = TRUE),
         N_controls_final = sum(code == "0", na.rm = TRUE))
})
phenotype_by_wave <- phenotype_by_wave %>%
  left_join(final_case_counts, by = "wave")

# More granular transition table for reviewer reproducibility.
transition_counts <- pheno_long %>%
  count(wave, original_class, reconstructed_class, name = "N") %>%
  arrange(match(wave, waves), original_class, reconstructed_class)

# -------------------------------------------------------------------------
# 5. Age reconstruction: date derived vs observed vs baseline offsets
# -------------------------------------------------------------------------
age_by_wave <- purrr::map_dfr(waves, function(w) {
  old_age <- suppressWarnings(as.numeric(as.character(data[[paste0("age_", w)]])))
  new_age_value <- suppressWarnings(as.numeric(
    as.character(new_age[[paste0("age_", w)]])))
  ids_data <- as.character(data$IID)
  ids_new <- as.character(new_age$IID)
  new_age_value <- new_age_value[match(ids_data, ids_new)]
  selected <- ids_data %in% ids(imputed_pheno)

  if (w == "W0") {
    source <- ifelse(!is.na(old_age), "Observed age", "Unavailable")
  } else {
    dates <- data[[paste0("p_psico_date1_", w)]]
    dob <- data[["birth_date"]]
    # Both variables should already be Date-class values in main preprocessing.
    age_calc <- as.numeric(difftime(as.Date(dates), as.Date(dob),
                                    units = "days")) / 365.25
    baseline <- suppressWarnings(as.numeric(as.character(data[["age_W0"]])))
    offset <- c(W1 = 3, W2 = 8, W3 = 12)[[w]]
    source <- case_when(
      !is.na(age_calc) ~ "Assessment-date derived",
      !is.na(old_age) ~ "Observed age",
      !is.na(baseline + offset) ~ "Baseline offset",
      TRUE ~ "Unavailable"
    )
  }

  tibble(IID = ids_data, wave = w, selected = selected,
         age_before = old_age, age_after = new_age_value,
         age_source = source)
})

age_summary <- age_by_wave %>%
  filter(selected) %>%
  group_by(wave) %>%
  summarise(
    N_participants = n(),
    N_age_originally_available = sum(!is.na(age_before)),
    N_age_after = sum(!is.na(age_after)),
    N_date_derived = sum(age_source == "Assessment-date derived"),
    N_original_age_used = sum(age_source == "Observed age"),
    N_baseline_offset_used = sum(age_source == "Baseline offset"),
    N_newly_filled_age = sum(is.na(age_before) & !is.na(age_after)),
    N_still_missing_age = sum(is.na(age_after)),
    .groups = "drop"
  ) %>%
  arrange(match(wave, waves))

if (any(age_summary$N_date_derived + age_summary$N_original_age_used +
        age_summary$N_baseline_offset_used + age_summary$N_still_missing_age !=
        age_summary$N_participants)) {
  warning("Age-source categories do not sum to the retained participant count")
}

# -------------------------------------------------------------------------
# 6. Print and export S13 + Figure 1C audit data
# -------------------------------------------------------------------------

cat("\n=== S13 A: Sequential participant filtering ===\n")
print(flow, n = Inf, width = Inf)
cat("\n=== S13 B: Phenotype selection groups ===\n")
print(group_counts, n = Inf, width = Inf)
cat("\n=== S13 C: Observed vs reconstructed ADHD by wave ===\n")
print(phenotype_by_wave, n = Inf, width = Inf)
cat("\n=== S13 D: Age data reconstruction ===\n")
print(age_summary, n = Inf, width = Inf)
cat("\n=== Source-to-final-object consistency ===\n")
print(checks, n = Inf, width = Inf)

stopifnot(nrow(proband_data) == flow$N_remaining[nrow(flow)])
if (nrow(proband_data) != 1463L) {
  warning("Final N is ", nrow(proband_data),
          ", not 1,463. Examine your current preprocessing inputs.")
}

write.csv(flow, file.path(output_dir, "S13_A_participant_filtering.csv"), row.names = FALSE)
write.csv(group_counts, file.path(output_dir, "S13_B_phenotype_eligibility.csv"), row.names = FALSE)
write.csv(phenotype_by_wave, file.path(output_dir, "S13_C_diagnosis_by_wave.csv"), row.names = FALSE)
write.csv(transition_counts, file.path(output_dir, "S13_C_detailed_transitions.csv"), row.names = FALSE)
write.csv(age_summary, file.path(output_dir, "S13_D_age_reconstruction.csv"), row.names = FALSE)
write.csv(checks, file.path(output_dir, "S13_object_checks.csv"), row.names = FALSE)
write.csv(genetic_qc_reference,
          file.path(output_dir, "S13_genetic_QC_manuscript_reference.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 7. Figure 1C flowchart (generated ONLY from the audited participant flow)
# -------------------------------------------------------------------------

figure1c_data <- flow %>%
  mutate(label = paste0(step, "\nN = ", format(N_remaining, big.mark = ",")),
         exclusion_label = ifelse(is.na(N_excluded), "",
                                  paste0("−", N_excluded, " excluded")),
         y = rev(seq_len(n())))

write.csv(figure1c_data,
          file.path(output_dir, "Figure1C_audited_nodes.csv"), row.names = FALSE)

p_flow <- ggplot(figure1c_data, aes(x = 0, y = y)) +
  geom_rect(aes(xmin = -2.2, xmax = 2.2, ymin = y - .32, ymax = y + .32),
            fill = "#F0F4F5", color = "#233B43", linewidth = .5) +
  geom_segment(data = head(figure1c_data, -1),
               aes(x = 0, xend = 0, y = y - .33, yend = y - .67),
               inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(.12, "inches")),
               linewidth = .45, color = "#828282") +
  geom_text(aes(label = label), size = 3.5, lineheight = .95) +
  geom_text(aes(x = 2.4, y = y + .52, label = exclusion_label),
            hjust = 0, color = "#555555", size = 3) +
  coord_cartesian(xlim = c(-2.5, 5.2), clip = "off") +
  scale_y_continuous(breaks = NULL) +
  theme_void(base_size = 11) +
  labs(caption = paste(
    "Unique-participant flow from the 2026 phenotype source.",
    "Genetic QC in the manuscript (2,511 to 2,190) is reported separately;",
    "no participant-level pre/post QC linkage was available for this audit.")) +
  theme(plot.caption = element_text(size = 8, hjust = 0, color = "#555555"),
        plot.margin = margin(15, 25, 15, 15))

ggsave(file.path(output_dir, "Figure1C_selection_flow.pdf"),
       p_flow, width = 10, height = 10.5)
ggsave(file.path(output_dir, "Figure1C_selection_flow.png"),
       p_flow, width = 10, height = 10.5, dpi = 300)

# -------------------------------------------------------------------------
# 8. SINGLE publication-style Supplementary Table S13 (sections A-D)
# -------------------------------------------------------------------------
fmt_s13_n <- function(x) {
  ifelse(is.na(x), "—", format(as.integer(x), big.mark = ",", trim = TRUE))
}

s13_section <- function(title) {
  tibble(Characteristic = title, Overall = "", W0 = "", W1 = "",
         W2 = "", W3 = "", .section = TRUE)
}

s13_overall <- function(label, text) {
  tibble(Characteristic = label, Overall = text, W0 = "", W1 = "",
         W2 = "", W3 = "", .section = FALSE)
}

s13_wave <- function(label, source, field) {
  stopifnot(field %in% names(source), "wave" %in% names(source))
  if (anyDuplicated(source$wave) || !all(waves %in% source$wave)) {
    stop("S13 source must contain exactly one row per wave for: ", label)
  }
  vals <- fmt_s13_n(source[[field]][match(waves, source$wave)])
  tibble(Characteristic = label, Overall = "", W0 = vals[1],
         W1 = vals[2], W2 = vals[3], W3 = vals[4], .section = FALSE)
}

# A: each arrow is before -> retained, with the number excluded in parentheses.
s13_A <- purrr::map_dfr(seq_len(nrow(flow)), function(i) {
  value <- if (is.na(flow$N_before[i])) {
    fmt_s13_n(flow$N_remaining[i])
  } else {
    paste0(fmt_s13_n(flow$N_before[i]), " → ",
           fmt_s13_n(flow$N_remaining[i]),
           " (−", fmt_s13_n(flow$N_excluded[i]), ")")
  }
  s13_overall(flow$step[i], value)
})

# B: each cell represents BEFORE / RETAINED AFTER RECONSTRUCTION / IN FINAL.
s13_B <- purrr::map_dfr(seq_len(nrow(group_counts)), function(i) {
  s13_overall(group_counts$participant_group[i], paste0(
    fmt_s13_n(group_counts$N_before_reconstruction[i]), " / ",
    fmt_s13_n(group_counts$N_retained_after_reconstruction[i]), " / ",
    fmt_s13_n(group_counts$N_in_final_analytic_sample[i])))
})

# C: all diagnosis summary variables, explicitly labeled by assessment wave.
s13_C_fields <- c(
  "Participants retained for phenotype reconstruction" = "N_retained",
  "Observed ADHD cases" = "N_observed_case",
  "Observed controls" = "N_observed_control",
  "Observed uncertain ADHD status" = "N_observed_uncertain",
  "Observed missing/other ADHD status" = "N_observed_missing_other",
  "ADHD cases after reconstruction" = "N_reconstructed_case",
  "Controls after reconstruction" = "N_reconstructed_control",
  "ADHD statuses changed" = "N_status_changed",
  "Cases carried forward from earlier waves" = "N_case_carried_forward",
  "Observed controls reclassified as cases" = "N_observed_control_to_case",
  "Uncertain/missing status assigned as cases" = "N_uncertain_missing_to_case",
  "Uncertain/missing status assigned as controls" = "N_uncertain_missing_to_control",
  "Participants in final analytical cohort" = "N_final",
  "ADHD cases in final analytical cohort" = "N_cases_final",
  "Controls in final analytical cohort" = "N_controls_final"
)
s13_C <- purrr::map2_dfr(names(s13_C_fields), unname(s13_C_fields),
                         ~ s13_wave(.x, phenotype_by_wave, .y))

# D: all eight age-source and completeness measures.
s13_D_fields <- c(
  "Participants retained for age reconstruction" = "N_participants",
  "Ages originally available" = "N_age_originally_available",
  "Ages available after reconstruction" = "N_age_after",
  "Ages derived from assessment dates" = "N_date_derived",
  "Original age values used" = "N_original_age_used",
  "Ages estimated from baseline offsets" = "N_baseline_offset_used",
  "Previously missing ages newly filled" = "N_newly_filled_age",
  "Ages still missing" = "N_still_missing_age"
)
s13_D <- purrr::map2_dfr(names(s13_D_fields), unname(s13_D_fields),
                         ~ s13_wave(.x, age_summary, .y))

s13_combined <- dplyr::bind_rows(
  s13_section("A. Participant inclusion and exclusions"), s13_A,
  s13_section("B. Longitudinal phenotype eligibility"), s13_B,
  s13_section("C. Observed and reconstructed ADHD status"), s13_C,
  s13_section("D. Age reconstruction"), s13_D
)

cat("\n=== Supplementary Table S13 - Combined ===\n")
print(dplyr::select(s13_combined, -.section), n = Inf, width = Inf)
write.csv(dplyr::select(s13_combined, -.section),
          file.path(output_dir, "S13_combined.csv"), row.names = FALSE,
          na = "")

# -------------------------------------------------------------------------
# 9. SINGLE native Word table (replaces the former four-subtable DOCX)
# -------------------------------------------------------------------------

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {

  section_rows <- which(s13_combined$.section)
  ft <- flextable::flextable(dplyr::select(s13_combined, -.section))
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = 8.5, part = "all")
  ft <- flextable::fontsize(ft, size = 9, part = "header")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::padding(ft, padding = 3, part = "all")
  ft <- flextable::align(ft, j = "Characteristic", align = "left", part = "all")
  ft <- flextable::align(ft, j = c("Overall", waves),
                         align = "center", part = "all")
  ft <- flextable::valign(ft, valign = "center", part = "all")

  # Make the section labels full-width rows inside the SAME table.
  for (i in section_rows) {
    ft <- flextable::merge_at(ft, i = i, j = 1:6, part = "body")
    ft <- flextable::bold(ft, i = i, bold = TRUE, part = "body")
    ft <- flextable::hline_top(
      ft, i = i, border = officer::fp_border(color = "666666", width = 0.7),
      part = "body")
  }
  ft <- flextable::set_table_properties(
    ft, layout = "fixed", opts_word = list(repeat_headers = TRUE)
  )
  ft <- flextable::width(ft, j = "Characteristic", width = 3.3)
  ft <- flextable::width(ft, j = "Overall", width = 2.5)
  ft <- flextable::width(ft, j = waves, width = 0.96)
  ft <- flextable::set_header_labels(
    ft, Characteristic = "Characteristic / processing step", Overall = "Overall")

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Supplementary Table S13",
                               style = "heading 1")
  doc <- officer::body_add_par(
    doc,
    "Participant filtering, longitudinal ADHD phenotype reconstruction, and age estimation",
    style = "Normal")
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_add_par(doc,
    paste(
      "Notes. A: before → remaining (excluded) at each consecutive step;",
      "the first row gives the starting N. B: before reconstruction / retained",
      "after reconstruction / present in the final analytical cohort. C and D",
      "refer to the same retained participants before final genetic/PCA filtering",
      "unless explicitly labeled 'final analytical cohort'. Observed ADHD",
      "diagnoses are distinct from reconstructed phenotypes. Case status persists",
      "after the first observed positive wave; eligible controls require four",
      "observed negative waves. Ages derived from dates or follow-up offsets are",
      "not new participant assessments. Published genetic QC: 2,511 before QC,",
      "2,190 after QC (321 excluded); these historical totals were not part of",
      "the participant-ID-level filtering audit."),
    style = "Normal"
  )
  doc <- officer::body_set_default_section(
    doc,
    officer::prop_section(
      page_size = officer::page_size(orient = "landscape"),
      page_margins = officer::page_mar(top = 0.65, bottom = 0.65,
                                     left = 0.65, right = 0.65)
    )
  )
  print(doc, target = file.path(output_dir, "sup_tab13.docx"))
  message("One-table Word document exported: ",
          file.path(output_dir, "sup_tab13.docx"))
} else {
  message("Word export skipped: install.packages(c('flextable','officer'))")
}

cat("\nOutputs saved in: ", normalizePath(output_dir), "\n", sep = "")
