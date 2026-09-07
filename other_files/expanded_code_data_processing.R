pacman::p_load(lubridate, dplyr, data.table, tidyr)

# -------------------------------------------------------------------------
# 2. Paths
# -------------------------------------------------------------------------
external_dir <- "C:/Users/cassi/Documents/work/0_external_files"
prs_file <- "C:/Users/cassi/Documents/work/PRS_database_imputed/PRSCS_ADHD_2023_Score.profile"
output_file <- "C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds"

# -------------------------------------------------------------------------
# 3. ID dictionary / mirror
# -------------------------------------------------------------------------
lucas_db <- readRDS(file.path(external_dir, "Lucas_MINI_BHRCS.rds"))

mi <-
  lucas_db %>%
  select(1, 4) %>%
  mutate(ident = as.numeric(ident)) %>%
  inner_join(.,
    fread(file.path(external_dir, "cass_1653_IDs_IID.tsv")),
    by = "IID")

str(mi)

# -------------------------------------------------------------------------
# 4. Load longitudinal BHRC data
#    New source contains demographics, dates, height and ADHD phenotype
# -------------------------------------------------------------------------
old_bhrc <-
  readRDS(file.path(external_dir, "Santoro_192BHRC_2024_08_17.rds"))

old_bhrc_list <-
  unclass(old_bhrc)

dob_info <-
  data.frame(
    ident = old_bhrc_list$ident,
    birth_date = as.Date(
      as.integer(unclass(old_bhrc_list$birth_date)),
      origin = "1970-01-01"),
    redcap_event_name = old_bhrc_list$redcap_event_name) %>%
  filter(redcap_event_name == "wave0_arm_1") %>%
  select(-redcap_event_name) %>%
  inner_join(., mi, by = "ident") %>%
  select(-ident) %>%
  distinct(ext_genid, .keep_all = TRUE)

str(dob_info)

bhrc_source <- readRDS(
  file.path(external_dir, "Santoro_278BHRC_2026_04_23.rds")) %>%
  inner_join(., dob_info, by = "ext_genid")

required_vars <- c(
  "ext_genid", "redcap_event_name", "site", "gender",
  "birth_date", # not in the new database -- gather from previous version
  "p_psico_date1", "age", "p_height", "dcanyhk")

missing_vars <- setdiff(required_vars, names(bhrc_source))

if (length(missing_vars) > 0) {
  stop(
    "Missing required variables in Santoro_278BHRC_2026_04_23.rds: ",
    paste(missing_vars, collapse = ", "))}

# -------------------------------------------------------------------------
# 5. Keep only true longitudinal waves W0-W3
#    COVID-specific W2 events are intentionally excluded
# -------------------------------------------------------------------------
bhrc_raw <-
  bhrc_source %>%
  filter(
    redcap_event_name %in% c(
      "wave0_arm_1",
      "wave1_arm_1",
      "wave2_arm_1",
      "wave3_arm_1")) %>%
  select(
    all_of(required_vars),
    any_of(c("p_height1", "p_height2", "dcadhdc", "dcadhdi", "dcadhdh", "dcadhdo"))) %>%
  mutate(
    wave = case_when(
      redcap_event_name == "wave0_arm_1" ~ "W0",
      redcap_event_name == "wave1_arm_1" ~ "W1",
      redcap_event_name == "wave2_arm_1" ~ "W2",
      redcap_event_name == "wave3_arm_1" ~ "W3",
      TRUE ~ NA_character_),
    p_psico_date1 = as.Date(p_psico_date1),
    birth_date = as.Date(birth_date)) %>%
  inner_join(mi, by = "ext_genid")

str(bhrc_raw)
# -------------------------------------------------------------------------
# 6. Participant-level demographics
# -------------------------------------------------------------------------
first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]}

demographics <-
  bhrc_raw %>%
  group_by(IID) %>%
  summarise(
    IID = first_nonmissing(IID),
    site = first_nonmissing(site),
    gender = first_nonmissing(gender),
    birth_date = first_nonmissing(birth_date),
    .groups = "drop") %>%
  mutate(
    site = case_when(
      site == 1 ~ "RS",
      site == 2 ~ "SP",
      is.na(site) ~ NA_character_,
      TRUE ~ "Other"),
    gender = case_when(
      gender == 2 ~ "Female",
      gender == 1 ~ "Male",
      is.na(gender) ~ NA_character_,
      TRUE ~ "Other"))

str(demographics)
# -------------------------------------------------------------------------
# 7. Longitudinal variables: one W0-W3 pivot
# -------------------------------------------------------------------------
longitudinal_wide <-
  bhrc_raw %>%
  select(IID, wave, p_psico_date1, age, p_height, dcanyhk) %>%
  pivot_wider(
    id_cols = IID,
    names_from = wave,
    values_from = c(
      p_psico_date1,
      age,
      p_height,
      dcanyhk),
    names_glue = "{.value}_{wave}") %>%
  rename(
      W0 = dcanyhk_W0,
      W1 = dcanyhk_W1,
      W2 = dcanyhk_W2,
      W3 = dcanyhk_W3) %>%
  select(-p_psico_date1_W0) %>%
  left_join(dob_info, by = "IID") %>%
  select(-birth_date, -ext_genid)

str(longitudinal_wide)
# -------------------------------------------------------------------------
# 8. Optional W3 QC variables
# -------------------------------------------------------------------------

qc_vars <- intersect(
  c("p_height1", "p_height2", "dcadhdc", "dcadhdi", "dcadhdh", "dcadhdo"),
  names(bhrc_raw))

if (length(qc_vars) > 0) {
  w3_qc <-
    bhrc_raw %>%
    filter(wave == "W3") %>%
    select(
      ext_genid,
      all_of(qc_vars))
} else {
  w3_qc <- NULL}

str(w3_qc)
# -------------------------------------------------------------------------
# 9. PRS
# -------------------------------------------------------------------------
prs <-
  fread(prs_file) %>%
  select(ident, PRSCS) %>%
  mutate(ident = as.numeric(ident)) %>%
  inner_join(mi, by = "ident") %>%
  rename(PRS = PRSCS)

str(prs)
# -------------------------------------------------------------------------
# 10. Build raw analysis database
# -------------------------------------------------------------------------
data <-
  reduce(
    list(demographics, longitudinal_wide, prs),
    inner_join, by = "IID") %>%
  mutate(
    site = factor(site, levels = c("SP", "RS")),
    gender = factor(gender, levels = c("Female", "Male")),
    across(W0:W3, factor)) %>%
  select(-ident, -ext_genid)

str(data)
# -------------------------------------------------------------------------
# 11. ADHD phenotype reconstruction / imputation
# -------------------------------------------------------------------------
t <-
  data %>%
  select(IID, W0, W1, W2, W3) %>%
  mutate(
    across(
      W0:W3,
      ~ as.numeric(as.character(.x))))
# Case first observed at W0
case_w0 <-
  t %>%
  filter(W0 == 2)
str(case_w0)

# Case first observed at W1
case_w1 <-
  t %>%
  filter(
    W1 == 2,
    !IID %in% case_w0$IID)
str(case_w1)

# Case first observed at W2
case_w2 <-
  t %>%
  filter(
    W2 == 2,
    !IID %in% c(case_w0$IID, case_w1$IID))
str(case_w2)

# Case first observed at W3
case_w3 <-
  t %>%
  filter(
    W3 == 2,
    !IID %in% c(case_w0$IID, case_w1$IID, case_w2$IID))
str(case_w3)

# Controls observed as 0 at every wave
control_all <-
  t %>%
  filter(W0 == 0, W1 == 0, W2 == 0, W3 == 0)
str(control_all)

imputed_pheno <-
  bind_rows(
    data.frame(IID = case_w0$IID, W0 = 2, W1 = 2, W2 = 2, W3 = 2),
    data.frame(IID = case_w1$IID, W0 = 0, W1 = 2, W2 = 2, W3 = 2),
    data.frame(IID = case_w2$IID, W0 = 0, W1 = 0, W2 = 2, W3 = 2),
    data.frame(IID = case_w3$IID, W0 = 0, W1 = 0, W2 = 0, W3 = 2),
    data.frame(IID = control_all$IID, W0 = 0, W1 = 0, W2 = 0, W3 = 0)) %>%
  distinct(IID, .keep_all = TRUE) %>%
  inner_join(., mi, by = "IID") %>%
  select(-ident, -ext_genid)

str(imputed_pheno)
# -------------------------------------------------------------------------
# 12. Age reconstruction
# -------------------------------------------------------------------------
new_age <- select(data,
    IID, birth_date, age_W0, age_W1, age_W2, age_W3,
    p_psico_date1_W1, p_psico_date1_W2, p_psico_date1_W3) %>%
  mutate(
    age_calc_W1 =  as.numeric(difftime(p_psico_date1_W1,birth_date,units = "days")) / 365.25,
    age_calc_W2 =  as.numeric(difftime(p_psico_date1_W2,birth_date,units = "days")) / 365.25,
    age_calc_W3 =  as.numeric(difftime(p_psico_date1_W3,birth_date,units = "days")) / 365.25,
    age_W1 = coalesce(age_calc_W1, age_W1, age_W0 + 3),
    age_W2 = coalesce(age_calc_W2, age_W2, age_W0 + 8),
    age_W3 = coalesce(age_calc_W3, age_W3, age_W0 + 12)) %>%
  select(IID, age_W0, age_W1, age_W2, age_W3)

str(new_age)
# -------------------------------------------------------------------------
# 13. Final sample for PCA
# -------------------------------------------------------------------------
fam <-
  fread(file.path(external_dir, "BHRC_Probands_Final.fam")) %>%
  rename(FID = 1, IID = 2)

pca_files_IDs <-
  inner_join(
    select(imputed_pheno, IID),
    select(demographics, IID, gender),
    by = "IID") %>%
  inner_join(., fam, by = "IID") %>%
  select(FID, IID, gender)

str(pca_files_IDs)

pca_files_IDs_fem <- filter(pca_files_IDs, gender == "Female")
str(pca_files_IDs_fem)

pca_files_IDs_man <- filter(pca_files_IDs, gender == "Male")
str(pca_files_IDs_man)

# fwrite(
#   pca_files_IDs,
#   file = "cass_BHRC_W3_IIDs_for_PCA.txt",
#   col.names = FALSE,
#   eol = "\n",
#   sep = "\t")

# fwrite(
#   pca_files_IDs_fem,
#   file = "cass_FEMALES_BHRC_W3_IIDs_for_PCA.txt",
#   col.names = FALSE,
#   eol = "\n",
#   sep = "\t")

# fwrite(
#   pca_files_IDs_man,
#   file = "cass_MALES_BHRC_W3_IIDs_for_PCA.txt",
#   col.names = FALSE,
#   eol = "\n",
#   sep = "\t")
# -------------------------------------------------------------------------
# 14. PCA eigenvectors
# -------------------------------------------------------------------------
all_pca <- fread("C:/Users/cassi/Documents/work/0_external_files/w3_pca/cass_1506_W3_PCA.eigenvec")
colnames(all_pca) <- c("FID", "IID", paste0("PC", 1:10))

fem_pca <- fread("C:/Users/cassi/Documents/work/0_external_files/w3_pca/cass_697_W3_PCA.eigenvec")
colnames(fem_pca) <- c("FID", "IID", paste0("PC", 1:10))

man_pca <- fread("C:/Users/cassi/Documents/work/0_external_files/w3_pca/cass_806_W3_PCA.eigenvec")
colnames(man_pca) <- c("FID", "IID", paste0("PC", 1:10))

# -------------------------------------------------------------------------
# 15. Family history
# -------------------------------------------------------------------------
family_hist <-
  lucas_db %>%
  select(IID, starts_with("imini_")) %>%
  mutate(
    parent_mania = ifelse(imini_cur_man != 0 | imini_lif_man != 0, 1, 0),
    # Retained exactly from original code
    parent_dep = ifelse(imini_cur_dep != 0 | imini_cur_rdep != 0, 1, 0),
    parent_panic = ifelse(imini_cur_panic != 0 | imini_lif_panic != 0, 1, 0),
    parent_psych = ifelse(imini_cur_psych != 0 | imini_lif_psych != 0, 1, 0),
    parent_adhd = ifelse(imini_cur_adhd != 0 | imini_child_adhd != 0, 1, 0),
    parent_alc = ifelse(imini_cur_alcdep != 0 | imini_cur_alcabus != 0, 1, 0),
    parent_drug = ifelse(imini_cur_drugdep != 0 |  imini_cur_drugabus != 0, 1, 0),
    any_hist = as.integer(if_any(
      c(parent_mania, parent_dep,
        parent_panic, parent_psych,
        parent_adhd, parent_alc,
        parent_drug), ~ .x == 1)),
    across(where(is.numeric), as.factor)) %>%
  rename(parent_anx = imini_cur_any_anx) %>%
  select(1, starts_with("parent_"), 12, any_hist)

# add the any_hist here

str(family_hist)
# -------------------------------------------------------------------------
# 16. Final article database
# -------------------------------------------------------------------------
proband_data <-
  list(
    data %>%
      filter(IID %in% all_pca$IID) %>%
      select(-starts_with("W"), -starts_with("age_W")),
    imputed_pheno,
    new_age) %>%
  reduce(inner_join, by = "IID") %>%
  select(IID, site, gender, PRS, starts_with("age_"), starts_with("W")) %>%
  mutate(W3 = ifelse(W3 == 2, 1, 0))

str(proband_data)

cass_BHRC <-
  list(proband_data = proband_data,
    family_history =
      family_hist %>%
      filter(IID %in% all_pca$IID),
    PCA_all_samples = all_pca,
    PCA_by_sex = rbind(fem_pca, man_pca),
    W3_QC = w3_qc)

# -------------------------------------------------------------------------
# 17. Basic checks
# -------------------------------------------------------------------------
cat("\nFinal proband sample:", nrow(cass_BHRC$proband_data), "\n")

cat("\nADHD phenotype counts:\n")
print(
  lapply(
    cass_BHRC$proband_data[c("W0", "W1", "W2", "W3")],
    table,
    useNA = "ifany"))

cat("\nAge missingness:\n")
print(
  colSums(
    is.na(
      cass_BHRC$proband_data[
        c("age_W0", "age_W1", "age_W2", "age_W3")])))

# -------------------------------------------------------------------------
# 18. Save final RDS
# -------------------------------------------------------------------------
saveRDS(cass_BHRC, output_file)

cat(
  "\nSaved final W0-W3 database to:\n",
  output_file,
  "\n")
