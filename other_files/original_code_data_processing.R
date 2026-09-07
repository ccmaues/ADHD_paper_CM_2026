pacman::p_load(lubridate, dplyr, data.table, tidyr)
# decode
mi <-
  readRDS("E:/0_external_files/Lucas_MINI_BHRCS.rds") %>%
  select(1, 4) %>%
  mutate(ident = as.numeric(ident))
str(mi)
# new wave

w3_vars <- c(
  # identifiers
  "ident",
  "probid",
  "redcap_event_name",
  "wave",

  # demographics
  "site",
  "gender",

  # dates / age
  "d_date",
  "birth_date",
  "age",

  # anthropometry
  "p_height",
  "p_height1",
  "p_height2",

  # primary ADHD phenotype
  "dcanyhk",

  # ADHD subtype diagnoses -- keep for QC/future analyses
  "dcadhdc",
  "dcadhdi",
  "dcadhdh",
  "dcadhdo")

bhrc <-
  readRDS("E:/0_external_files/Santoro_192BHRC_2024_08_17.rds") %>%
  select(ident, gender, birth_date, d_date, age, site) %>%
  inner_join(., mi, by = "ident") %>%
  mutate(
    site = case_when( # changed as requested
      site == 1 ~ "RS",
      site == 2 ~ "SP",
      is.na(site) ~ "NA",
      .default = "Outro"),
    gender = case_when( # changed as requested
      gender == 2 ~ "Female",
      gender == 1 ~ "Male",
      is.na(gender) ~ "NA",
      .default = "Outro"),
    d_date = as.Date(d_date),
    year = year(d_date),
    wave = case_when(
      year < 2012 ~ "W0",
      year >= 2012 & year < 2015 ~ "W1",
      year >= 2015 ~ "W2",
      is.na(year) ~ NA_character_)) %>%
  select(-year) %>%
  filter(!is.na(wave)) %>% # remove covid waves
  pivot_wider(
    names_from = wave,
    values_from = c("age", "d_date"),
    names_glue = "{.value}_{wave}") %>%
  data.frame()
head(bhrc)
# height
h <-
  readRDS("E:/objects_R/Antrop.rds") %>%
  select(ident, redcap_event_name, p_height) %>%
  inner_join(., mi, by = "ident") %>%
  mutate(
    redcap_event_name = case_when(
      redcap_event_name == "wave0_arm_1" ~ "W0",
      redcap_event_name == "wave1_arm_1" ~ "W1",
      redcap_event_name == "wave2_arm_1" ~ "W2")) %>%
  rename(ident = 1, wave = 2, height = 3, IID = 4) %>%
  pivot_wider(
    names_from = wave,
    values_from = height,
    id_cols =  "IID") %>%
  rename(height = 2)
str(h)
# phenotype should be in the second object
pheno <-
  readRDS("E:/0_external_files/dawba_20200526.rds") %>%
  select(subjectid, redcap_event_name, dcanyhk) %>%
  rename(IID = 1, wave = 2, ADHD = 3) %>%
  mutate(wave = case_when(
    wave == "wave0_arm_1" ~ "W0",
    wave == "wave1_arm_1" ~ "W1",
    wave == "wave2_arm_1" ~ "W2",
    TRUE ~ NA_character_)) %>%
  data.frame() %>%
  pivot_wider(
    names_from = wave,
    values_from = ADHD) %>%
  select(IID, W0, W1, W2)
pheno$IID <- gsub("^", "C", pheno$IID)
str(pheno)
# PRS remains unchanged
prs <-
fread("E:/PRS_database_imputed/PRSCS_ADHD_2023_Score.profile") %>%
  select(ident, PRSCS) %>%
  inner_join(., mi, by = "ident") %>%
  select(-ident) %>%
  rename(PRS = 1)
str(prs)
# Joining all
data <- plyr::join_all(
  list(pheno, h, bhrc, prs),
  by = "IID",
  type = "full") %>%
  mutate(
   site = factor(site, levels = c("SP", "RS")),
   gender = factor(gender, levels = c("Female", "Male")),
   across(W0:W2, factor))
str(data)

#saveRDS(data, "/media/santorolab/C207-3566/objects_R/cass_BHRC_without_imputation_21112024.RDS")
t <- data %>%
  select(IID, W0, W1, W2) %>%
  mutate(
    W0 = as.numeric(as.character(W0)),
    W1 = as.numeric(as.character(W1)),
    W2 = as.numeric(as.character(W2))) %>%
  mutate(
    sum1 = rowSums(select(., W0, W1, W2), na.rm = TRUE), # aqui
    sum2 = rowSums(select(., W0, W1, W2), na.rm = FALSE)) # aqui
str(t) # data used

case1 <- filter(t, sum1 == 6) # case in all waves

imputed_pheno <-
  data.frame(IID = case1$IID, W0 = 2, W1 = 2, W2 = 2)

case2 <- # case from w0
  filter(t, W0 == 2) %>%
  filter(!IID %in% imputed_pheno$IID)

imputed_pheno <- bind_rows(
  imputed_pheno,
  data.frame(IID = case2$IID, W0 = 2, W1 = 2, W2 = 2))

case3 <- # case in W1
  filter(t, W1 == 2) %>%
  filter(!IID %in% imputed_pheno$IID)

imputed_pheno <- bind_rows(
  imputed_pheno,
  data.frame(IID = case3$IID, W0 = 0, W1 = 2, W2 = 2))

case4 <- # case in W2
  filter(t, W2 == 2) %>%
  filter(!IID %in% imputed_pheno$IID)

imputed_pheno <- bind_rows(
  imputed_pheno,
  data.frame(IID = case4$IID, W0 = 0, W1 = 0, W2 = 2))

# controls that came in all waves
case5 <- filter(t, sum2 == 0)

imputed_pheno <- bind_rows(
  imputed_pheno,
  data.frame(IID = case5$IID, W0 = 0, W1 = 0, W2 = 0)) %>%
  inner_join(., mi, by = "IID") %>%
  select(-ident)

str(imputed_pheno)
# did not estimated the age at W0 with birthday
# used the one put in redcap
new_age <- 
  select(data, IID, birth_date, age_W0, d_date_W0, d_date_W1, d_date_W2) %>%
  mutate(
    age_W1 = as.numeric(difftime(d_date_W1, birth_date, units = "days") / 365.25),
    age_W2 = as.numeric(difftime(d_date_W2, birth_date, units = "days") / 365.25),
   # first eval 2011 and last 2014
    age_W1 = case_when(is.na(age_W1) ~ age_W0 + 3, TRUE ~ age_W1),
   # first eval 2011 and last 2019
    age_W2 = case_when(is.na(age_W2) ~ age_W0 + 8, TRUE ~ age_W2)) %>%
  select(IID, age_W0, age_W1, age_W2)
pca_file_IDs <- 
  select(data, -W0, -W1, -W2, -age_W0, -age_W1, -age_W2) %>%
  inner_join(., imputed_pheno, new_age, by = "IID") %>%
  filter(!is.na(PRS)) %>%
  select(IID)
fwrite(pca_file_IDs, file = "cass_BHRC_IIDs_for_PCA.txt", col.names = FALSE, eol = "\n")
all_pca <- fread("E:/PCA_files_cass/cass_final_PCA/all_samples_PCA.eigenvec")
fem_pca <- fread("E:/PCA_files_cass/cass_final_PCA/all_females_PCA.eigenvec")
man_pca <- fread("E:/PCA_files_cass/cass_final_PCA/all_males_PCA.eigenvec")
colnames(all_pca) <- c("FID", "IID", paste0("PC", 1:10, sep = ""))
colnames(fem_pca) <- c("FID", "IID", paste0("PC", 1:10, sep = ""))
colnames(man_pca) <- c("FID", "IID", paste0("PC", 1:10, sep = ""))
family_hist <-
  readRDS("E:/0_external_files/Lucas_MINI_BHRCS.rds") %>%
  select(IID, starts_with("imini_")) %>%
  mutate(
      parent_mania = ifelse(imini_cur_man != 0 | imini_lif_man != 0, 1, 0),
      parent_dep = ifelse(imini_cur_dep != 0 | imini_cur_rdep != 0, 1, 0),
      parent_panic = ifelse(imini_cur_panic != 0 | imini_lif_panic != 0, 1, 0), # should this not enter the anx_any?
      parent_psych = ifelse(imini_cur_psych != 0 | imini_lif_psych != 0, 1, 0),
      parent_adhd = ifelse(imini_cur_adhd != 0 | imini_child_adhd != 0, 1, 0),
      parent_alc = ifelse(imini_cur_alcdep != 0 | imini_cur_alcabus != 0, 1, 0),
      parent_drug = ifelse(imini_cur_drugdep != 0 | imini_cur_drugabus != 0, 1, 0),
      across(is.numeric, as.factor)) %>%
  select(1, starts_with("parent_"), 12) %>%
  rename(parent_anx = 9)
str(family_hist)
cass_BHRC <- list(
  proband_data =
    plyr::join_all(
    list(
      select(data, -W0, -W1, -W2, -age_W0, -age_W1, -age_W2, -ident) %>% filter(IID %in% all_pca$IID),
      imputed_pheno,
      new_age,
      prs),
    by = "IID", type = "inner"),
  family_history = filter(family_hist, IID %in% all_pca$IID),
  PCA_all_samples = all_pca,
  PCA_by_sex = rbind(fem_pca, man_pca))

str(cass_BHRC)
saveRDS(cass_BHRC, "E:/cass_BHRC_28042025_ARTICLE.RDS")
