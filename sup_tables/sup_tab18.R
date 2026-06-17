# PRS interaction verification with Odds ratio model
pacman::p_load(dplyr, flextable, tidyr, broom, purrr)

# Article dataset
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# Diagnosis standardtization
database <-
	data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	mutate(across(c(W0, W1, W2), as.factor)) %>%
	inner_join(., hist, by = "IID")

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))

wd <-
  cbind(select(database, -PRS), PRS = new_PRS) %>%
  rename(
    diagnosis_W0 = W0,
    diagnosis_W1 = W1,
    diagnosis_W2 = W2) %>%
  select(
    IID, gender, site,
    starts_with("diagnosis"),
    starts_with("age_"),
    any_hist, PRS)

# PRS gender interaction test
# ---------------
# all samples
# ---------------
m_w0 <- glm(
  diagnosis_W0 ~ PRS * gender + age_W0 + site + any_hist,
  data = wd,
  family = "binomial")

m_w1 <- glm(
  diagnosis_W1 ~ PRS * gender + age_W1 + site + any_hist,
  data = wd,
  family = "binomial")

m_w2 <- glm(
  diagnosis_W2 ~ PRS * gender + age_W2 + site + any_hist,
  data = wd,
  family = "binomial")

mods <- list(
  W0 = m_w0,
  W1 = m_w1,
  W2 = m_w2)

# Controlar efeito aleatorio
# glmmTMB(
#   percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid)
#    data = nanosight_intersect,
#    family = beta_family(link = logit))
# dispformula = age

interactions <-
  imap_dfr(mods, ~
    tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "PRS:genderMale") %>%
    mutate(wave = .y)) %>%
  select(wave, estimate, conf.low, conf.high, p.value) %>%
	mutate(
		estimate = sprintf("%.2f", estimate),
		conf.low = sprintf("%.2f", conf.low),
		conf.high = sprintf("%.2f", conf.high),
  	p.value = format.pval(p.value, digits = 3, eps = .001))

# final table
final <- interactions %>%
	flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::align(part = "all", align = "center") %>%
	flextable::autofit()

# Fix table
# flextable::save_as_docx(
#   "Supplementary Table S18" = final,
#   path = "sup_tab18.docx")
