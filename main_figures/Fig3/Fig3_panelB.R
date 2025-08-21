pacman::p_load(dplyr, tidyr, survival, envalysis, survminer, ggplot2, ggthemr, broom)
# Survplot for number os diagnosis | cox HR
data <- readRDS("E:/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS), by = "IID") %>%
  select(-FID)

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <-
	cbind(select(database, -PRS), PRS = new_PRS) %>%
  mutate(
    risk = ntile(PRS, 100),
    percentile = case_when(
      risk >= 90 ~ "90th",
      risk <= 10 ~ "10th",
      TRUE ~ "else"),
    percentile = factor(
      percentile,
      levels = c("10th", "else", "90th"))) %>%
  select(IID, site, W0, W1, W2, age_W0, age_W1, age_W2, percentile, gender, site)

# Mirror
mi <-
	readRDS("E:/0_external_files/Lucas_MINI_BHRCS.rds") %>%
	select(ident, IID) %>%
	mutate(ident = as.numeric(ident))

# ANX phenotype (PGC)
anx_PGC <-
	readRDS("E:/0_external_files/Santoro_192BHRC_2024_08_17.rds") %>%
	select(ident, dcanyanx_pgc, redcap_event_name) %>%
	inner_join(mi, by = "ident")

# Unprocessed phenotypes (includes dcanyhk)
# Warning message:
# In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#   Loglik converged before variable  1,2 ; coefficient may be infinite.
phenotype <-
	readRDS("E:/0_external_files/dawba_20200526.rds") %>%
	select(-dcanyanx) %>%
	mutate(subjectid = gsub("^", "C", subjectid)) %>%
	rename(IID = 1) %>%
	inner_join(anx_PGC, by = c("IID", "redcap_event_name")) %>%
	select(-ident) %>%
	select(
		IID, redcap_event_name, dcptsd, dcocd,
		dcanyanx_pgc, dcmadep, dcanyhk, dcodd,
		dceat, dcpsych) %>%
	mutate(
		across(dcptsd:dcpsych, ~ifelse(. == 0, 0, 1)),
		rx_number = rowSums(across(dcptsd:dcpsych), na.rm = TRUE),
		across(dcptsd:rx_number, ~as.factor(.)),
		redcap_event_name = case_when(
			redcap_event_name == "wave0_arm_1" ~ "W0",
			redcap_event_name == "wave1_arm_1" ~ "W1",
			redcap_event_name == "wave2_arm_1" ~ "W2",
			TRUE ~ NA)) %>%
	rename(wave = 2) %>%
	data.frame()

# Diagnosis number
rx_only <-
	select(phenotype, IID, wave, rx_number) %>%
	mutate(
    wave = factor(wave, levels = c("W0", "W1", "W2")),
    rx_number = as.numeric(as.character(rx_number))) %>%
  group_by(IID) %>%
  summarise(max_rx = max(rx_number)) %>%
  filter(IID %in% database$IID)
str(rx_only)
# get the highest one

## keep controls the same
without_entry <-
  filter(database, W2 == 0) %>%
  select(IID, age_W2, W2) %>%
  rename(time = 2, status = 3)
str(without_entry)

## with the first occurance
temp1 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, W0, W1, W2) %>%
  pivot_longer(
    cols = starts_with("W"),
    names_to = "wave",
    values_to = "diagnosis")
str(temp1)

## Age data
temp2 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, age_W0, age_W1, age_W2) %>%
  pivot_longer(
    cols = starts_with("age_W"),
    names_to = "wave",
    values_to = "age") %>%
  mutate(wave = gsub("age_", "", wave))
str(temp2)

with_entry <-
  inner_join(temp1, temp2, by = c("IID", "wave")) %>%
  filter(diagnosis == 1) %>% # any time diagnosis
  group_by(IID) %>%
  filter(age == min(age)) %>%
  ungroup() %>%
  select(-wave) %>%
  rename(status = 2, time = 3)
str(with_entry)

# add here the rx_number
survival_data <-
  rbind(with_entry, without_entry) %>%
  inner_join(., select(database, IID, site, percentile, gender, site), by = "IID") %>%
  data.frame() %>%
  inner_join(., rx_only, by = "IID") %>%
  select(-IID)
str(survival_data)

# the only difference, is that I have put the percentile
for_plot <-
  rbind(
    tidy(
      coxph(Surv(time, status) ~ percentile + gender + site + max_rx, data = survival_data),
      exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(tag = "all"),
    tidy(
      coxph(Surv(time, status) ~ percentile + gender + site, data = filter(survival_data, max_rx == 2)),
      exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(tag = "2"),
    tidy(
      coxph(Surv(time, status) ~ percentile + gender + site, data = filter(survival_data, max_rx >= 3)),
      exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(tag = ">=3")) %>%
  filter(!term %in% c("percentileelse", "siteRS")) %>%
  mutate(
    term = case_when(
      term == "genderMale" ~ "Gender",
      term == "percentile90th" ~ "90th",
      term == "percentile10th" ~ "10th",
      term == "max_rx" ~ "RX_number",
      TRUE ~ term),
  stars = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""),
  estimate = round(estimate, 2),
  CI = paste0(round(conf.low, 2), "—", round(conf.high, 2), stars),
  tag = factor(tag, levels = c("all", "2", ">=3")))

# o número de diagnósticos poderia ter retirado o ADHD
# pode ser um viés
ggthemr("grape")
final <-
  ggplot(for_plot, aes(term, estimate, fill = term, color = term)) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.5, position = position_dodge(width = 1)) +
    geom_col() +
    geom_text(
      aes(label = estimate, hjust = 0.5, vjust = 1.5),
      color = "white", size = 3,
			position = position_dodge(width = 1)) +
    geom_text(
      aes(y = conf.high + 1, label = stars),
      position = position_dodge(width = 1),
      vjust = 0.7, size = 5, show.legend = FALSE, angle = 90) +
    labs(y = "Harzard Ratio", x = "") +
    facet_wrap(~ tag, nrow = 1, strip.position = "top", scales = "free_x") +
    theme_publish(base_family = 7) +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_blank(),
      legend.title = element_blank())

ggsave(
  "Fig3_panelB.png", final, device = "png",
  width = 100, height = 60, units = "mm",
  dpi = 300, bg = "white")
