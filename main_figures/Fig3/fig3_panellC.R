# geom_point do PRS por grupo (talvez mudar pra painel B
# ao invés de C)
options(scipen = 999) # disable scientific notation
pacman::p_load(ggplot2, envalysis, dplyr, ggthemr)
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
  select(IID, site, W0, W1, W2, age_W0, age_W1, age_W2, PRS, percentile, gender, site)

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

# Unprocessed phenotypes (dcanyhk excluded)
phenotype <-
	readRDS("E:/0_external_files/dawba_20200526.rds") %>%
	select(-dcanyanx) %>%
	mutate(subjectid = gsub("^", "C", subjectid)) %>%
	rename(IID = 1) %>%
	inner_join(anx_PGC, by = c("IID", "redcap_event_name")) %>%
	select(-ident) %>%
	select(
		IID, redcap_event_name, dcptsd, dcocd, #dcanyhk,
		dcanyanx_pgc, dcmadep, dcodd, dceat, dcpsych) %>%
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

# 90th associated w/ XXX etc
for_plot <-
	inner_join(database, rx_only, by = "IID") %>%
	select(PRS, percentile, max_rx)

glm(PRS ~ max_rx*percentile, family = "gaussian", data = for_plot) %>%
broom::tidy()

# plot

# ggplot(for_plot, aes())