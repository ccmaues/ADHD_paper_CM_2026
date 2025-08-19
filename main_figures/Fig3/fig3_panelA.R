# spaguetti plot of RX number per wave
data <- readRDS("/home/santorolab/Desktop/cassia/pendrive_cass_BK/cass_BHRC_28042025_ARTICLE.RDS")

mi <-
	readRDS("/media/santorolab/C207-3566/0_external_files/Lucas_MINI_BHRCS.rds") %>%
	select(ident, IID) %>%
	mutate(ident = as.numeric(ident))

anx_PGC <-
	readRDS("/media/santorolab/C207-3566/0_external_files/Santoro_192BHRC_2024_08_17.rds") %>%
	select(ident, dcanyanx_pgc, redcap_event_name) %>%
	inner_join(mi, by = "ident")

phenotype <-
	readRDS("/media/santorolab/C207-3566/0_external_files/dawba_20200526.rds") %>%
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

rx_only <-
  select(phenotype, IID, wave, rx_number) %>%
  pivot_wider(names_from = "wave", values_from = "rx_number") %>%
  rename(rx_W0 = W0, rx_W1 = W1, rx_W2 = W2) %>%
	filter(IID %in% data$proband_data$IID) %>%
	pivot_longer(
		cols = starts_with("rx_W"),
		names_to = "wave",
		values_to = "rx") %>%
	mutate(wave = gsub("rx_", "", wave))

for_plot <-
	rx_only %>%
	group_by(wave, rx) %>%
	summarise(N = n())

ggthemr("fresh")
ggplot(for_plot, aes(wave, N, color = rx, group = rx)) +
geom_point() +
geom_path() +
scale_y_continuous(n.breaks = 10) +
theme_publish()