pacman::p_load(dplyr, tidyr, ggthemr, envalysis, ggplot2)
# Survplot for number os diagnosis | cox HR
data <- readRDS("E:/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

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

# Unprocessed phenotypes
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
	mutate(wave = factor(wave, levels = c("W0", "W1", "W2"))) %>%
	group_by(wave, rx_number) %>%
	summarise(N = n()) %>%
	filter(!rx_number == 0)
str(rx_only)

# Plot
# make so it is stacked by wave
ggthemr("grape")
p1 <-
	ggplot(rx_only, aes(wave, N, fill = rx_number)) +
		geom_bar(stat = "identity", position = "dodge") +
		scale_y_continuous(n.breaks = 20) +
		theme_publish(base_size = 7) +
		theme(
			legend.position = "right",
			panel.grid.major.x = element_line(
				linewidth = 0.2,
				color = "#a1a1a1",
				linetype = "dashed")) +
		coord_flip()

ggsave(
  "Fig3_panelA.png", p1, device = "png",
  width = 150, height = 40, units = c("mm"),
  dpi = 300, bg = "white")
