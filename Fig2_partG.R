pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, broom)
options(scipen = 999) # disable scientific notation

# Dataset
data <- readRDS("/home/santorolab/Desktop/cassia/pendrive_cass_BK/cass_BHRC_28042025_ARTICLE.RDS")

# Proband data
database <-
  data$proband_data %>%
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS), by = "IID") %>%
  select(-FID)

# CM-PRS
# https://www.repository.cam.ac.uk/items/0cf1dc7c-0842-434f-b425-3bc5a9f33e05
mi <-
	readRDS("/home/santorolab/Desktop/cassia/pendrive_cass_BK/0_external_files/Lucas_MINI_BHRCS.rds") %>%
	select(ident, IID) %>%
	mutate(ident = as.numeric(ident))

cm_prs <-
	fread("/home/santorolab/Desktop/cassia/pendrive_cass_BK/PRS_database_imputed/PRSCS_CMT_2021_Score.profile") %>%
	select(ident, PRSCS) %>%
	inner_join(., mi, by = "ident") %>%
	select(-ident) %>%
	inner_join(., select(all_pcs, -PRS), by = "IID") %>%
	rename(PRS = PRSCS)

# PRS correction
shapiro.test(all_pcs$PRS)
new_ADHD <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
new_CM <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = cm_prs))
database <- cbind(select(database, -PRS), PRS_CM = new_CM, PRS_ADHD = new_ADHD)

# Family history
parents <- data$family_history

# Working data
# test_CM <-
# 	inner_join(database, parents, by = "IID") %>%
# 	select(PRS_CM, starts_with("parent_"))

# test_ADHD <-
# 	inner_join(database, parents, by = "IID") %>%
# 	select(PRS_ADHD, starts_with("parent_"))

inner_join(database, parents, by = "IID") %>%
	mutate(
		diagnosis = case_when(
			parent_psych == 1 ~ "psych",
			parent_mania == 1 ~ "mania",
			parent_dep == 1 ~ "dep",
			parent_panic == 1 ~ "panic",
			parent_adhd == 1 ~ "adhd",
			parent_alc == 1 ~ "alc",
			parent_drug == 1 ~ "drug",
			parent_anx == 1 ~ "anx",
			.default = no_diagnosis))
			# falta a comorbidade

# PRS load: PRS association testing of groups
# based on the parent diagnosis
plot_object <-
	rbind(
		glm(PRS_ADHD ~ ., data = test_ADHD, family = "gaussian") %>%
		tidy(conf.int = TRUE) %>%
		mutate(phenotype = "ADHD"),
		glm(PRS_CM ~ ., data = test_CM, family = "gaussian") %>%
		tidy(conf.int = TRUE) %>%
		mutate(phenotype = "CM")) %>%
	filter(!term == "(Intercept)") %>%
	mutate(
		term = gsub("parent_", "",  term),
		term = gsub("1", "",  term),
		term = toupper(term))

# Plot
ggthemr("earth")
final <-
	ggplot(plot_object, aes(estimate, term, fill = term)) +
		geom_errorbar(
			aes(xmin = conf.low, xmax = conf.high),
			color = "#3f3f3f",
			width = 0.3,
			linewidth = 0.2) +
		geom_vline(
			xintercept = 0,
			color = "#4b4b4b",
			linetype = "dashed",
			linewidth = 0.2) +
		geom_col(width = 0.7) +
		facet_wrap(~ phenotype, nrow = 2) +
		labs(x = "\u03B2 (95% CI)", y = "", fill = "") +
		theme_publish() +
		theme(legend.position = "none")

# Save plot
ggsave(
	"Fig2_panelG.png", device = "png", units = "mm",
	width = 80, height = 200, bg = "white")
