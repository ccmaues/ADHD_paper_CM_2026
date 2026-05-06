# PRS prediction over time (R2, AUROC, and AUCPR)
pacman::p_load(dplyr, data.table, tidyr, DescTools, nsROC, PRROC, envalysis, ggplot2, ggthemr, patchwork)

# Article dataset
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Diagnosis standardtization
database <-
	data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	mutate(across(c(W0, W1, W2), as.factor))

females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

sex_stratified_pcs_F <-
	data$PCA_by_sex %>%
	inner_join(., select(females, IID, PRS), by = "IID")

sex_stratified_pcs_M <-
	data$PCA_by_sex %>%
	inner_join(., select(males, IID, PRS), by = "IID")

# PRS correction
# shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(
	PRS ~ PC1 + PC2 + PC3 + PC4,
	family = "gaussian",
	data = all_pcs))

# Why no gender here? The GWAS in which this PRS is derived
# has not shown any statistical relevance regarding gender
# differences. Thus, for values derived from it, not correction.

# shapiro.test(sex_stratified_pcs_F$PRS)
new_PRS_stratified_fem <-	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = sex_stratified_pcs_F))

# shapiro.test(sex_stratified_pcs_M$PRS)
new_PRS_stratified_man <-	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = sex_stratified_pcs_M))

# plotting obj
fixed_prs <- lapply(
	list(
		cbind(select(database, -PRS), PRS = new_PRS),
		cbind(select(females, -PRS), PRS = new_PRS_stratified_fem),
		cbind(select(males, -PRS), PRS = new_PRS_stratified_man)),
	select, PRS, W0, W1, W2, gender) %>%
	lapply(., mutate, n.risk = factor(ntile(PRS, 10))) %>%
	setNames(., c("all", "fem", "man"))

library(broom)
ggthemr("fresh")
# W1 has peak prediction + onset age

# -----------------------
# All dataset
# -----------------------
model <- glm(W1 ~ n.risk + gender, data = fixed_prs$all, family = binomial)
or_df <-
	tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(!term %in% c("(Intercept)", "genderMale")) %>%
  mutate(strata = as.numeric(gsub("n.risk", "", term))) %>%
	bind_rows(
		tibble(
			term = "n.risk1", estimate = 1,
			std.error = 0, statistic = 0,
			p.value = 0, conf.low = 0,
			conf.high = 0, strata = 1)) %>%
  arrange(strata)

ggplot(or_df, aes(x = strata, y = estimate)) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
		geom_line(color = "#4e4e4e") +
		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, linewidth = 0.3) +
		geom_point(size = 2) +
		scale_x_continuous(
			breaks = 1:10,
			labels = c(
				"1st", "2nd", "3rd", "4th", "5th",
				"6th", "7th", "8th", "9th", "10th")) +
		theme_publish(base_size = 10) +
		labs(x = "PRS risk strata", y = "Odds Ratio") +
		theme(panel.grid = element_line(size = 0.2))

ggsave(
	"fig2_panelD.png",
	device = "png",
	units = "cm",
	width = 8,
	height = 7,
	dpi = 400,
	bg = "white")

# # -----------------------
# # female dataset
# # -----------------------
# model <- glm(W1 ~ n.risk, data = fixed_prs$fem, family = binomial)
# or_df <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(q = as.numeric(gsub("n.risk", "", term)))

# p2 <-
# 	ggplot(or_df, aes(x = q, y = estimate)) +
# 		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
# 		geom_line(color = "#4e4e4e") +
# 		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, linewidth = 0.3) +
# 		geom_point(size = 3) +
# 		theme_publish(base_size = 10) +
# 		labs(x = "", y = "Odds ratio") +
# 		theme(
# 			panel.grid = element_line(size = 0.2),
# 			axis.line.x = element_blank(),
# 			axis.title.x = element_blank(),
#     	axis.text.x  = element_blank(),
#     	axis.ticks.x = element_blank())

# # -----------------------
# # male dataset
# # -----------------------
# model <- glm(W1 ~ n.risk, data = fixed_prs$man, family = binomial)
# or_df <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(q = as.numeric(gsub("n.risk", "", term)))

# # weird predic
# p3 <-
# 	ggplot(or_df, aes(x = q, y = estimate)) +
# 		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
# 		geom_line(color = "#4e4e4e") +
# 		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, linewidth = 0.3) +
# 		geom_point(size = 3) +
# 		theme_publish(base_size = 10) +
# 		labs(x = "PRS risk strata", y = "") +
# 		theme(panel.grid = element_line(size = 0.2))

# library(patchwork)

# p1 / p2 / p3

# ggsave(
# 	"fig2_panelB.png",
# 	device = "png",
# 	units = "cm",
# 	width = 7,
# 	height = 20,
# 	dpi = 400,
# 	bg = "white")

# Good place to talk WHY with nagelkerke I am using
# age and gender as covariables
# 1. Age: I intent to show ONLY all time diagnosis (atm), so age can be a bias
# 2. Gender: Well... Theres a bunch of info showing that there IS bias
# 3. Why no PCs for R2: correction if applied to PRS (which makes sense),
# not to the diagnosis
# note: for stratified R², gender wont be added to the equation

