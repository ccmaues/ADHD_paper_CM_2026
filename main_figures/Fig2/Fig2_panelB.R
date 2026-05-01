# PRS prediction over time (R2, AUROC, and AUCPR)
pacman::p_load(dplyr, data.table, tidyr, DescTools, nsROC, PRROC, envalysis, ggplot2, ggthemr, patchwork)

# Article dataset
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Diagnosis standardtization
database <-
	data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, "Case", "Control"),
		W1 = ifelse(W1 == 2, "Case", "Control"),
		W2 = ifelse(W2 == 2, "Case", "Control")) %>%
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

# Prediction dataset
prediction <- lapply(
	list(
		cbind(select(database, -PRS), PRS = new_PRS),
		cbind(select(females, -PRS), PRS = new_PRS_stratified_fem),
		cbind(select(males, -PRS), PRS = new_PRS_stratified_man)),
	select, PRS, W0, W1, W2, gender, starts_with("age_W")) %>%
	setNames(., c("all", "fem", "man"))

# Good place to talk WHY with nagelkerke I am using
# age and gender as covariables
# 1. Age: I intent to show ONLY all time diagnosis (atm), so age can be a bias
# 2. Gender: Well... Theres a bunch of info showing that there IS bias
# 3. Why no PCs for R2: correction if applied to PRS (which makes sense),
# not to the diagnosis
# note: for stratified R², gender wont be added to the equation

# -----------------------------
# pseudo R2
# -----------------------------
nagelkerke <- data.frame(
	# all dataset
	all_W0 = PseudoR2(glm(W0 ~ PRS + gender + age_W0, family = "binomial", data = prediction$all), which = "Nagelkerke") * 100,
	all_W1 = PseudoR2(glm(W1 ~ PRS + gender + age_W1, family = "binomial", data = prediction$all), which = "Nagelkerke") * 100,
	all_W2 = PseudoR2(glm(W2 ~ PRS + gender + age_W2, family = "binomial", data = prediction$all), which = "Nagelkerke") * 100,

	# female
	fem_W0 = PseudoR2(glm(W0 ~ PRS + age_W0, family = "binomial", data = prediction$fem), which = "Nagelkerke") * 100,
	fem_W1 = PseudoR2(glm(W1 ~ PRS + age_W1, family = "binomial", data = prediction$fem), which = "Nagelkerke") * 100,
	fem_W2 = PseudoR2(glm(W2 ~ PRS + age_W2, family = "binomial", data = prediction$fem), which = "Nagelkerke") * 100,

	# males
	man_W0 = PseudoR2(glm(W0 ~ PRS + age_W0, family = "binomial", data = prediction$man), which = "Nagelkerke") * 100,
	man_W1 = PseudoR2(glm(W1 ~ PRS + age_W1, family = "binomial", data = prediction$man), which = "Nagelkerke") * 100,
	man_W2 = PseudoR2(glm(W2 ~ PRS + age_W2, family = "binomial", data = prediction$man), which = "Nagelkerke") * 100)

# -----------------------------
# AUROC
# -----------------------------
auroc <- data.frame(
	# all
	all_W0 = gROC(prediction$all$PRS, prediction$all$W0, side = "auto")$auc * 100,
	all_W1 = gROC(prediction$all$PRS, prediction$all$W1, side = "auto")$auc * 100,
	all_W2 = gROC(prediction$all$PRS, prediction$all$W2, side = "auto")$auc * 100,

	# females
	fem_W0 = gROC(prediction$fem$PRS, prediction$fem$W0, side = "auto")$auc * 100,
	fem_W1 = gROC(prediction$fem$PRS, prediction$fem$W1, side = "auto")$auc * 100,
	fem_W2 = gROC(prediction$fem$PRS, prediction$fem$W2, side = "auto")$auc * 100,

	# males
	man_W0 = gROC(prediction$man$PRS, prediction$man$W0, side = "auto")$auc * 100,
	man_W1 = gROC(prediction$man$PRS, prediction$man$W1, side = "auto")$auc * 100,
	man_W2 = gROC(prediction$man$PRS, prediction$man$W2, side = "auto")$auc * 100)

# -----------------------------
# AUCPR
# -----------------------------
aucpr <- data.frame(
	# all
	all_W0 = pr.curve(prediction$all$PRS, prediction$all$W0, sorted = FALSE)$auc.integral * 100,
	all_W1 = pr.curve(prediction$all$PRS, prediction$all$W1, sorted = FALSE)$auc.integral * 100,
	all_W2 = pr.curve(prediction$all$PRS, prediction$all$W2, sorted = FALSE)$auc.integral * 100,

	# females
	fem_W0 = pr.curve(prediction$fem$PRS, prediction$fem$W0, sorted = FALSE)$auc.integral * 100,
	fem_W1 = pr.curve(prediction$fem$PRS, prediction$fem$W1, sorted = FALSE)$auc.integral * 100,
	fem_W2 = pr.curve(prediction$fem$PRS, prediction$fem$W2, sorted = FALSE)$auc.integral * 100,

	# males
	man_W0 = pr.curve(prediction$man$PRS, prediction$man$W0, sorted = FALSE)$auc.integral * 100,
	man_W1 = pr.curve(prediction$man$PRS, prediction$man$W1, sorted = FALSE)$auc.integral * 100,
	man_W2 = pr.curve(prediction$man$PRS, prediction$man$W2, sorted = FALSE)$auc.integral * 100)

# plotting object
for_plot <-
	rbind(R2 = nagelkerke, AUROC = auroc, PR = aucpr) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_time") %>%
	separate(group_time, into = c("group", "time"), sep = "_") %>%
	  pivot_longer(
    cols = c(R2, AUROC, PR),
    names_to = "metric",
    values_to = "value")

ggplot(for_plot, aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.7) +
	facet_grid(metric ~ time, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#a9b8c6", color = "black"),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    legend.position = "none") +
	theme_publish()


ggsave("urgh_facet_kill_me.png", bg = "white")
