# Model performance over time (all)
pacman::p_load(dplyr, tidyr, ggthemr, envalysis, ggplot2, nsROC, PRROC, DescTools)

# Article dataset
data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Diagnosis standardtization
database <-
	data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0),
		across(c(W0, W1, W2), as.numeric)) %>%
	inner_join(., data$PCA_all_samples, by = "IID")

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
# shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = all_pcs))

# working dataset
for_prediction <- cbind(select(database, -PRS), PRS = new_PRS, group = "all")

# Prediction estimate per wave
for_plot <- data.frame(
	wave = factor(c(rep("W0", 3), rep("W1", 3), rep("W2", 3)), levels = c("W0", "W1", "W2")),
	predictor = rep(c("R2", "AUROC", "AUCPR"), 3),
	value = c(
		## W0
		PseudoR2(glm(W0 ~ PRS + gender + age_W0, family = "binomial", data = for_prediction), which = "Nagelkerke"),
		gROC(for_prediction$PRS, for_prediction$W0, pvac.auc = TRUE, side = "auto")$auc,
		pr.curve(scores.class0 = for_prediction$PRS, weights.class0 = for_prediction$W0, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)$auc.integral,
		## W1
		PseudoR2(glm(W1 ~ PRS + gender + age_W1, family = "binomial", data = for_prediction), which = "Nagelkerke"),
		gROC(for_prediction$PRS, for_prediction$W1, pvac.auc = TRUE, side = "auto")$auc,
		pr.curve(scores.class0 = for_prediction$PRS, weights.class0 = for_prediction$W1, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)$auc.integral,
		## W2
		PseudoR2(glm(W2 ~ PRS + gender + age_W2, family = "binomial", data = for_prediction), which = "Nagelkerke"),
		gROC(for_prediction$PRS, for_prediction$W2, pvac.auc = TRUE, side = "auto")$auc,
		pr.curve(scores.class0 = for_prediction$PRS, weights.class0 = for_prediction$W2, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)$auc.integral)) %>%
	mutate(value = value * 100)

# Plot dataset
ggthemr("fresh")
p <-
	ggplot(for_plot, aes(x = wave, y = value, group = 1, color = wave)) +
		geom_line(size = 1.2, color = "#e0e0e0") +
		geom_point(size = 1.5) +
		geom_text(
			aes(label = sprintf("%.2f", value)),
			vjust = -0.8,
			hjust = -0.2,
			size = 3.5,
			fontface = "bold") +
		facet_wrap(~predictor, scales = "free_y", ncol = 1) +
		scale_x_discrete(
			expand = expansion(mult = c(0.05, 0.20))) +
		scale_y_continuous(
			breaks = function(x) seq(min(x), max(x), length.out = 4),
			labels = \(x) sprintf("%.2f", x),
			expand = expansion(mult = c(0.08, 0.15))) +
		coord_cartesian(clip = "off") +
		labs(x = "Wave", y = "Predictor") +
		theme_publish(base_size = 12) +
		theme(
			panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.2),
				axis.line = element_line(linewidth = 0.2),
				plot.title = element_text(face = "bold", size = 14),
				plot.subtitle = element_text(size = 11, color = "grey40"),
				plot.margin = margin(10, 30, 10, 10),
				legend.position = "none")

# save panel A file
ggsave(
	"fig3_panelA.png",
	p,
	device = "png",
	units = "cm",
	width = 10,
	height = 10,
	dpi = 400,
	bg = "white")
