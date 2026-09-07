pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, DescTools)

# Explained variance of diagnosis
data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender, W0, W1, W2, W3), by = "IID")

all_pcs$PRS_adj <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = all_pcs))

r2_ref_w3 <-
	PseudoR2(
		glm(W3 ~ PRS_adj, family = "binomial", data = all_pcs),
		which = "Nagelkerke")

all <-
	inner_join(
		select(all_pcs, IID, PC1, PC2, PC3, PC4),
		select(
			database,
			IID, PRS, gender,
			age_W0, age_W1, age_W2, age_W3,
			W0, W1, W2, W3),
		by = "IID")

all_w0 <- data.frame(
	R2 = c(
  	PseudoR2(glm(W0 ~ PRS, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + PC1, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + PC2, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + PC3, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + PC4, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + age_W0, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + gender, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0 + gender, family = "binomial", data = all), which = "Nagelkerke")),
		var = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W0", "Sex", "Full model")) %>%
	mutate(var = factor(var, levels = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W0", "Sex", "Full model")))

all_w1 <- data.frame(
		R2 = c(
	  PseudoR2(glm(W1 ~ PRS, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + PC1, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + PC2, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + PC3, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + PC4, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + age_W1, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + gender, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1 + gender, family = "binomial", data = all), which = "Nagelkerke")),
		var = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W1", "Sex", "Full model")) %>%
	mutate(var = factor(var, levels = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W1", "Sex", "Full model")))

all_w2 <-
	data.frame(
		R2 = c(
	  PseudoR2(glm(W2 ~ PRS, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + PC1, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + PC2, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + PC3, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + PC4, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + age_W2, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + gender, family = "binomial", data = all), which = "Nagelkerke"),
  	PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2 + gender, family = "binomial", data = all), which = "Nagelkerke")),
		var = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W2", "Sex", "Full model")) %>%
	mutate(var = factor(var, levels = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W2", "Sex", "Full model")))

all_w3 <-
	data.frame(
		R2 = c(
			PseudoR2(glm(W3 ~ PRS, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + PC1, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + PC2, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + PC3, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + PC4, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + age_W3, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + gender, family = "binomial", data = all), which = "Nagelkerke"),
			PseudoR2(glm(W3 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W3 + gender, family = "binomial", data = all), which = "Nagelkerke")),
		var = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W3", "Sex", "Full model")) %>%
	mutate(
		var = factor(
			var,
			levels = c("PRS", "PC1", "PC2", "PC3", "PC4", "age_W3", "Sex", "Full model")))

ggthemr("grape")

p1 <-
	ggplot(all_w0, aes(var, R2, fill = var)) +
		geom_col() +
		geom_text(label = paste0(round(all_w0$R2 * 100, 2), "%"), vjust = -1) +
    scale_fill_grey(start = 0.2, end = 0.8) +
		geom_hline(
			yintercept = all_w0$R2[1],
			linetype = "dashed",
			color = "red") +
		scale_y_continuous(limits = c(0, 0.05), n.breaks = 5, expand = expansion(mult = c(0.05, 0.15))) +
		labs(y = "", x = "") +
		theme_publish() +
		theme(
			legend.position = "none",
			axis.title = element_text(size = 10),
			axis.text = element_text(size = 10))

p2 <-
	ggplot(all_w1, aes(var, R2, fill = var)) +
	  geom_col() +
	  geom_text(label = paste0(round(all_w1$R2 * 100, 2), "%"), vjust = -1) +
    scale_fill_grey(start = 0.2, end = 0.8) +
		geom_hline(
			yintercept = all_w1$R2[1],
			linetype = "dashed",
			color = "red") +
	  scale_y_continuous(limits = c(0, 0.05), n.breaks = 5, expand = expansion(mult = c(0.05, 0.15))) +
	  labs(y = "\nExplained Variance [Nagelkerke]\n", x = "") +
	  theme_publish() +
	  theme(
			legend.position = "none",
			axis.title = element_text(size = 10),
			axis.text = element_text(size = 10))

p3 <-
	ggplot(all_w2, aes(var, R2, fill = var)) +
    geom_col() +
    geom_text(label = paste0(round(all_w2$R2 * 100, 2), "%"), vjust = -1) +
    scale_y_continuous(limits = c(0, 0.05), n.breaks = 5, expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_grey(start = 0.2, end = 0.8) +
		geom_hline(
			yintercept = all_w2$R2[1],
			linetype = "dashed",
			color = "red") +
    labs(y = "", x = "") +
    theme_publish() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

p4 <-
	ggplot(all_w3, aes(var, R2, fill = var)) +
	geom_col() +
	geom_text(label = paste0(round(all_w3$R2 * 100, 2), "%"), vjust = -1) +
	scale_fill_grey(start = 0.2, end = 0.8) +
		geom_hline(
			yintercept = all_w3$R2[1],
			linetype = "dashed",
			color = "red") +
	scale_y_continuous(
		limits = c(0, 0.05),
		n.breaks = 5,
		expand = expansion(mult = c(0.05, 0.15))) +
	labs(y = "", x = "") +
	theme_publish() +
	theme(
		legend.position = "none",
		axis.title = element_text(size = 10),
		axis.text = element_text(size = 10))

library(patchwork)

final <- p1 / p2 / p3 / p4 +
	plot_annotation(tag_levels = c("A", "B", "C", "D"))

ggsave("fig6_sup.png",final, device = "png", width = 200, height = 300, units = "mm", dpi = 300, bg = "white")