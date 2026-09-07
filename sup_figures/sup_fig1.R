pacman::p_load(dplyr, data.table, ggplot2, ggdist, tidyquant, ggthemr, envalysis, tidyr, patchwork)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Samples
database <- data$proband_data # latest version
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

pcs_females <-
	data$PCA_by_sex %>%
	inner_join(., select(females, IID, PRS), by = "IID")

pcs_males <-
	data$PCA_by_sex %>%
	inner_join(., select(males, IID, PRS), by = "IID")

# PRS adjustment
# shapiro.test(all_pcs$PRS)
new_PRS <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = all_pcs))

# shapiro.test(pcs_females$PRS)
new_PRS_fem <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = pcs_females))

# shapiro.test(pcs_males$PRS)
new_PRS_man <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = pcs_males))

# Plotting object
for_plot <- list()
for_plot$panelA <- data.frame(
  PRS = c(new_PRS, new_PRS_fem, new_PRS_man),
  dataset =
		rep(c("all", "female", "male"),
		times = c(length(new_PRS),
							length(new_PRS_fem),
							length(new_PRS_man)))) %>%
	mutate(dataset = factor(dataset, levels = c("all", "female", "male")))

for_plot$panelB <-
	select(database, gender, starts_with("age_")) %>%
	pivot_longer(cols = starts_with("age_"),
							 names_to = "wave",
							 values_to = "age") %>%
	mutate(wave = gsub("age_", "", wave))

ggthemr("grape")

# Panel A: violin plot (all | fem | male)
p1 <-
	ggplot(for_plot$panelA, aes(dataset, PRS, fill = dataset)) +
	scale_y_continuous(n.breaks = 7) +
	geom_violin(alpha = 0.5, width = 0.7) +
	geom_boxplot(
		color = "white",
		outlier.colour = "red",
		outlier.size = 2,
		width = 0.1,
		alpha = 0.5) +
	scale_x_discrete(
	labels = c(
		sprintf("All (N=%d)", length(new_PRS)),
		sprintf("Females (N=%d)", length(new_PRS_fem)),
		sprintf("Males (N=%d)", length(new_PRS_man)))) +
	labs(x = "") +
	theme_publish() +
	theme(
		legend.position = "none",
		panel.grid.major.y = element_line(
			linetype = "dashed",
			linewidth = 0.2,
			color = "#cfcfcf"))

ggthemr_reset()
ggthemr("fresh")

# Panel B: raincloud plot - horizontal (age x wave)
p2 <-
	ggplot(for_plot$panelB, aes(x = age, y = wave, fill = wave)) +
	ggdist::stat_halfeye(
		adjust = 1,
		justification = -0.2,
		.width = 0,
		point_colour = NA,
		scale = 0.5,
		alpha = 0.9) +
	geom_boxplot(
		width = 0.1,
		alpha = 0.5,
		outlier.colour = "red") +
	ggdist::stat_dots(
		side = "left",
		justification = 1.05,
		binwidth = 0.15,
		alpha = 0.4,
		dotsize = 0.2) +
	tidyquant::scale_fill_tq() +
	tidyquant::theme_tq() +
	scale_x_continuous(n.breaks = 15) +
	labs(y = "", fill = "", x = "Age (yr)") +
	theme_publish() +
	theme(
		legend.position = "none",
		panel.grid.major.x = element_line(
			color = "#cfcfcf",
			linetype = "dashed",
			linewidth = 0.2))

# Final plot export
final <- p1 / p2 +  plot_annotation(tag_levels = 'A')
ggsave("fig1_sup.png", device = "png", height = 300, width = 200, units = "mm")

