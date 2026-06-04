# Violin plot PRS per diagnosis (W0, W1, W2) all samples
pacman::p_load(dplyr, tidyr, ggthemr, envalysis, ggplot2, patchwork)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Samples & PCs
database <-
	data$proband_data %>% 
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	mutate(across(starts_with("W"), ~ ifelse(.x == 2, "Case", "Control")))

females <-
	filter(database, gender == "Female") %>%
	select(!c(starts_with("PC"), "FID")) %>%
	inner_join(., data$PCA_by_sex, by = "IID")
	
males <- 
	filter(database, gender == "Male") %>%
	select(!c(starts_with("PC"), "FID")) %>%
	inner_join(., data$PCA_by_sex, by = "IID")

# PRS adjustment
# shapiro.test(all_pcs$PRS)
new_PRS <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = database))

# shapiro.test(pcs_females$PRS)
new_PRS_fem <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = females))

# shapiro.test(pcs_males$PRS)
new_PRS_man <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = males))

# Plotting object
for_plot <-
	data.frame(
		W0 = database$W0,
		W1 = database$W1,
		W2 = database$W2,
		PRS = new_PRS) %>%
	pivot_longer(.,
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "status") %>%
	mutate(group = factor(
      paste(wave, status),
      levels = c(
        "W0 Case", "W0 Control",
        "W1 Case", "W1 Control",
        "W2 Case", "W2 Control")))

ggthemr("fresh")
# I wanted to make each for wave (the same color scheme i was using)

# Violin plot PRS per wave
final <- ggplot(for_plot, aes(x = group, y = PRS, fill = wave)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(
	width = 0.1,
    alpha = 0.5,
	outlier.colour = "red",
	outlier.size = 2) +
  	scale_fill_manual(
		values = c(
			W0 = "#65ADC2",
			W1 = "#233B43",
			W2 = "#E84646")) +
	scale_x_discrete(
		labels = c(
			"W0 Case" = "Case",
			"W0 Control" = "Control",
			"W1 Case" = "Case",
			"W1 Control" = "Control",
			"W2 Case" = "Case",
			"W2 Control" = "Control")) +
	labs(x = "") +
	theme_publish(base_size = 10) +
	theme(
		panel.grid.major.y = element_line(
			color = "#cfcfcf",
			linetype = "dashed",
			linewidth = 0.2))

# "#65ADC2"
# "#233B43"
# "#E84646"
# ggthemr("fresh")$palette$swatch

ggsave(
	"Fig2_panelB.png",
	final,
	device = "png",
	units = "cm",
	width = 10,
	height = 10,
	dpi = 400,
	bg = "white")