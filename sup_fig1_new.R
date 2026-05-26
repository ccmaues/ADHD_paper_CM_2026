pacman::p_load(dplyr, data.table, ggplot2, raincloudplots, ggthemr, envalysis, tidyr, patchwork)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

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
for_plot <- data.frame(
  PRS = c(new_PRS, new_PRS_fem, new_PRS_man),
  dataset =
		rep(c("all", "female", "male"),
		times = c(length(new_PRS),
							length(new_PRS_fem),
							length(new_PRS_man)))) %>%
	mutate(dataset = factor(dataset, levels = c("all", "female", "male")))

for_plot2 <-
	select(database, gender, starts_with("age_")) %>%
	pivot_longer(cols = starts_with("age_"),
							 names_to = "wave",
							 values_to = "age") %>%
	mutate(wave = gsub("age_", "", wave))

ggthemr("grape")

# Panel A: violin plot (all | fem | male)
p1 <-
	ggplot(for_plot, aes(dataset, PRS, fill = dataset)) +
	geom_violin() +
	geom_boxplot(color = "white", outlier.colour = "red", outlier.size = 2) +
	theme_publish()

ggthemr("fresh")

# Panel B: raincloud plot - horizontal (age x wave)
ggplot(for_plot2, aes(x = age, y = wave, fill = wave)) +
  geom_half_violin(side = "l", adjust = 2, trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.alpha = 0, alpha = 0.7) +
  geom_point(position = position_jitter(height = 0.1), alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "Raincloud Plot of Age by Wave",
    	 x = "Age", y = "Wave")
# note: we might have some bias in the age for females and males, but this
# might be more evident in the tables...?
# Final plot export