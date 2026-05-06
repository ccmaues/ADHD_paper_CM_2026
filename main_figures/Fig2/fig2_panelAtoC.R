# Density plot PRS per diagnosis
pacman::p_load(dplyr, tidyr, ggthemr, envalysis, ggplot2, patchwork)

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

# Plotting object
for_plot <- rbind(
		cbind(select(database, -PRS), PRS = new_PRS, group = "all"),
		cbind(select(females, -PRS), PRS = new_PRS_stratified_fem, group = "females"),
		cbind(select(males, -PRS), PRS = new_PRS_stratified_man, group = "males")) %>%
	mutate(group = factor(group, level = c("all", "females", "males"))) %>%
	select(PRS, W0, W1, W2, group)

means_df <- for_plot %>%
  group_by(group, W2) %>%
  summarise(mean_PRS = mean(PRS), .groups = "drop")

ggthemr("fresh") # find out which theme I have used b4
# okay... At WHICH timepoint should I choose? The latest?
# maybe the BEST prediction point will be used for density plot
# And in the supp we can put all density plots

# get all toguether in the same plot so I can use the same scale
# Panel A:plot 1 (All PRS by case'n'control)
ggplot(for_plot, aes(x = PRS, fill = W2)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(
    data = means_df,
    aes(xintercept = mean_PRS, color = W2),
    linetype = "dashed",
		linewidth = 0.3) +
  facet_wrap(~group, ncol = 1, scales = "free_x") +
  labs(
		caption = paste(
			"N = 1,553; F = 692; M = 861",
			"\nPRS ~ PC1 + PC2 + PC3 + PC4"),
		y = "Density") +
	scale_fill_manual(values = c("#B07AA1", "#E08D3C")) +
	scale_color_manual(values = c("#B07AA1", "#E08D3C")) +
	scale_x_continuous(n.breaks = 6) +
  theme_publish(base_size = 10) +
	theme(
		legend.title = element_blank(),
		panel.grid = element_line(size = 0.2))

#"#A69F98"
# save panel A file
ggsave(
	"fig2_panelAtoC.png",
	device = "png",
	units = "cm",
	width = 7,
	height = 20,
	dpi = 400,
	bg = "white")
