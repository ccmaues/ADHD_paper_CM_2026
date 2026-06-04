# PCA with genetic data (showing risk and gender)
pacman::p_load(dplyr, data.table, tidyr, envalysis, ggplot2, ggthemr)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# Article dataset
database <-
	data$proband_data %>% # latest version
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	mutate(across(c(W0, W1, W2), as.factor))

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

# adjusted PRS values (risk strata info)
new_PRS <- residuals(glm(
	PRS ~ PC1 + PC2 + PC3 + PC4,
	family = "gaussian",
	data = all_pcs))

# plotting object
for_plot <-
	cbind(select(database, IID, gender), PRS = new_PRS) %>%
	mutate(n.risk = ntile(PRS, 10)) %>%
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	select(gender, n.risk, PC1, PC2)

# main theme
ggthemr("grape")

# Panel A: PC1 x PC2
p1 <-
	ggplot(for_plot, aes(PC1, PC2, shape = gender, color = n.risk)) +
		geom_point(size = 1.5, alpha = 0.5) +
		labs(y = "PC2", x = "PC1", shape = "", color = "") +
		theme_publish(base_size = 10) +
		theme(
			legend.position = "bottom",
			legend.margin = margin(t = 0, r = 10, b = 0, l = 0),
			panel.grid = element_line(size = 0.2)) +
		scale_color_viridis_c(
			option = "magma",
			direction = -1,
			breaks = 1:10,
			labels = paste0(1:10))

ggsave(
	"fig2_panelA.png",
	p1,
	device = "png",
	units = "cm",
	width = 9,
	height = 20,
	dpi = 400,
	bg = "white")