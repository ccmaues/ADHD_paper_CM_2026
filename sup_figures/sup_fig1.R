pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, patchwork)

data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

fem_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(females, IID, PRS), by = "IID")

man_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(males, IID, PRS), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4 + gender, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS = new_PRS)

shapiro.test(fem_pcs$PRS)
new_PRS_fem <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = fem_pcs))
females <- cbind(select(females, -PRS), PRS = new_PRS_fem)

shapiro.test(man_pcs$PRS)
new_PRS_man <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = man_pcs))
males <- cbind(select(males, -PRS), PRS = new_PRS_man)

ggthemr("greyscale")

# Panel A: all corrected PRS values
fig1 <-
	ggplot(database, aes(sample = PRS)) +
	geom_qq_line(alpha = 0.5, color = "red") +
	geom_qq(size = 1, color = "black") +
	labs(y = "PRS distribution", x = "Theorical distribution") +
	theme_publish()

# Panel B: all corrected PRS stratified by sex
fig2 <-
	ggplot(database, aes(gender, PRS)) +
	geom_boxplot() +
	labs(x = "") +
	theme_publish()

# Panel C: all corrected PRS stratified by site
fig3 <-
	ggplot(database, aes(site, PRS)) +
	geom_boxplot() +
	labs(x = "") +
	theme_publish()

final <- fig1 / (fig2 + fig3) + plot_annotation(tag_levels = 'A')

ggsave("fig1_sup.png", final, device = "png", height = 300, width = 200, units = "mm")