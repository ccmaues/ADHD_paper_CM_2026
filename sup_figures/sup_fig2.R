pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, patchwork)

val_10 <- fread("/media/santorolab/C207-3566/PCA_files_cass/cass_final_PCA/all_samples_PCA.eigenval")

data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
pca_all <- data$PCA_all_samples # latest version

for_plot <- inner_join(pca_all, select(data$proband_data, IID, gender, site), by = "IID")

var_exp10 <- val_10 / sum(val_10)

scree_data <-
	rbind(data.frame(PC = 1:10, var_exp = var_exp10)) %>%
	rename(PC = 1, var_exp = 2) %>%
	mutate(var_exp_pct = var_exp * 100)

ggthemr("grape")

# Panel A: PC1 x PC2
fig1 <-
	ggplot(for_plot, aes(PC1, PC2, shape = gender, color = site)) +
	geom_point(size = 3, alpha = 0.5) +
  labs(y = "PC2", x = "PC1", shape = "", color = "") +
  theme_publish() +
	theme(legend.position = "top")

# Panel B: PC1 x PC3
fig2 <-
	ggplot(for_plot, aes(PC1, PC3, shape = gender, color = site)) +
	geom_point(size = 3, alpha = 0.5) +
  labs(y = "PC3", x = "PC1") +
  theme_publish() +
	theme(legend.position = "none")

# Panel C: PC1 x PC4
fig3 <-
	ggplot(for_plot, aes(PC1, PC4, shape = gender, color = site)) +
	geom_point(size = 3, alpha = 0.5) +
  labs(y = "PC4", x = "PC1") +
  theme_publish() +
	theme(legend.position = "none")

# Panel D: scree plot for PCs
fig4 <-
	ggplot(scree_data, aes(x = PC, y = var_exp_pct)) +
		geom_col() +
		geom_line(color = "black", linewidth = 1, alpha = 0.5) +
		geom_point(color = "black", size = 2) +
		geom_text(aes(label = paste0(round(var_exp_pct, 2), "%")), vjust = -1.5, color = "black", size = 5, angle = 25) +
		scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
		scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
		labs(
			x = "Principal Component",
			y = "% explained variance") +
		theme_publish() +
		theme(
			legend.position = "none",
			strip.text.x = element_blank(),
			axis.text = element_text(size = 10),
			axis.title = element_text(size = 10))

final <- fig1 / (fig2 + fig3) / fig4 + plot_annotation(tag_levels = 'A')
final

ggsave("fig2_sup.png", final, device = "png", height = 300, width = 200, units = "mm")