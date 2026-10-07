# PCA with genetic data (showing risk and gender)
pacman::p_load(dplyr, data.table, tidyr, envalysis, ggplot2, ggthemr)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

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
new_PRS <-
	all_pcs %>%
	mutate(
		PRS = residuals(glm(
			PRS ~ PC1 + PC2 + PC3 + PC4,
			family = "gaussian",
			data = .))) %>%
	select(IID, PRS)

# plotting object
for_plot <-
    database %>%
    select(IID, gender) %>%
    inner_join(new_PRS, by = "IID") %>%
    mutate(
        n.risk = ntile(PRS, 10),
        n.risk = factor(n.risk, levels = 1:10)) %>%
    inner_join(data$PCA_all_samples, by = "IID") %>%
    select(gender, n.risk, PC1, PC2)

risk_colors <-
    colorRampPalette(c(
        "#65ADC2",
        "#E84646",
        "#9B59B6"))(10)

# Panel A: PC1 x PC2
p1 <-
    ggplot(for_plot, aes(PC1, PC2, shape = gender, color = n.risk)) +
        geom_point(size = 1.5, alpha = 0.45) +
        labs(
            y = "PC2",
            x = "PC1",
            shape = "",
            color = "PGS risk strata") +
		guides(color = guide_legend(nrow = 1)) +
        theme_publish(base_size = 10) +
        theme(
            legend.position = "top",
            legend.margin = margin(t = 0, r = 10, b = 0, l = 0),
            panel.grid.major = element_line(
                color = "#cfcfcf",
                linetype = "dashed",
                linewidth = 0.2),
            panel.grid.minor = element_blank()) +
        scale_color_manual(
            values = risk_colors,
            breaks = factor(1:10, levels = 1:10),
            labels = 1:10,
            drop = FALSE)

ggsave(
	"fig3_panelA.png",
	p1,
	device = "png",
	units = "cm",
	width = 20,
	height = 6,
	dpi = 400,
	bg = "white")
