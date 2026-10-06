# select only those that converted at W3
pacman::p_load(dplyr, ggplot2, envalysis)
data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")
database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
	mutate(w3_conversion = ifelse(W2 == 0 & W3 == 1, TRUE,  FALSE)) %>%
	inner_join(., data$PCA_all_samples, by = "IID")

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

for_plot <-
	cbind(select(database, IID, gender, w3_conversion), PRS = new_PRS$PRS) %>%
	mutate(n.risk = ntile(PRS, 10)) %>%
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	select(gender, n.risk, PC1, PC2, w3_conversion, PRS)

# p1 <-
	ggplot(for_plot, aes(PC1, PC2, shape = gender, color = w3_conversion)) +
		geom_point(size = 4, alpha = 0.5) +
		labs(y = "PC2", x = "PC1", shape = "", color = "") +
		scale_color_manual(values = c("#c6c7c717", "#ff6600")) +
		theme_publish(base_size = 10) +
		theme(
			legend.position = "bottom",
			legend.margin = margin(t = 0, r = 10, b = 0, l = 0),
			panel.grid = element_line(linewidth = 0.2))


ggplot(for_plot, aes(y = PRS, color = w3_conversion)) +
geom_histogram() +
coord_flip()

