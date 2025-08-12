pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, broom)
options(scipen = 999) # disable scientific notation

# Make a table in this...

# Dataset
data <- readRDS("/home/santorolab/Desktop/cassia/pendrive_cass_BK/cass_BHRC_28042025_ARTICLE.RDS")

# Proband data
database <-
  data$proband_data %>%
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS), by = "IID") %>%
  select(-FID)

# CM-PRS
# https://www.repository.cam.ac.uk/items/0cf1dc7c-0842-434f-b425-3bc5a9f33e05
mi <-
	readRDS("/home/santorolab/Desktop/cassia/pendrive_cass_BK/0_external_files/Lucas_MINI_BHRCS.rds") %>%
	select(ident, IID) %>%
	mutate(ident = as.numeric(ident))

# PRS correction
shapiro.test(all_pcs$PRS)
new_ADHD <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS_CM = new_CM, PRS_ADHD = new_ADHD)

# Family history
parents <- data$family_history

# Working data
test_ADHD <-
	inner_join(database, parents, by = "IID") %>%
	select(PRS_ADHD, starts_with("parent_")) %>%
	rename(
		PRS = 1, Mania = 2, Depression = 3,
		Panic = 4, Psychosis = 5, ADHD = 6,
		Alcohol_abuse = 7, Drug_abuse = 8, Anxiety = 9)

# PRS load: PRS association testing of groups
# based on the parent diagnosis
# plot_object <-
	rbind(
		glm(PRS ~ ., data = test_ADHD, family = "gaussian") %>%
		tidy(conf.int = TRUE)) %>%
	# filter(!term == "(Intercept)") %>%
	# mutate(
	# 	term = gsub("_", " ",  term),
	# 	term = gsub("1", "",  term)) %>%

temp <- report::report_table(glm(PRS ~ ., data = test_ADHD, family = "gaussian"))
rempsyc::nice_table(temp, highlight = TRUE, report = "lm", short = TRUE)

# Plot
ggthemr("earth")
final <-
	ggplot(plot_object, aes(estimate, term, fill = term)) +
		geom_errorbar(
			aes(xmin = conf.low, xmax = conf.high),
			color = "#3f3f3f",
			width = 0.3,
			linewidth = 0.2) +
		geom_vline(
			xintercept = 0,
			color = "#4b4b4b",
			linetype = "dashed",
			linewidth = 0.2) +
		geom_col(width = 0.7) +
		scale_x_continuous(limits = c(-0.25, 0.25)) +
		labs(x = "\u03B2 (95% CI)", y = "", fill = "") +
		theme_publish() +
		theme(
			legend.position = "none",
			panel.grid.major.x = element_line(linewidth = 0.2))

# Save plot
ggsave(
	"Fig3_panelA.png", device = "png", units = "mm",
	width = 100, height = 80, bg = "white")
