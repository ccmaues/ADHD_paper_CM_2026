# Model performance over time (all)
pacman::p_load(dplyr, tidyr, nsROC, PRROC, DescTools, purrr, patchwork, ggplot2, envalysis, ggthemr)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Family history
hist <-
	data$family_history %>%
	select(IID, any_hist)

# Proband data
database <-
  data$proband_data %>% # latest version
  inner_join(., hist, by = "IID") %>%
 	mutate(across(c(W0, W1, W2), ~ifelse(.x == 2, 1, 0)))

datasets <- list(
  all = database %>%
    inner_join(data$PCA_all_samples, by = "IID"),
  females = database %>%
    inner_join(data$PCA_by_sex, by = "IID") %>%
    filter(gender == "Female"),
  males = database %>%
    inner_join(data$PCA_by_sex, by = "IID") %>%
    filter(gender == "Male"))

# Functions ----------------------------
# PRS correction
correct_prs <- function(df) {
  new_prs <- residuals(glm(
      PRS ~ PC1 + PC2 + PC3 + PC4,
      family = "gaussian",
      data = df))
  cbind(select(df, -PRS), PRS = new_prs)}

# evaluation
get_metrics <- function(wave, df, sex_adjust = TRUE) {
  age_var <- paste0("age_", wave)
  formula <- if (sex_adjust) {
    as.formula(paste0(wave, " ~ PRS + any_hist + gender + ", age_var))
  } else {
    as.formula(paste0(wave, " ~ PRS + any_hist + ", age_var))}
  model <- glm(formula, family = "binomial", data = df)
  r2 <- as.numeric(PseudoR2(model, which = "Nagelkerke"))
  auroc <- as.numeric(gROC(X = df$PRS, D = df[[wave]], pvac.auc = TRUE, side = "auto")$auc)
  aucpr <- as.numeric(
    pr.curve(scores.class0 = df$PRS, weights.class0 = df[[wave]],
      curve = TRUE, sorted = FALSE, max.compute = TRUE,
      min.compute = TRUE, rand.compute = TRUE)$auc.integral)
  tibble(wave = wave, predictor = c("R2", "AUROC", "AUCPR"), value = c(r2, auroc, aucpr))}
# -----------------------

waves <- c("W0", "W1", "W2", "W3")
tabs <- map(datasets, correct_prs)

for_plot <- imap_dfr(
  tabs,
  ~ map_dfr(waves, get_metrics, df = .x, sex_adjust = (.y == "all")) %>%
  mutate(subset = .y)) %>%
  mutate(
		predictor = recode(predictor, "R2" = "R²"),
    wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
    predictor = factor(predictor, levels = c("R²", "AUROC", "AUCPR")))

# Plot dataset
ggthemr("fresh")

p1 <-
	filter(for_plot, subset == "all") %>%
	ggplot(aes(x = wave, y = value * 100, group = 1, color = wave)) +
		geom_line(size = 1.2, color = "#c7c7c7") +
		geom_point(size = 3) +
	scale_color_manual(
		values = c(
			W0 = "#65ADC2",
			W1 = "#233B43",
			W2 = "#E84646",
			W3 = "#9B59B6")) +
		facet_wrap(~predictor, scales = "free_y", ncol = 1) +
		scale_x_discrete(
			expand = expansion(mult = c(0.05, 0.20))) +
		scale_y_continuous(
			breaks = function(x) seq(min(x), max(x), length.out = 4),
			labels = \(x) sprintf("%.2f", x),
			expand = expansion(mult = c(0.08, 0.15))) +
		coord_cartesian(clip = "off") +
		labs(x = "", y = "Predictor") +
		theme_publish(base_size = 12) +
		theme(
			panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.4),
				axis.line = element_line(linewidth = 0.2),
				plot.title = element_text(face = "bold", size = 14),
				plot.subtitle = element_text(size = 11, color = "grey40"),
				plot.margin = margin(10, 30, 10, 10),
				legend.position = "none")

p2 <-
	filter(for_plot, subset == "males") %>%
	ggplot(aes(x = wave, y = value * 100, group = 1, color = wave)) +
		geom_line(size = 1.2, color = "#c7c7c7") +
		geom_point(size = 3) +
		scale_color_manual(
			values = c(
				W0 = "#65ADC2",
				W1 = "#233B43",
				W2 = "#E84646",
				W3 = "#9B59B6")) +
		facet_wrap(~predictor, scales = "free_y", ncol = 1) +
		scale_x_discrete(
			expand = expansion(mult = c(0.05, 0.20))) +
		scale_y_continuous(
			breaks = function(x) seq(min(x), max(x), length.out = 4),
			labels = \(x) sprintf("%.2f", x),
			expand = expansion(mult = c(0.08, 0.15))) +
		coord_cartesian(clip = "off") +
		labs(x = "Wave", y = "") +
		theme_publish(base_size = 12) +
		theme(
			panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.4),
				axis.line = element_line(linewidth = 0.2),
				plot.title = element_text(face = "bold", size = 14),
				plot.subtitle = element_text(size = 11, color = "grey40"),
				plot.margin = margin(10, 30, 10, 10),
				legend.position = "none")
p3 <-
	filter(for_plot, subset == "females") %>%
	ggplot(aes(x = wave, y = value * 100, group = 1, color = wave)) +
		geom_line(size = 1.2, color = "#c7c7c7") +
		geom_point(size = 3) +
		scale_color_manual(
			values = c(
				W0 = "#65ADC2",
				W1 = "#233B43",
				W2 = "#E84646",
				W3 = "#9B59B6")) +
		facet_wrap(~predictor, scales = "free_y", ncol = 1) +
		scale_x_discrete(
			expand = expansion(mult = c(0.05, 0.20))) +
		scale_y_continuous(
			breaks = function(x) seq(min(x), max(x), length.out = 4),
			labels = \(x) sprintf("%.2f", x),
			expand = expansion(mult = c(0.08, 0.15))) +
		coord_cartesian(clip = "off") +
		labs(x = "", y = "") +
		theme_publish(base_size = 12) +
		theme(
			panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.4),
				axis.line = element_line(linewidth = 0.2),
				plot.title = element_text(face = "bold", size = 14),
				plot.subtitle = element_text(size = 11, color = "grey40"),
				plot.margin = margin(10, 30, 10, 10),
				legend.position = "none")

final <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')

# save panel A file
ggsave(
	"fig2.png",
	final,
	device = "png",
	units = "cm",
	width = 30,
	height = 12,
	dpi = 400,
	bg = "white")
