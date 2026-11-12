# Model performance over time (all)
pacman::p_load(dplyr, tidyr, nsROC, PRROC, DescTools, purrr, patchwork, ggplot2, envalysis, ggthemr)

data <- readRDS("C:/Users/Cass/Documents/work/cass_12102026_ARTICLE.rds")

# Family history
hist <-
	data$family_history %>%
	select(IID, any_hist)

# Proband data
database <-
  data$proband_data %>% # latest version
  inner_join(., hist, by = "IID") %>%
 	mutate(across(c(W0, W1, W2), ~ifelse(.x == 2, 1, 0))) # W3 is already as 1/0

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
  observed <- model.response(model.frame(model))
  predicted <- as.numeric(predict(model, type = "response"))

  r2 <- as.numeric(PseudoR2(model, which = "Nagelkerke"))

  auroc <- as.numeric(PRROC::roc.curve(
    scores.class0 = predicted,
    weights.class0 = observed,
    curve = FALSE)$auc)

  aucpr <- as.numeric(PRROC::pr.curve(
    scores.class0 = predicted,
    weights.class0 = observed,
    curve = FALSE)$auc.integral)

  tibble(wave = wave,
         predictor = c("R2", "AUROC", "AUCPR"),
         value = c(r2, auroc, aucpr))}
# -----------------------
nrow(data)
waves <- c("W0", "W1", "W2", "W3")
tabs <- map(datasets, correct_prs)

for_plot <- imap_dfr(
  tabs,
  ~ map_dfr(waves, get_metrics, df = .x, sex_adjust = (.y == "all")) %>%
  mutate(subset = .y)) %>%
  mutate(
		predictor = recode(predictor, "R2" = "R²"),
    wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
    predictor = factor(predictor, levels = c("R²", "AUROC", "AUCPR")),
    subset = factor(subset, levels = c("all", "males", "females")))

# Plot dataset (put into one plotting object)
ggthemr("fresh")

final <-
    ggplot(for_plot, aes(x = wave, y = value * 100, group = predictor, linetype = predictor, shape = predictor)) +
    geom_line(linewidth = 0.7, color = "#828282") +
    geom_point(aes(color = wave), size = 3) +
    scale_color_manual(
        values = c(
            W0 = "#65ADC2",
            W1 = "#233B43",
            W2 = "#E84646",
            W3 = "#9B59B6")) +
    scale_linetype_manual(
        values = c(
            "R²" = "solid",
            "AUROC" = "dashed",
            "AUCPR" = "dotted")) +
    scale_shape_manual(values = c("R²" = 16, "AUROC" = 17, "AUCPR" = 15)) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.20))) +
    scale_y_continuous(
        limits = c(0, 70),
        breaks = seq(0, 70, by = 5),
        labels = \(x) ifelse(x %% 10 == 0, sprintf("%.2f", x), ""),
        expand = expansion(mult = c(0.02, 0.05))) +
    facet_wrap(
        vars(subset),
        nrow = 1,
        labeller = as_labeller(c(
            all = "All samples (N=1 463)",
            males = "Males (N=789)",
            females = "Females (N=674)"))) +
    coord_cartesian(clip = "off") +
    labs(x = "Wave",
         y = "Predictor",
         linetype = "Metric",
         shape = "Metric") +
    theme_publish(base_size = 12) +
    theme(
        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed", linewidth = 0.4),
        axis.line = element_line(linewidth = 0.2),
        strip.background = element_rect(fill = "grey80", color = NA),
        strip.text = element_text(size = 10),
        plot.margin = margin(10, 30, 10, 10),
        legend.position = "bottom") +
    guides(color = "none")

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
