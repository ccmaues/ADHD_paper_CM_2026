# Prediction over time
pacman::p_load(dplyr, data.table, tidyr, DescTools, nsROC, PRROC, envalysis, ggplot2, ggthemr, patchwork)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")
database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

# Family history
hist <- data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# PCA per subset
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

females_pcs <-
	filter(database, gender == "Female") %>%
	inner_join(., data$PCA_by_sex, by = "IID")

males_pcs <-
	filter(database, gender == "Male") %>%
	inner_join(., data$PCA_by_sex, by = "IID")

# PRS correction per subset
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <-
  cbind(select(database, -PRS), PRS = new_PRS) %>%
  mutate(risk = ntile(PRS, 5)) %>%
  inner_join(., hist, by = "IID")

shapiro.test(females_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = females_pcs))
females <-
  cbind(select(females_pcs, -PRS), PRS = new_PRS) %>%
  mutate(risk = ntile(PRS, 5)) %>%
  inner_join(., hist, by = "IID")

shapiro.test(males_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = males_pcs))
males <-
  cbind(select(males_pcs, -PRS), PRS = new_PRS) %>%
  mutate(risk = ntile(PRS, 5)) %>%
  inner_join(., hist, by = "IID")

# Pseudo-R2 longitudinally
# I could do a ANOVA for checking the R2 diffs
# per risk
ggthemr("fresh")

# add familty history
waves <- c("W0", "W1", "W2", "W3")
calc_r2 <- function(data, stratified = FALSE) {
	purrr::map_dfr(waves, \(w)
		purrr::map_dfr(1:5, \(r) {
			df <- filter(data, risk == r)
			formula <-
				if(!stratified)
					as.formula(paste0(w, " ~ PRS + gender + any_hist + age_", w))
				else
					as.formula(paste0(w, " ~ PRS + any_hist + age_", w))
			tibble(
				risk = r,
				wave = w,
				R2 = PseudoR2(
					glm(formula, family = "binomial", data = df),
					which = "Nagelkerke"))}))}

fp1 <- calc_r2(database)
fp2 <- calc_r2(males, stratified = TRUE)
fp3 <- calc_r2(females, stratified = TRUE)

fp1 <- fp1 %>% mutate(risk = factor(risk, levels = 1:5), wave = factor(wave, levels = waves))
fp2 <- fp2 %>% mutate(risk = factor(risk, levels = 1:5), wave = factor(wave, levels = waves))
fp3 <- fp3 %>% mutate(risk = factor(risk, levels = 1:5), wave = factor(wave, levels = waves))

new_x_axis <- c("1st", "2nd", "3rd", "4th", "5th")

p1 <-
  ggplot(fp1, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 8, limits = c(0, 18)) +
    scale_x_discrete(labels = new_x_axis) +
    scale_fill_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    scale_color_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    labs(y = "Pseudo-R²", x = "PRS quintile") +
    theme_publish(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

p2 <-
  ggplot(fp2, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 8, limits = c(0, 18)) +
    scale_x_discrete(labels = new_x_axis) +
    scale_fill_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    scale_color_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    labs(y = "", x = "PRS quintile") +
    theme_publish(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

p3 <-
  ggplot(fp3, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 8, limits = c(0, 18)) +
    scale_x_discrete(labels = new_x_axis) +
    scale_fill_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    scale_color_manual(
      values = c(
        W0 = "#65ADC2",
        W1 = "#233B43",
        W2 = "#E84646",
        W3 = "#9B59B6")) +
    labs(y = "", x = "PRS quintile") +
    theme_publish(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

final <- p1 + p2 + p3 + plot_annotation(tag_levels = "A")

ggsave(
  "Fig5_part2.png",
  final,
  device = "png",
  width = 30,
  height = 7,
  units = "cm",
  dpi = 400,
  bg = "white")
