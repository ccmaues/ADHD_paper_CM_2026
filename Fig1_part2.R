# Prediction over time
pacman::p_load(dplyr, data.table, tidyr, DescTools, nsROC, PRROC, envalysis, ggplot2, ggthemr, patchwork)

data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
  inner_join(., data$PCA_all_samples, by = "IID")

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
# apply te correction AT the glm of R2
# shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))

database <-
  cbind(database, PRS_for_risk = new_PRS) %>%
  mutate(risk = ntile(PRS_for_risk, 5))

# Pseudo-R2 longitudinally
# I could do a ANOVA for checking the R2 diffs
# per risk
ggthemr("fresh")

fp1 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W0, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W0, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W0, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W0, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W0, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W1, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W1, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W1, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W1, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W1, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W2, family = "binomial", data = filter(database, risk == 1)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W2, family = "binomial", data = filter(database, risk == 2)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W2, family = "binomial", data = filter(database, risk == 3)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W2, family = "binomial", data = filter(database, risk == 4)), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + gender + age_W2, family = "binomial", data = filter(database, risk == 5)), which = "Nagelkerke")))

fp2 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 1 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 2 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 3 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 4 & gender == "Male")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 5 & gender == "Male")), which = "Nagelkerke")))

fp3 <-
  data.frame(
    risk = factor(c(1, 2, 3, 4, 5), levels = c(1, 2, 3, 4, 5)),
    wave = c("W0", "W0", "W0", "W0", "W0", "W1", "W1", "W1", "W1", "W1", "W2", "W2", "W2", "W2", "W2"),
    R2 = c(
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W0 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W0, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W1 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W1, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 1 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 2 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 3 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 4 & gender == "Female")), which = "Nagelkerke"),
      PseudoR2(glm(W2 ~ PRS + PC1 + PC2 + PC3 + PC4 + age_W2, family = "binomial", data = filter(database, risk == 5 & gender == "Female")), which = "Nagelkerke")))

new_x_axis <- c("1st", "2nd", "3rd", "4th", "5th")

p1 <-
  ggplot(fp1, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 10, limits = c(0, 17)) +
    scale_x_discrete(labels = new_x_axis) +
    labs(y = "Pseudo-R²", x = "PRS quintile") +
    theme_publish(base_size = 7) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

p2 <-
  ggplot(fp2, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 10, limits = c(0, 17)) +
    scale_x_discrete(labels = new_x_axis) +
    labs(y = "", x = "PRS quintile") +
    theme_publish(base_size = 7) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

p3 <-
  ggplot(fp3, aes(risk, R2 * 100, fill = wave, color = wave, group = wave)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6) +
    scale_y_continuous(n.breaks = 10, limits = c(0, 17)) +
    scale_x_discrete(labels = new_x_axis) +
    labs(y = "", x = "PRS quintile") +
    theme_publish(base_size = 7) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

p1 + p2 + p3 + plot_annotation(tag_levels = "A")

ggsave(
  "Fig1_part2.png", device = "png",
  width = 180, height = 50, units = "mm",
  dpi = 300, bg = "white")
