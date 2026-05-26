pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, broom)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")
database <-
  data$proband_data %>% # latest version
	filter(gender == "Female") %>%
	inner_join(., data$PCA_by_sex, by = "IID")

# change for the females subset
val_10 <- fread("D:/cass_HD/DD_CM_backup/PCA_files_cass/cass_final_PCA/all_females_PCA.eigenval")

var_exp10 <- val_10 / sum(val_10)

scree_data <-
	rbind(data.frame(PC = 1:10, var_exp = var_exp10)) %>%
	rename(PC = 1, var_exp = 2) %>%
	mutate(var_exp_pct = var_exp * 100)

# PRS correction
shapiro.test(database$PRS)

model <- glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = database)

model_diagnostics <- augment(model)

model_diagnostics$.fitted <- predict(model, type = "response")

partial_residuals <- residuals(model, type = "partial")

partial_residuals_df <- as.data.frame(partial_residuals)

partial_residuals_df$.fitted <- model_diagnostics$.fitted

ggthemr("grape")

# Panel A: Residuals vs. Fitted Values Plot
p1 <-
	ggplot(model_diagnostics, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted Values", y = "Residuals") +
  theme_publish()

# Panel B: Scale-Location Plot
p2 <- ggplot(model_diagnostics, aes(sample = .std.resid)) +
  geom_qq_line(color = "red") +
  geom_qq(alpha = 0.6, color = "black") +
  labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
  theme_publish()

# Panel C: QQ Plot of Residuals
p3 <-
	ggplot(model_diagnostics, aes(x = .fitted, y = sqrt(abs(.std.resid)))) +
  geom_point(alpha = 0.6, color = "black") +
  geom_smooth(se = FALSE, color = "red") +
  labs(x = "Fitted Values", y = "sqrt(|Standardized Residuals|)") +
  theme_publish()

# Panel D: Cook's Distance Plot
p4 <-
	ggplot(model_diagnostics, aes(x = seq_along(.cooksd), y = .cooksd)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(x = "Observation Index", y = "Cook's Distance") +
  theme_publish()

# Panel E: PC% explained variance
p5 <-
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

# Combine plots
library(patchwork)
final <- (p1 + p2) / (p3 + p4) + p5 + plot_annotation(tag_levels = 'A')
final

ggsave("fig4_alternative1.png", final, device = "png", height = 300, width = 200, units = "mm")