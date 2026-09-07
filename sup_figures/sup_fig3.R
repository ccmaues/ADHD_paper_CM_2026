pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, broom)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")
database <- data$proband_data # latest version

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
model <- glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs)
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

# Combine plots
library(patchwork)
final <- (p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

ggsave("fig3_sup.png", final, device = "png", height = 300, width = 200, units = "mm")