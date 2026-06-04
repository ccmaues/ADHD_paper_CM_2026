# Odds ratio per PRS strata on W0, W1 and W2 (all samples)
pacman::p_load(dplyr, broom, ggplot2, envalysis, ggthemr)

data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0)) %>%
	select(IID, gender, W0, W1, W2)

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  select(IID, PC1, PC2, PC3, PC4) %>%
	inner_join(., select(data$proband_data, IID, PRS), by = "IID")

new_PRS <- residuals(glm(
	PRS ~ PC1 + PC2 + PC3 + PC4,
	family = "gaussian",
	data = all_pcs)) %>%
	as.data.frame() %>%
	cbind(data$proband_data$IID, .) %>%
	rename(IID = 1, PRS = 2)

# Family history
hist <- data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# Working data
wd <-
	plyr::join_all(
		list(database, new_PRS, hist),
		by = "IID", type = "inner") %>%
  	mutate(decile = ntile(PRS, 10)) %>%
	select(!c(PRS, IID))

# Models
models <-
	list(w0 = glm(W0 ~ factor(decile) + gender + any_hist, family = binomial, data = wd),
		 w1 = glm(W1 ~ factor(decile) + gender + any_hist, family = binomial, data = wd),
		 w2 = glm(W2 ~ factor(decile) + gender + any_hist, family = binomial, data = wd)) %>%
  lapply(function(mod) {
    tidy(mod,
         exponentiate = TRUE,
         conf.int = TRUE) %>%
      filter(grepl("decile", term)) %>%
      mutate(decile = 2:10) %>%
      bind_rows(
        tibble(
          decile = 1,
          estimate = 1,
          conf.low = 1,
          conf.high = 1)) %>%
      arrange(decile)})

ggthemr("fresh")
# -----------------------
# All dataset
# -----------------------

ylims <- c(
  min(
    models$w0$conf.low,
    models$w1$conf.low,
    models$w2$conf.low,
    na.rm = TRUE),
  max(
    models$w0$conf.high,
    models$w1$conf.high,
    models$w2$conf.high,
    na.rm = TRUE))

p1 <-
	ggplot(models$w0, aes(x = decile, y = estimate)) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
		geom_line(color = "#4e4e4e") +
		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "#65acc2a1", width = 0, linewidth = 0.3) +
		geom_point(size = 2, color = "#65ADC2") +
		coord_cartesian(ylim = ylims) +
		scale_y_continuous(n.breaks = 7) +
		scale_x_continuous(
			breaks = 1:10,
			labels = c(
				"1st", "2nd", "3rd", "4th", "5th",
				"6th", "7th", "8th", "9th", "10th")) +
		theme_publish(base_size = 10) +
		labs(x = "", y = "") +
		theme_publish()

p2 <-
	ggplot(models$w1, aes(x = decile, y = estimate)) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
		geom_line(color = "#4e4e4e") +
		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "#233b4394", width = 0, linewidth = 0.3) +
		geom_point(size = 2, color = "#233B43") +
		coord_cartesian(ylim = ylims) +
		scale_y_continuous(n.breaks = 7) +
		scale_x_continuous(
			breaks = 1:10,
			labels = c(
				"1st", "2nd", "3rd", "4th", "5th",
				"6th", "7th", "8th", "9th", "10th")) +
		theme_publish(base_size = 10) +
		labs(x = "", y = "Odds Ratio") +
		theme_publish()

p3 <-
	ggplot(models$w2, aes(x = decile, y = estimate)) +
		geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 0.3) +
		geom_line(color = "#4e4e4e") +
		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "#e84646a2", width = 0, linewidth = 0.3) +
		geom_point(size = 2, color = "#E84646") +
		coord_cartesian(ylim = ylims) +
		scale_y_continuous(n.breaks = 7) +
		scale_x_continuous(
			breaks = 1:10,
			labels = c(
				"1st", "2nd", "3rd", "4th", "5th",
				"6th", "7th", "8th", "9th", "10th")) +
		theme_publish(base_size = 10) +
		labs(x = "PRS risk strata", y = "") +
		theme_publish()

final <- p1 / p2 / p3 + plot_annotation(tag_levels = "A")

ggsave(
	"fig2_panelC.png",
	device = "png",
	units = "cm",
	width = 10,
	height = 17,
	dpi = 400,
	bg = "white")
