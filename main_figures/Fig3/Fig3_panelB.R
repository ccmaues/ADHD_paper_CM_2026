# Odds ratio per PRS strata on W0, W1 and W2 (all samples)
pacman::p_load(dplyr, broom, ggplot2, envalysis, ggthemr, patchwork)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
	select(IID, gender, W0, W1, W2, W3, starts_with("age_"))

# PCA
all_pcs <-
  	data$PCA_all_samples %>%
  	select(IID, PC1, PC2, PC3, PC4) %>%
	inner_join(., select(data$proband_data, IID, PRS), by = "IID")

new_PRS <-
	all_pcs %>%
	mutate(
		PRS = residuals(glm(
			PRS ~ PC1 + PC2 + PC3 + PC4,
			family = "gaussian",
			data = .))) %>%
	select(IID, PRS)

# Family history
hist <-
	data$family_history %>%
	mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
	select(IID, any_hist)

# Working data
wd <-
	plyr::join_all(
		list(database, new_PRS, hist),
		by = "IID",
		type = "inner") %>%
	mutate(decile = ntile(PRS, 10)) %>%
	select(!c(PRS, IID))

# Models
models <-
	list(
		w0 = glm(W0 ~ factor(decile) + gender + any_hist + age_W0, family = binomial, data = wd),
		w1 = glm(W1 ~ factor(decile) + gender + any_hist + age_W1, family = binomial, data = wd),
		w2 = glm(W2 ~ factor(decile) + gender + any_hist + age_W2, family = binomial, data = wd),
		w3 = glm(W3 ~ factor(decile) + gender + any_hist + age_W3, family = binomial, data = wd)) %>%
	lapply(function(mod) {
		tidy(mod, exponentiate = TRUE, conf.int = TRUE) %>%
		filter(grepl("decile", term)) %>%
		mutate(
			decile = 2:10,
			sig = p.value < 0.05) %>%
		bind_rows(
			tibble(
				decile = 1,
				estimate = 1,
				conf.low = 1,
				conf.high = 1,
				p.value = NA,
				sig = FALSE)) %>%
		arrange(decile)})

# Case counts per decile -------------------------------------------------
# Plot dataset ----------------------------------------------------------
plot_data <-
    bind_rows(
        W0 = models$w0,
        W1 = models$w1,
        W2 = models$w2,
        W3 = models$w3,
        .id = "wave") %>%
    mutate(wave = factor(wave, levels = c("W0", "W1", "W2", "W3")))

# Case counts per decile ------------------------------------------------
case_counts <-
    bind_rows(
        W0 = wd %>%
            group_by(decile) %>%
            summarise(n_cases = sum(W0 == 1, na.rm = TRUE), .groups = "drop"),
        W1 = wd %>%
            group_by(decile) %>%
            summarise(n_cases = sum(W1 == 1, na.rm = TRUE), .groups = "drop"),
        W2 = wd %>%
            group_by(decile) %>%
            summarise(n_cases = sum(W2 == 1, na.rm = TRUE), .groups = "drop"),
        W3 = wd %>%
            group_by(decile) %>%
            summarise(n_cases = sum(W3 == 1, na.rm = TRUE), .groups = "drop"),
        .id = "wave") %>%
    mutate(wave = factor(wave, levels = c("W0", "W1", "W2", "W3")))

ggthemr("fresh")

# -----------------------
# All dataset
# -----------------------
ylims <- c(
    0,
    max(plot_data$conf.high, na.rm = TRUE))

max_cases <- max(case_counts$n_cases, na.rm = TRUE)

scale_factor <- (ylims[2] * 0.75) / max_cases

p_facet <-
    ggplot(plot_data, aes(x = decile, y = estimate, color = wave)) +
        geom_col(
            data = case_counts,
            aes(x = decile, y = n_cases * scale_factor, fill = wave),
            inherit.aes = FALSE,
            width = 0.72,
            alpha = 0.15) +
        geom_hline(
            yintercept = 1,
            linetype = "dashed",
            color = "grey",
            linewidth = 0.3) +
        geom_line(
            aes(group = wave),
            color = "#4e4e4e",
            linewidth = 0.5) +
        geom_errorbar(
            aes(ymin = conf.low, ymax = conf.high),
            width = 0,
            linewidth = 0.3,
            alpha = 0.6) +
        geom_point(size = 2) +
        geom_point(
            data = filter(plot_data, sig),
            shape = 21,
            size = 3.2,
            stroke = 0.8,
            fill = NA,
            color = "#4e4e4e") +
        coord_cartesian(ylim = ylims) +
        scale_y_continuous(
            n.breaks = 7,
            sec.axis = sec_axis(
                ~ . / scale_factor,
                name = "N cases")) +
        scale_x_continuous(
            breaks = 1:10,
            labels = c(
                "1st", "2nd", "3rd", "4th", "5th",
                "6th", "7th", "8th", "9th", "10th")) +
        scale_color_manual(
            values = c(
                W0 = "#65ADC2",
                W1 = "#233B43",
                W2 = "#E84646",
                W3 = "#9B59B6")) +
        scale_fill_manual(
            values = c(
                W0 = "#65ADC2",
                W1 = "#233B43",
                W2 = "#E84646",
                W3 = "#9B59B6")) +
        facet_grid(wave ~ .) +
        labs(
            x = "PGS risk strata",
            y = "Odds Ratio",
            color = "Wave") +
        theme_publish(base_size = 10) +
        theme(
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            panel.spacing.y = grid::unit(0.5, "lines"),
            legend.position = "bottom") +
        guides(fill = "none")

ggsave(
    "fig3_panelB.png",
    p_facet,
    device = "png",
    units = "cm",
    width = 10,
    height = 17,
    dpi = 400,
    bg = "white")