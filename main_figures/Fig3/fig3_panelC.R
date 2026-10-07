# Violin plot PRS per diagnosis (W0, W1, W2) all samples
pacman::p_load(dplyr, tidyr, ggthemr, envalysis, ggplot2, patchwork, ggpattern)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Samples & PCs
database <-
	data$proband_data %>% 
	inner_join(., data$PCA_all_samples, by = "IID") %>%
	mutate(across(c(W0, W1, W2, W3), ~ case_when(
		.x %in% c(1, 2) ~ "Case",
		.x == 0 ~ "Control",
		TRUE ~ NA_character_)))

females <-
	filter(database, gender == "Female") %>%
	select(!c(starts_with("PC"), "FID")) %>%
	inner_join(., data$PCA_by_sex, by = "IID")
	
males <- 
	filter(database, gender == "Male") %>%
	select(!c(starts_with("PC"), "FID")) %>%
	inner_join(., data$PCA_by_sex, by = "IID")

# PRS adjustment
# shapiro.test(all_pcs$PRS)
new_PRS <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = database))

# shapiro.test(pcs_females$PRS)
new_PRS_fem <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = females))

# shapiro.test(pcs_males$PRS)
new_PRS_man <-
	residuals(glm(
		PRS ~ PC1 + PC2 + PC3 + PC4,
		family = "gaussian",
		data = males))

# Plotting object
for_plot <-
    data.frame(
        W0 = database$W0,
        W1 = database$W1,
        W2 = database$W2,
        W3 = database$W3,
        PRS = new_PRS) %>%
    pivot_longer(.,
        cols = starts_with("W"),
        names_to = "wave",
        values_to = "status") %>%
    mutate(
        status = factor(status, levels = c("Case", "Control")),
        wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
        group = factor(
            paste(wave, status),
            levels = c(
                "W0 Case", "W0 Control",
                "W1 Case", "W1 Control",
                "W2 Case", "W2 Control",
                "W3 Case", "W3 Control")),
        xpos = c(
            0.86, 1.14,
            1.41, 1.69,
            1.96, 2.24,
            2.51, 2.79)[match(
            paste(wave, status),
            c(
                "W0 Case", "W0 Control",
                "W1 Case", "W1 Control",
                "W2 Case", "W2 Control",
                "W3 Case", "W3 Control"))])

ggthemr("fresh")
# I wanted to make each for wave (the same color scheme i was using)

final <-
    ggplot(for_plot, aes(x = xpos, y = PRS, group = group, fill = wave, pattern = status)) +
        geom_violin_pattern(
            width = 0.24,
            alpha = 0.25,
            color = NA,
            trim = FALSE,
            pattern_angle = 45,
            pattern_density = 0.08,
            pattern_spacing = 0.03,
            pattern_alpha = 0.25,
            pattern_colour = "#4e4e4e") +
        geom_boxplot(
            width = 0.06,
            fill = NA,
            color = "#4e4e4e",
            linewidth = 0.4,
            outlier.shape = 21,
            outlier.fill = "#ff4d4d66",
            outlier.colour = "#ff000099",
            outlier.stroke = 0.4,
            outlier.size = 1) +
        scale_pattern_manual(
            values = c(
                Case = "stripe",
                Control = "none")) +
        scale_fill_manual(
            values = c(
                W0 = "#65ADC2",
                W1 = "#233B43",
                W2 = "#E84646",
                W3 = "#9B59B6")) +
        scale_x_continuous(
            breaks = c(1.00, 1.55, 2.10, 2.65),
            labels = c("W0", "W1", "W2", "W3"),
            expand = expansion(mult = c(0.05, 0.05))) +
        labs(
            x = "",
            y = "PGS",
            fill = "",
            pattern = "") +
        theme_publish(base_size = 14) +
        theme(
            panel.grid.major.y = element_line(
                color = "#cfcfcf",
                linetype = "dashed",
                linewidth = 0.2),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            legend.position = "top") +
    guides(
        fill = guide_legend(
            order = 1,
            override.aes = list(
                pattern = "none")),
        pattern = guide_legend(
            order = 2,
            override.aes = list(
                fill = "grey80",
                colour = NA)))

ggsave(
    "Fig3_panelC.png",
    final,
    device = "png",
    units = "cm",
    width = 13,
    height = 10,
    dpi = 400,
    bg = "white")
