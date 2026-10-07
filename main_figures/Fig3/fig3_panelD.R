pacman::p_load(tidyverse, envalysis)

# Supplementary Table 6 — whole-sample logistic models
table6 <-
    tribble(
        ~variable,        ~wave, ~OR,  ~lower, ~upper,
        "PGS",            "W0",  2.52, 0.93,   6.86,
        "PGS",            "W1",  3.63, 1.43,   9.31,
        "PGS",            "W2",  4.55, 1.83,  11.40,
        "PGS",            "W3",  4.26, 1.81,  10.10,
        "Male vs female", "W0",  1.79, 1.34,   2.41,
        "Male vs female", "W1",  1.88, 1.43,   2.47,
        "Male vs female", "W2",  1.98, 1.52,   2.60,
        "Male vs female", "W3",  1.94, 1.51,   2.49,
        "Age",            "W0",  0.94, 0.87,   1.01,
        "Age",            "W1",  0.91, 0.85,   0.97,
        "Age",            "W2",  0.94, 0.88,   1.00,
        "Age",            "W3",  0.91, 0.86,   0.97,
        "Family history", "W0",  2.18, 1.62,   2.92,
        "Family history", "W1",  2.01, 1.53,   2.66,
        "Family history", "W2",  2.05, 1.56,   2.68,
        "Family history", "W3",  1.78, 1.38,   2.31,
        "RS vs SP",       "W0",  1.48, 1.10,   2.00,
        "RS vs SP",       "W1",  1.41, 1.07,   1.86,
        "RS vs SP",       "W2",  1.39, 1.06,   1.82,
        "RS vs SP",       "W3",  1.29, 1.00,   1.67) %>%
    mutate(
        wave = factor(
            wave,
            levels = c("W0", "W1", "W2", "W3")),
        variable = factor(
            variable,
            levels = rev(c(
                "PGS",
                "Male vs female",
                "Age",
                "Family history",
                "RS vs SP"))),
        significant = case_when(
            variable == "PGS" ~ wave != "W0",
            variable == "Age" ~ wave %in% c("W1", "W3"),
            TRUE ~ TRUE))

# Pastel whiskers
wave_colors <- c(
    W0 = "#78B4C5",
    W1 = "#455A64",
    W2 = "#E98088",
    W3 = "#A574BD")

# Plotting
dodge <- position_dodge(width = 0.6)

pA <-
    ggplot(
        table6,
        aes(
            y = variable,
            x = OR,
            color = wave,
            group = wave)) +

    geom_vline(
        xintercept = 1,
        linetype = "dashed",
        color = "grey60",
        linewidth = 0.4) +

    geom_errorbar(
        aes(
            xmin = lower,
            xmax = upper),
        orientation = "y",
        position = dodge,
        width = 0.3,
        linewidth = 0.6) +

    # Outer ring for significant estimates
    geom_point(
        aes(alpha = if_else(significant, 1, 0)),
        position = dodge,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.6,
        size = 2.8,
        show.legend = FALSE) +

    scale_alpha_identity() +

    # Smaller colored center
    geom_point(
        aes(fill = wave),
        position = dodge,
        shape = 21,
        stroke = 0,
        size = 1.8) +

    scale_color_manual(
        values = wave_colors,
        guide = "none") +

    scale_fill_manual(
        name = "Wave",
        values = c(
            W0 = "#5DAABD",
            W1 = "#263F48",
            W2 = "#EF4444",
            W3 = "#9B59B6")) +

    scale_x_log10(
        breaks = c(0.5, 1, 2, 4, 8, 12),
        labels = c("0.5", "1", "2", "4", "8", "12")) +

    labs(
        x = "Adjusted odds ratio (95% CI)",
        y = NULL) +

    theme_publish(base_size = 12) +
    theme(
        legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank())

ggsave(
    "Fig3_panelD.png",
    plot = pA,
    device = "png",
    units = "cm",
    width = 10,
    height = 13,
    dpi = 400,
    bg = "white")