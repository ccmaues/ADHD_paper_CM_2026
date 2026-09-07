pacman::p_load(dplyr, data.table, ggthemr, ggplot2, envalysis, broom, survival, tidyr)
# make HZ col plot with stratification for high risk strata 90th>=

source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

wd <-
	inner_join(survival_data, hist, by = "IID") %>%
	select(-IID)

# the only difference, is that I have put the percentile
# out of the strata function and added the family_history
top_risk <-
  tidy(
    coxph(Surv(time, status) ~ gender + site + any_hist, data = filter(wd, percentile == "90th")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "High risk")
low_risk <-
  tidy(
    coxph(Surv(time, status) ~ gender + site + any_hist, data = filter(wd, percentile == "10th")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Low risk")

for_plot_HR <-
  rbind(top_risk, low_risk) %>%
  mutate(
    term = recode(
      term,
      "siteRS" = "Site",
      "genderMale" = "Gender",
      "any_hist" = "Family\nhistory"),
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""),
    data = "all",
    estimate = round(estimate, 2),
    CI = paste0(round(conf.low, 2), "—", round(conf.high, 2), stars)) %>%
  filter(term %in% c("Gender", "Family\nhistory")) %>%
  mutate(term = factor(term, levels = c("Family\nhistory", "Gender")))

ggthemr("grape")

final <-
  ggplot(for_plot_HR, aes(term, estimate, fill = term, color = term)) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.5,
      position = position_dodge(width = 1)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(y = conf.high + 0.5, label = stars),
      position = position_dodge(width = 0.7),
      vjust = 0.7,
      size = 7,
      show.legend = FALSE,
      angle = 90) +
    scale_y_continuous(n.breaks = 5) +
    labs(y = "Hazard Ratio", x = "") +
    theme_publish(base_size = 15) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.4),
      axis.line.x = element_line(linewidth = 0.3),
      axis.line.y = element_line(linewidth = 0.3),
    strip.background = element_rect(fill = "#c4c4c4", linewidth = 0),
    strip.text = element_text(face = "bold", color = "#424141", size = 10)) +
    facet_wrap(~ group, nrow = 1, scales = "free_x")

ggsave(
  "Fig6_panelF.png",
  final,
  device = "png",
  width = 10,
  height = 7,
  units = "cm",
  dpi = 400,
  bg = "white")
