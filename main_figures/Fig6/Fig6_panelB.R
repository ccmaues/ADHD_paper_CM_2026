pacman::p_load(dplyr, data.table, ggthemr, ggplot2, envalysis, broom, survival, tidyr)
# make HZ col plot  for Overall HZ, 90, 10, 90  w/ family, 10 w/ family, 90 n/family 10 n/family
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

wd <-
	inner_join(survival_data, hist, by = "IID") %>%
	select(-IID)

# the only difference, is that I have put the percentile
# out of the strata function and added the family_history
all <-
  tidy(
    coxph(Surv(time, status) ~ percentile + gender + site + any_hist, data = wd),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Overall")

for_plot_HR <-
  all %>%
  mutate(
    term = recode(
      term,
      "genderMale" = "Gender",
      "siteRS" = "Site",
      "any_hist" = "Family history",
      "percentile90th" = "90th"),
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""),
    estimate = round(estimate, 2),
    CI = paste0(round(conf.low, 2), "—", round(conf.high, 2), stars)) %>%
  filter(term %in% c("Gender", "90th", "Family history")) %>%
  mutate(term = factor(term, levels = c("90th", "Family history", "Gender")))

ggthemr("grape")

final <-
  ggplot(for_plot_HR, aes(term, estimate, fill = term, color = term)) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.5,
      position = position_dodge(width = 1)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(y = conf.high + 0.2, label = stars),
      position = position_dodge(width = 0.7),
      vjust = 0.7,
      size = 7,
      show.legend = FALSE,
      angle = 90) +
    scale_y_continuous(n.breaks = 8) +
    labs(y = "Hazard Ratio", x = "") +
    theme_publish(base_size = 15) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(
				color = "grey90",
				linetype = "dashed",
				linewidth = 0.4),
      axis.line.x = element_line(linewidth = 0.5),
      axis.line.y = element_line(linewidth = 0.5))

ggsave(
  "Fig6_panelB.png",
  final,
  device = "png",
  width = 10,
  height = 10,
  units = "cm",
  dpi = 400,
  bg = "white")
