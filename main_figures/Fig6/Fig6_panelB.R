pacman::p_load(dplyr, data.table, ggthemr, ggplot2, envalysis, broom, survival, tidyr)
# make HZ col plot  for Overall HZ, 90, 10, 90  w/ family, 10 w/ family, 90 n/family 10 n/family
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

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
females <-
  tidy(
    coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(wd, gender == "Female")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Females")
males <-
  tidy(
    coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(wd, gender == "Male")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Males")

for_plot_HR <-
  rbind(all, females, males) %>%
  mutate(
    term = recode(
      term,
      "genderMale" = "Gender",
      "siteRS" = "Site",
      "any_hist" = "W/ history",
      "percentile90th" = "90th"),
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""),
    data = "all",
    estimate = round(estimate, 2),
    CI = paste0(round(conf.low, 2), "—", round(conf.high, 2), stars)) %>%
  filter(term %in% c("Gender", "90th", "W/ history")) %>%
  mutate(
    term = factor(term, levels = c("90th", "W/ history", "Gender")),
    group = factor(group, levels = c("Overall", "Males", "Females")))

ggthemr("grape")

final <-
  ggplot(for_plot_HR, aes(term, estimate, fill = term, color = term)) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.5,
      position = position_dodge(width = 1)) +
    geom_col() +
    geom_text(
      aes(
        label = estimate,
        hjust = 0.5,
        vjust = 1.5),
      color = "white",
      size = 3) +
    geom_text(
      aes(y = conf.high + 1, label = stars),
      position = position_dodge(width = 0.7),
      vjust = 0.7,
      size = 5,
      show.legend = FALSE,
      angle = 90) +
    scale_y_continuous(limits = c(0, 6)) +
    # scale_x_discrete(drop = TRUE) +
    labs(y = "Harzard Ratio", x = "") +
    theme_publish(base_family = 7) +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_blank(),
      legend.title = element_blank()) +
    facet_wrap(~ group, nrow = 1, strip.position = "top", scales = "free_x")

ggsave(
  "Fig6_panelD.png", final, device = "png",
  width = 100, height = 60, units = "mm",
  dpi = 300, bg = "white")
