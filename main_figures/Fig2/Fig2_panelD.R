pacman::p_load(dplyr, data.table, ggthemr, ggplot2, envalysis, broom, survival, tidyr)
# make HZ col plot  for Overall HZ, 90, 10, 90  w/ family, 10 w/ family, 90 n/family 10 n/family
data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <-
  data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS), by = "IID") %>%
  select(-FID)

# Family history
hist <-
  data$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <-
	cbind(select(database, -PRS), PRS = new_PRS) %>%
  mutate(
    risk = ntile(PRS, 100),
    percentile = case_when(
      risk >= 90 ~ "90th",
      risk <= 10 ~ "10th",
      TRUE ~ "else"),
    percentile = factor(
      percentile,
      levels = c("10th", "else", "90th"))) %>%
	inner_join(., select(hist, IID, any_hist), by = "IID") %>%
  select(IID, site, W0, W1, W2, age_W0, age_W1, age_W2, percentile, gender, site, any_hist)

## keep controls the same
without_entry <-
  filter(database, W2 == 0) %>%
  select(IID, age_W2, W2) %>%
  rename(time = 2, status = 3)
str(without_entry)

## with the first occurance
temp1 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, W0, W1, W2) %>%
  pivot_longer(
    cols = starts_with("W"),
    names_to = "wave",
    values_to = "diagnosis")
str(temp1)

## Age data
temp2 <-
  filter(database, !IID %in% without_entry$IID) %>%
  select(IID, age_W0, age_W1, age_W2) %>%
  pivot_longer(
    cols = starts_with("age_W"),
    names_to = "wave",
    values_to = "age") %>%
  mutate(wave = gsub("age_", "", wave))
str(temp2)

with_entry <-
  inner_join(temp1, temp2, by = c("IID", "wave")) %>%
  filter(diagnosis == 1) %>% # any time diagnosis
  group_by(IID) %>%
  filter(age == min(age)) %>%
  ungroup() %>%
  select(-wave) %>%
  rename(status = 2, time = 3)
str(with_entry)

survival_data <-
  rbind(with_entry, without_entry) %>%
  inner_join(., select(database, IID, site, percentile, gender, site, any_hist), by = "IID") %>%
  select(-IID) %>%
  data.frame()

# the only difference, is that I have put the percentile
# out of the strata function and added the family_history
all <-
  tidy(
    coxph(Surv(time, status) ~ percentile + gender + site + any_hist, data = survival_data),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Overall")
females <-
  tidy(
    coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(survival_data, gender == "Female")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Females")
males <-
  tidy(
    coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(survival_data, gender == "Male")),
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
  "Fig2_panelD.png", final, device = "png",
  width = 100, height = 60, units = "mm",
  dpi = 300, bg = "white")
