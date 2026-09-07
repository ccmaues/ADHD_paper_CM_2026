pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, patchwork, stringr)
options(scipen = 999) # disable scientific notation

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")
database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS = new_PRS)

# Stratify by sex
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

## We used the original PRS for the risk estimation here
## Higher sample size should be used for quantile
## In this one we just want to see the pattern of the first column
## thus, we use its adjusted PRS
calc_prev <- function(data, n, column_name, wave) {
  df <-
    select(data, all_of(column_name), wave) %>%
    mutate(quantile_n = ntile(!!sym(column_name), n)) %>%
    select(quantile_n, wave)
  p <-
    as.data.frame(cbind(1:n, table(df$quantile_n, df[[wave]]))) %>%
    rename("ntile" = 1, "Control" = 2, "Case" = 3) %>%
    mutate(
      prev = Case / (Control + Case),
      ntile = factor(ntile)) %>%
    select(-Control, -Case)
  return(p)}

## for all datas
p1 <- calc_prev(database, 5, "PRS", "W0") %>% rename(W0 = 2)
p2 <- calc_prev(database, 5, "PRS", "W1") %>% rename(W1 = 2)
p3 <- calc_prev(database, 5, "PRS", "W2") %>% rename(W2 = 2)
p4 <- calc_prev(database, 5, "PRS", "W3") %>% rename(W3 = 2)

## Prevalence calculation (Female)
p1_fem <- calc_prev(females, 5, "PRS", "W0") %>% rename(W0 = 2)
p2_fem <- calc_prev(females, 5, "PRS", "W1") %>% rename(W1 = 2)
p3_fem <- calc_prev(females, 5, "PRS", "W2") %>% rename(W2 = 2)
p4_fem <- calc_prev(females, 5, "PRS", "W3") %>% rename(W3 = 2)

## Prevalence calculation (Male)
p1_male <- calc_prev(males, 5, "PRS", "W0") %>% rename(W0 = 2)
p2_male <- calc_prev(males, 5, "PRS", "W1") %>% rename(W1 = 2)
p3_male <- calc_prev(males, 5, "PRS", "W2") %>% rename(W2 = 2)
p4_male <- calc_prev(males, 5, "PRS", "W3") %>% rename(W3 = 2)

# Delta plot preparation
all_plot_delta <-
	inner_join(p1, p4, by = "ntile") %>%
	mutate(deltaW0W3 = W3 - W0) %>%
	select(1, 4) %>%
	rename(delta = 2)

female_plot_delta <-
	inner_join(p1_fem, p4_fem, by = "ntile") %>%
	mutate(deltaW0W3 = W3 - W0) %>%
	select(1, 4) %>%
	rename(delta = 2)

male_plot_delta <-
	inner_join(p1_male, p4_male, by = "ntile") %>%
	mutate(deltaW0W3 = W3 - W0) %>%
	select(1, 4) %>%
	rename(delta = 2)
# Plot data preparation
for_plot_overall <-
	plyr::join_all(list(p1, p2, p3, p4), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

for_plot_female <-
	plyr::join_all(list(p1_fem, p2_fem, p3_fem, p4_fem), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

for_plot_male <-
	plyr::join_all(list(p1_male, p2_male, p3_male, p4_male), by = "ntile", type = "inner") %>%
	tidyr::pivot_longer(cols = starts_with("W"), values_to = "prevalence", names_to = "wave")

database_long <-
  select(database, W0, W1, W2, W3, PRS) %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

females_long <-
  select(database, gender, W0, W1, W2, W3, PRS) %>%
  filter(gender == "Female") %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

males_long <-
  select(database, gender, W0, W1, W2, W3, PRS) %>%
  filter(gender == "Male") %>%
  mutate(risk = ntile(PRS, 5)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
		risk = factor(risk, levels = c(1, 2, 3, 4, 5)))

new_x_axis <- c("1st", "2nd", "3rd", "4th", "5th")
ggthemr("fresh")

p1 <-
	ggplot(for_plot_overall, aes(ntile, prevalence * 100, color = wave, group = wave)) +
		geom_line(linetype = "solid", linewidth = 1, alpha = 0.5) +
		geom_point(size = 4) +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 10, limits = c(5, 35)) +
		scale_color_manual(
			values = c(
				W0 = "#65ADC2",
				W1 = "#233B43",
				W2 = "#E84646",
				W3 = "#9B59B6")) +
		theme_publish(base_size = 12) +
		labs(
			y = "Prevalence",
			x = "") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
	  panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.3))

p2 <-
	ggplot(for_plot_male, aes(ntile, prevalence * 100, color = wave, group = wave)) +
		geom_line(linetype = "solid", linewidth = 1, alpha = 0.5) +
		geom_point(size = 4) +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 10, limits = c(5, 35)) +
		scale_color_manual(
			values = c(
				W0 = "#65ADC2",
				W1 = "#233B43",
				W2 = "#E84646",
				W3 = "#9B59B6")) +
		theme_publish(base_size = 12) +
		labs(y = "", x = "", color = "") +
    theme(
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
	  panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.3))

p3 <-
	ggplot(for_plot_female, aes(ntile, prevalence * 100, color = wave, group = wave)) +
		geom_line(linetype = "solid", linewidth = 1, alpha = 0.5) +
		geom_point(size = 4) +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 10, limits = c(5, 35)) +
		scale_color_manual(
			values = c(
				W0 = "#65ADC2",
				W1 = "#233B43",
				W2 = "#E84646",
				W3 = "#9B59B6")) +
		theme_publish(base_size = 12) +
		labs(
			y = "",
			x = "") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
	  panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.3))

## Delta plots
p6 <-
	ggplot(all_plot_delta, aes(ntile, delta * 100)) +
		geom_col(fill = "#129990") +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 7, limits = c(0, 20)) +
    labs(x = "", y = "\u0394 Prevalence") +
    theme_publish(base_size = 12)

p7 <-
	ggplot(male_plot_delta, aes(ntile, delta * 100)) +
		geom_col(fill = "#129990") +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 7, limits = c(0, 20)) +
    labs(x = "PRS quintile", y = "") +
    theme_publish(base_size = 12)

p8 <-
	ggplot(female_plot_delta, aes(ntile, delta * 100)) +
		geom_col(fill = "#129990") +
		scale_x_discrete(labels = new_x_axis) +
		scale_y_continuous(n.breaks = 7, limits = c(0, 20)) +
    labs(y = "", x = "") +
    theme_publish(base_size = 12)

final <-
  ((p1 / p6 + plot_layout(heights = c(1, 0.5))) |
   (p2 / p7 + plot_layout(heights = c(1, 0.5))) |
   (p3 / p8 + plot_layout(heights = c(1, 0.5)))) +
  plot_annotation(tag_levels = "A")

ggsave(
  "Fig5_part1.png",
  final,
  device = "png",
  width = 30,
  height = 15,
  units = "cm",
  dpi = 300,
  bg = "white")
