pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, patchwork)

data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

ggthemr("greyscale")

for_plot <-
	select(database, site, gender, starts_with("age_")) %>%
	pivot_longer(cols = starts_with("age_"), names_to = "wave", values_to = "age") %>%
	mutate(wave = gsub("age_", "", wave))

# Panel A: age per wave
fig1 <-
	ggplot(for_plot, aes(age, wave)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(outliers = TRUE) +
	scale_x_continuous(n.breaks = 10) +
  labs(y = "", x = "Age (yr)") +
  theme_publish()

# Panel B: age by sex and wave
fig2 <-
	ggplot(for_plot, aes(age, wave, fill = gender)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outliers = TRUE) +
	scale_x_continuous(n.breaks = 10) +
  labs(y = "", x = "Age (yr)", fill = "") +
  theme_publish()

# Panel C: age by sex, wave and site
fig3 <-
	ggplot(for_plot, aes(age, wave, fill = site)) +
	stat_boxplot(geom = "errorbar") +
	geom_boxplot(outliers = TRUE) +
	scale_x_continuous(n.breaks = 10) +
	labs(y = "", x = "Age (yr)", fill = "") +
	theme_publish()

final <- fig1 / (fig2 + fig3) + plot_annotation(tag_levels = 'A')
final

ggsave("fig2_sup.png", final, device = "png", height = 300, width = 200, units = "mm")