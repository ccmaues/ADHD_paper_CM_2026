pacman::p_load(
	dplyr, data.table, ggplot2, ggthemr,
	envalysis, tidyr, patchwork)

# prevalence per age
data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

## Prepare data for function usage
temp1 <-
	select(
		database,
		IID, gender,
		age_W0, age_W1, age_W2, age_W3) %>%
	pivot_longer(
		cols = starts_with("age_W"),
		names_to = "wave",
		values_to = "age")

temp1$wave <- gsub("age_", "", temp1$wave)

temp2 <-
	select(
		database,
		IID, gender,
		W0, W1, W2, W3) %>%
	pivot_longer(
		cols = c(W0, W1, W2, W3),
		names_to = "wave",
		values_to = "diagnosis")

all <-
	inner_join(
		temp1,
		temp2,
		by = c("IID", "wave", "gender"))

fem <- filter(all, gender == "Female")
man <- filter(all, gender == "Male")

calc_prev_by_age <- function(data, wave, new_column_name) {
	data %>%
		filter(wave == {{wave}}) %>%
		mutate(age = round(age, 0)) %>%
		group_by(diagnosis, age) %>%
		summarise(freq = n(), .groups = "drop") %>%
		spread(age, freq, fill = 0) %>%
		t() %>%
		as.data.frame() %>%
		slice(-1) %>%
		tibble::rownames_to_column(var = "age_diagnosis") %>%
		rename(age = 1, control = 2, case = 3) %>%
		mutate(prev = case / (control + case)) %>%
		filter(control >= 20) %>%
		select(-control) %>%
		rename({{new_column_name}} := 3, N = case)
}

# Overall data
overall_for_plot <-
	rbind(
		cbind(calc_prev_by_age(all, "W3", "prevalence"), wave = "W3"),
		cbind(calc_prev_by_age(all, "W2", "prevalence"), wave = "W2"),
		cbind(calc_prev_by_age(all, "W1", "prevalence"), wave = "W1"),
		cbind(calc_prev_by_age(all, "W0", "prevalence"), wave = "W0")) %>%
	mutate(
		age = as.numeric(age),
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")))

# Female data
female_for_plot <-
	rbind(
		cbind(calc_prev_by_age(fem, "W3", "prevalence"), wave = "W3"),
		cbind(calc_prev_by_age(fem, "W2", "prevalence"), wave = "W2"),
		cbind(calc_prev_by_age(fem, "W1", "prevalence"), wave = "W1"),
		cbind(calc_prev_by_age(fem, "W0", "prevalence"), wave = "W0")) %>%
	mutate(
		age = as.numeric(age),
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")))

# Male data
male_for_plot <-
	rbind(
		cbind(calc_prev_by_age(man, "W3", "prevalence"), wave = "W3"),
		cbind(calc_prev_by_age(man, "W2", "prevalence"), wave = "W2"),
		cbind(calc_prev_by_age(man, "W1", "prevalence"), wave = "W1"),
		cbind(calc_prev_by_age(man, "W0", "prevalence"), wave = "W0")) %>%
	mutate(
		age = as.numeric(age),
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")))

wave_colors <- c(
	"W0" = "#65ADC2",
	"W1" = "#233B43",
	"W2" = "#E84646",
	"W3" = "#9B59B6")

ggthemr("fresh")

p1 <-
	ggplot(
		overall_for_plot,
		aes(
			x = age,
			y = prevalence * 100,
			color = wave,
			group = wave)) +
	geom_line(linewidth = 1, alpha = 0.4) +
	geom_point(aes(size = N)) +
	stat_smooth(
		method = "lm",
		formula = y ~ x,
		geom = "smooth",
		se = FALSE,
		linetype = "dashed",
		linewidth = 1) +
	scale_color_manual(values = wave_colors) +
	scale_x_continuous(n.breaks = 20) +
	scale_y_continuous(
		n.breaks = 8,
		limits = c(0, 30)) +
	labs(y = "\n", x = "") +
	guides(size = guide_legend(nrow = 1)) +
	theme_publish() +
	theme(
		legend.position = "none",
		text = element_text(size = 15),
		axis.text = element_text(size = 15),
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank(),
		axis.line.y = element_line(
			color = "black",
			linewidth = 0.5))

p2 <-
	ggplot(
		male_for_plot,
		aes(
			x = age,
			y = prevalence * 100,
			color = wave,
			group = wave)) +
	geom_line(linewidth = 1, alpha = 0.4) +
	geom_point(aes(size = N)) +
	stat_smooth(
		method = "lm",
		formula = y ~ x,
		geom = "smooth",
		se = FALSE,
		linetype = "dashed",
		linewidth = 1) +
	scale_color_manual(values = wave_colors) +
	scale_x_continuous(n.breaks = 20) +
	scale_y_continuous(
		n.breaks = 8,
		limits = c(0, 30)) +
	labs(y = "Prevalence\n", x = "") +
	theme_publish() +
	theme(
		legend.position = "none",
		text = element_text(size = 15),
		axis.text = element_text(size = 15),
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank(),
		axis.line.y = element_line(
			color = "black",
			linewidth = 0.5))

p3 <-
	ggplot(
		female_for_plot,
		aes(
			x = age,
			y = prevalence * 100,
			color = wave,
			group = wave)) +
	geom_line(linewidth = 1, alpha = 0.4) +
	geom_point(aes(size = N)) +
	stat_smooth(
		method = "lm",
		formula = y ~ x,
		geom = "smooth",
		se = FALSE,
		linetype = "dashed",
		linewidth = 1) +
	scale_color_manual(values = wave_colors) +
	scale_x_continuous(n.breaks = 20) +
	scale_y_continuous(
		n.breaks = 8,
		limits = c(0, 30)) +
	labs(y = "\n", x = "Age (yr)", color = "") +
	theme_publish() +
	theme(
		axis.title = element_text(size = 15),
		axis.text = element_text(size = 15),
		axis.line.y = element_line(
			color = "black",
			linewidth = 0.5))

final <-
	p1 / p2 / p3 +
	plot_annotation(tag_levels = c("A", "B", "C"))

ggsave(
	"fig10_sup.png",
	final,
	device = "png",
	width = 200,
	height = 300,
	units = "mm",
	dpi = 300,
	bg = "white")