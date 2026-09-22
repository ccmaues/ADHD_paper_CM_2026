# figure 1 and 2 (old ones) into a unique panel
pacman::p_load(data.table, ggplot2, ggthemr, envalysis, tidyverse, broom, patchwork)
options(scipen = 999) # disable scientific notation

# MAKE VERY SIMILAR TO THE https://www.nature.com/articles/s41380-023-02293-8
data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")
database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x))) %>%
	select(IID, gender, W0, W1, W2, W3, starts_with("age_"), PRS)
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
	data$PCA_all_samples %>%
	inner_join(., select(database, IID, PRS, gender), by = "IID")

females_pcs <-
	filter(database, gender == "Female") %>%
	inner_join(., data$PCA_by_sex, by = "IID")

males_pcs <-
	filter(database, gender == "Male") %>%
	inner_join(., data$PCA_by_sex, by = "IID")

# PRS correction per subset
shapiro.test(all_pcs$PRS)

new_PRS <-
	all_pcs %>%
	mutate(
		PRS = residuals(
			glm(PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian", data = .))) %>%
	select(IID, PRS)

database <-
	select(database, -PRS) %>%
	inner_join(., new_PRS, by = "IID")

shapiro.test(females_pcs$PRS)

new_PRS_fem <-
	females_pcs %>%
	mutate(
		PRS = residuals(
			glm(PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian", data = .))) %>%
	select(IID, PRS)

females <-
	filter(database, gender == "Female") %>%
	select(-PRS) %>%
	inner_join(., new_PRS_fem, by = "IID")

shapiro.test(males_pcs$PRS)

new_PRS_man <-
	males_pcs %>%
	mutate(
		PRS = residuals(
			glm(PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian", data = .))) %>%
	select(IID, PRS)

males <-
	filter(database, gender == "Male") %>%
	select(-PRS) %>%
	inner_join(., new_PRS_man, by = "IID")

# plot object
database_long <-
  select(database, W0, W1, W2, W3, PRS, gender) %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")),
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

females_long <-
  select(females, gender, W0, W1, W2, W3, PRS) %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")),
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

males_long <-
  select(males, gender, W0, W1, W2, W3, PRS) %>%
  mutate(risk = ntile(PRS, 100)) %>%
	pivot_longer(
		cols = starts_with("W"),
		names_to = "wave",
		values_to = "diagnosis") %>%
  mutate(
		wave = factor(wave, levels = c("W0", "W1", "W2", "W3")),
		diagnosis = factor(diagnosis, levels = c(0, 1)),
    type1 = case_when(
      risk >= 90 ~ "high",
      risk <= 10 ~ "low",
      TRUE ~ "else"),
    type2 = ifelse(risk <= 10, "low", "else"),
    type3 = ifelse(risk >= 90, "high", "else"),
    type1 = factor(type1, levels = c("low", "else", "high")),
    type2 = factor(type2, levels = c("else", "low")),
    type3 = factor(type3, levels = c("else", "high")))

# Testing if the inclusion of wave into the model changes anything (IT DOES)
# Type 1: 90th x 10th
full_model <- glm(diagnosis ~ type1 * wave + gender, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type1 + gender, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: with significance

# Type 2: 10th x else
full_model <- glm(diagnosis ~ type2 * wave + gender, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type2 + gender, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: with significance

# Type 3: 90th x else
full_model <- glm(diagnosis ~ type3 * wave + gender, family = "binomial", data = database_long)
reduced_model <- glm(diagnosis ~ type3 + gender, family = "binomial", data = database_long)
anova(reduced_model, full_model, test = "LRT")
# note: with significance

# We want to see wave-specific odds ratios of diagnosis (ORs) for each risk group
# so we do stratify this

################### ALL SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
	database_long %>%
	filter(type1 %in% c("high", "low")) %>%
	group_by(wave) %>%
	do(tidy(
		glm(
			diagnosis ~ type1 + gender,
			family = "binomial",
			data = .),
		exponentiate = TRUE,
		conf.int = TRUE)) %>%
	filter(term == "type1high") %>%
	transmute(
		wave,
		group = "high",
		OR = estimate,
		CI_low = conf.low,
		CI_high = conf.high,
		p.value)

# Low risk vs. else stratified by wave
type2_results <-
	database_long %>%
	group_by(wave) %>%
	do(tidy(
		glm(
			diagnosis ~ type2 + gender,
			family = "binomial",
			data = .),
		exponentiate = TRUE,
		conf.int = TRUE)) %>%
	filter(term == "type2low") %>%
	transmute(
		wave,
		group = "low",
		OR = estimate,
		CI_low = conf.low,
		CI_high = conf.high,
		p.value)

# High risk vs. else stratified by wave
type3_results <-
	database_long %>%
	group_by(wave) %>%
	do(tidy(
		glm(
			diagnosis ~ type3 + gender,
			family = "binomial",
			data = .),
		exponentiate = TRUE,
		conf.int = TRUE)) %>%
	filter(term == "type3high") %>%
	transmute(
		wave,
		group = "high",
		OR = estimate,
		CI_low = conf.low,
		CI_high = conf.high,
		p.value)

fp1 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "all"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "all"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "all"))

################### FEMALE SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
  females_long %>%
  filter(type1 %in% c("high", "low")) %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type1, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "low", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# Low risk vs. else stratified by wave
type2_results <-
  females_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type2, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "low"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "low") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# High risk vs. else stratified by wave
type3_results <-
  females_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type3, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

fp2 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "females"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "females"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "females"))

################### MALES SAMPLES #####################
# High vs. Low risk stratified by wave
type1_results <-
  males_long %>%
  filter(type1 %in% c("high", "low")) %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type1, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "low", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# Low risk vs. else stratified by wave
type2_results <-
  males_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type2, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "low"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "low") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

# High risk vs. else stratified by wave
type3_results <-
  males_long %>%
  group_by(wave) %>%
  do(tidy(glm(diagnosis ~ type3, family = "binomial", data = .), exponentiate = TRUE, conf.int = TRUE)) %>%
  mutate(
    group = ifelse(term == "(Intercept)", "else", "high"),
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high) %>%
  filter(group == "high") %>%
  select(wave, group, OR, CI_low, CI_high, p.value)

fp3 <- bind_rows(
  type1_results %>% mutate(comparison = "90th vs. 10th", strata = "males"),
  type2_results %>% mutate(comparison = "Else vs. 10th", strata = "males"),
  type3_results %>% mutate(comparison = "90th vs. Else", strata = "males"))

#### Plotting
ggthemr("fresh")
for_plot <-
	rbind(fp1, fp2, fp3) %>%
	mutate(
		strata = factor(strata, levels = c("all", "males", "females")),
		comparison = factor(
			comparison,
			levels = c(
				"90th vs. 10th",
				"Else vs. 10th",
				"90th vs. Else")),
		signif = case_when(
			p.value < 0.001 ~ "***",
			p.value < 0.01 ~ "**",
			p.value < 0.05 ~ "*",
			TRUE ~ ""),
		star_y = CI_high + 0.25)

## odds ratio over wave plots
# https://www.epirhandbook.com/en/new_pages/regression.html#forest-plot
dodge <- position_dodge(width = 0.7)

final <-
	ggplot(
		for_plot,
		aes(
			y = OR,
			x = comparison,
			color = wave,
			group = wave)) +

	geom_hline(
		yintercept = 1,
		linetype = "dashed",
		color = "black",
		linewidth = 0.2,
		alpha = 0.5) +

	geom_errorbar(
		aes(
			ymin = CI_low,
			ymax = CI_high),
		position = dodge,
		width = 0.3,
		linewidth = 0.5,
		alpha = 0.6) +

	geom_point(
		position = dodge,
		size = 2) +

  geom_text(
    aes(
      y = star_y,
      label = signif),
    position = dodge,
    angle = 90,
    vjust = 0.8,
    hjust = 0,
    size = 5,
    show.legend = FALSE) +

	labs(
		color = "",
		x = "",
		y = "Odds ratio") +

	scale_color_manual(
		values = c(
			W0 = "#65ADC2",
			W1 = "#233B43",
			W2 = "#E84646",
			W3 = "#9B59B6")) +

	theme_publish(base_size = 10) +
	theme(legend.position = "top") +

	facet_wrap(
		~ strata,
		labeller = as_labeller(
			c(
				"all" = "A",
				"males" = "B",
				"females" = "C")))

ggsave(
  "Figure_4.png",
  device = "png",
  width = 20,
  height = 7,
  units = "cm",
  dpi = 400,
  bg = "white")
