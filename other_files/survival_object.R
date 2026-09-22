pacman::p_load(dplyr, data.table, tidyr)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

database <-
	data$proband_data %>%
	mutate(across(c(W0, W1, W2, W3), ~ ifelse(.x == 2, 1, .x)))

# Family history
hist <-
	data$family_history %>%
	mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

# PCA per subset
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
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4 + gender,
				family = "gaussian",
				data = .))) %>%
	select(IID, PRS)

database <-
	select(database, -PRS) %>%
	inner_join(., new_PRS, by = "IID") %>%
	mutate(
		risk = ntile(PRS, 100),
		percentile = case_when(
			risk >= 90 ~ "90th",
			risk <= 10 ~ "10th",
			TRUE ~ "else"),
		percentile = factor(
			percentile,
			levels = c("10th", "else", "90th")))

shapiro.test(females_pcs$PRS)

new_PRS_fem <-
	females_pcs %>%
	mutate(
		PRS = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = .))) %>%
	select(IID, PRS)

females <-
	filter(database, gender == "Female") %>%
	select(-PRS, -risk, -percentile) %>%
	inner_join(., new_PRS_fem, by = "IID") %>%
	mutate(
		risk = ntile(PRS, 100),
		percentile = case_when(
			risk >= 90 ~ "90th",
			risk <= 10 ~ "10th",
			TRUE ~ "else"),
		percentile = factor(
			percentile,
			levels = c("10th", "else", "90th")))

shapiro.test(males_pcs$PRS)

new_PRS_man <-
	males_pcs %>%
	mutate(
		PRS = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = .))) %>%
	select(IID, PRS)

males <-
	filter(database, gender == "Male") %>%
	select(-PRS, -risk, -percentile) %>%
	inner_join(., new_PRS_man, by = "IID") %>%
	mutate(
		risk = ntile(PRS, 100),
		percentile = case_when(
			risk >= 90 ~ "90th",
			risk <= 10 ~ "10th",
			TRUE ~ "else"),
		percentile = factor(
			percentile,
			levels = c("10th", "else", "90th")))

make_survival <- function(data) {

	# keep controls the same
	without_entry <-
		filter(data, W3 == 0) %>%
		select(IID, age_W3, W3) %>%
		rename(time = 2, status = 3)

	# with the first occurrence
	temp1 <-
		filter(data, !IID %in% without_entry$IID) %>%
		select(IID, W0, W1, W2, W3) %>%
		pivot_longer(
			cols = starts_with("W"),
			names_to = "wave",
			values_to = "diagnosis")

	# Age data
	temp2 <-
		filter(data, !IID %in% without_entry$IID) %>%
		select(IID, age_W0, age_W1, age_W2, age_W3) %>%
		pivot_longer(
			cols = starts_with("age_W"),
			names_to = "wave",
			values_to = "age") %>%
		mutate(wave = gsub("age_", "", wave))

	with_entry <-
		inner_join(temp1, temp2, by = c("IID", "wave")) %>%
		filter(diagnosis == 1) %>%
		group_by(IID) %>%
		filter(age == min(age)) %>%
		ungroup() %>%
		select(-wave) %>%
		rename(status = 2, time = 3)

	survival_data <-
		rbind(with_entry, without_entry) %>%
		inner_join(
			.,
			select(data, IID, site, percentile, gender),
			by = "IID") %>%
		data.frame()

	return(survival_data)}

survival_data <- make_survival(database)
survival_females <- make_survival(females)
survival_males <- make_survival(males)