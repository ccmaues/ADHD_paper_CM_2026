# # Odds ratio of ADHD diagnosis per wave predicted by the calculated PRS
# # females and males
# # ATTENTION CASSIA FROM FUTURE: we stopped this one
# # because the results and figure 5 points out for maybe a
# # non linear association between PRS and gender for females
# pacman::p_load(dplyr, broom, flextable, gtsummary)

# data <- readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")

# # Family history
# hist <- data$family_history %>%
#   mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0)) %>%
# 	select(IID, any_hist)

# females <-
# 	data$proband_data %>% # latest version
# 	inner_join(., hist, by = "IID") %>%
# 	mutate(
# 		W0 = ifelse(W0 == 2, 1, 0),
# 		W1 = ifelse(W1 == 2, 1, 0),
# 		W2 = ifelse(W2 == 2, 1, 0),
# 		across(c(W0, W1, W2), as.numeric)) %>%
# 	inner_join(., data$PCA_by_sex, by = "IID") %>%
# 	filter(gender == "Female")

# males <-
# 	data$proband_data %>% # latest version
# 	inner_join(., hist, by = "IID") %>%
# 	mutate(
# 		W0 = ifelse(W0 == 2, 1, 0),
# 		W1 = ifelse(W1 == 2, 1, 0),
# 		W2 = ifelse(W2 == 2, 1, 0),
# 		across(c(W0, W1, W2), as.numeric)) %>%
# 	inner_join(., data$PCA_by_sex, by = "IID") %>%
# 	filter(gender == "Male")

# # PCA + PRS
# new_PRS_fem <- residuals(glm(
# 	PRS ~ PC1 + PC2 + PC3 + PC4,
# 	family = "gaussian",
# 	data = females))

# new_PRS_man <- residuals(glm(
# 	PRS ~ PC1 + PC2 + PC3 + PC4,
# 	family = "gaussian",
# 	data = males))

# # Working data
# wd <- list(
# 	fem = cbind(select(females, -PRS), PRS = new_PRS_fem),
# 	man = cbind(select(males, -PRS), PRS = new_PRS_man))

# str(wd)

# # models
# fw0 <- glm(W0 ~ PRS + age_W0 + any_hist + site, family = "binomial", data = wd$fem)
# fw1 <- glm(W1 ~ PRS + age_W1 + any_hist + site, family = "binomial", data = wd$fem)
# fw2 <- glm(W2 ~ PRS + age_W2 + any_hist + site, family = "binomial", data = wd$fem)

# mw0 <- glm(W0 ~ PRS + age_W0 + any_hist + site, family = "binomial", data = wd$man)
# mw1 <- glm(W1 ~ PRS + age_W1 + any_hist + site, family = "binomial", data = wd$man)
# mw2 <- glm(W2 ~ PRS + age_W2 + any_hist + site, family = "binomial", data = wd$man)

# final <- tbl_merge(
#   tbls = list(
#     tbl_regression(fw0, intercept = TRUE),
#     tbl_regression(fw1, intercept = TRUE),
#     tbl_regression(fw2, intercept = TRUE)),
#   tab_spanner = c("**W0**", "**W1**", "**W2**")) %>%
#   as_flex_table()

# save_as_docx(
#   "Supplementary Table S19_part1" = final,
#   path = "sup_tab19_part1.docx")

# final <- tbl_merge(
#   tbls = list(
#     tbl_regression(mw0, intercept = TRUE),
#     tbl_regression(mw1, intercept = TRUE),
#     tbl_regression(mw2, intercept = TRUE)),
#   tab_spanner = c("**W0**", "**W1**", "**W2**")) %>%
#   as_flex_table()

# save_as_docx(
#   "Supplementary Table S19_part2" = final,
#   path = "sup_tab19_part2.docx")
