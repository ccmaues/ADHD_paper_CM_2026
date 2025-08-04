pacman::p_load(dplyr, data.table, ggplot2, envalysis, tidyr, patchwork, ggsurvfit, survminer, survival)

# https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf
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
  select(IID, site, W0, W1, W2, age_W0, age_W1, age_W2, percentile, gender, site)

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
  inner_join(., select(database, IID, site, percentile, gender, site), by = "IID") %>%
  filter(gender == "Female") %>%
  select(-IID, -gender) %>%
  data.frame()

cox <- coxph(Surv(time, status) ~ strata(percentile) + site, data = survival_data)
fit <- survfit(cox)

# Event probability
# all data
surv <-
  ggsurvplot(
    fit,                                     # objeto com a função
    data = survival_data,                     # objeto criador da função
    fun = "event",                            # função de transformação da curva de sobrevivência
    xlab = "Age (yr)",                        # titulo do eixo x
    risk.table = FALSE,                       # tabela de risco
    ggtheme = theme_publish(base_size = 7),   # tema
    risk.table.y.text = FALSE,                # usar legenda de linha
    censor.size = 2.2,                        # tamanho do censor
    size = 0.6)                               # tamanho da linha

p1 <-
  surv$plot +
  scale_color_manual(values = c("90th" = "#670D2F", "else" = "#c4c4c470", "10th" = "#129990")) +
  scale_y_continuous(n.breaks = 8, limits = c(0, 0.2)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 5, 10, 12, 15, 20, 25)) +
  geom_segment(aes(x = 12, xend = 12, y = 0, yend = 0.12), color = "black", linetype = "solid", size = 0.2) +
  geom_point(aes(x = 12, y = 0.113), color = "#670D2F", size = 1) +
  geom_point(aes(x = 12, y = 0.12), color = "#129990", size = 1) +
  geom_text(aes(x = 12.5, y = 0.09, label = "11.3%"), color = "#670D2F", size = 2.5, hjust = -0.1) +
  geom_text(aes(x = 9.3, y = 0.12, label = "12%"), color = "#129990", size = 2.5, hjust = -0.1) +
  labs(y = "Cumulative event probability", x = "Age (yr)") +
  theme_publish(base_size = 7) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

ggsave(
  "Fig2_panelC.png", p1, device = "png",
  width = 90, height = 50, units = c("mm"),
  dpi = 300, bg = "white")
