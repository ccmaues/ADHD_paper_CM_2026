pacman::p_load(dplyr, ggplot2, envalysis, tidyr, patchwork, ggsurvfit, survminer, survival)

source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

#---------------------------------
# All samples
#---------------------------------
wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + gender + site, data = wd)
fit1 <- survfit(cox)

# Get event prob at 12 yr
s1 <- summary(fit1, times = 12)

bt1 <- sprintf("%.2f", (1 - s1$surv[1]) * 100) # 10th
tp1 <- sprintf("%.2f", (1 - s1$surv[3]) * 100) # 90th

bt1_y <- 1 - s1$surv[1]
tp1_y <- 1 - s1$surv[3]

#---------------------------------
# females
#---------------------------------
cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist, data =  filter(wd, gender == "Female"))
fit2 <- survfit(cox)

s2 <- summary(fit2, times = 12)

bt2 <- sprintf("%.2f", (1 - s2$surv[1]) * 100) # 10th
tp2 <- sprintf("%.2f", (1 - s2$surv[3]) * 100) # 90th

bt2_y <- 1 - s2$surv[1]
tp2_y <- 1 - s2$surv[3]

#---------------------------------
# males
#---------------------------------
cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = filter(wd, gender == "Male"))
fit3 <- survfit(cox)

s3 <- summary(fit3, times = 12)

bt3 <- sprintf("%.2f", (1 - s3$surv[1]) * 100) # 10th
tp3 <- sprintf("%.2f", (1 - s3$surv[3]) * 100) # 90th

bt3_y <- 1 - s3$surv[1]
tp3_y <- 1 - s3$surv[3]

# Panel A
surv <-
  ggsurvplot(
    fit1,                                     # objeto com a função
    data = wd,                                # objeto criador da função
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
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2), limits = c(0, 0.20), labels = function(y) sprintf("%.2f", y)) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 12, 15, 20, 30)) +
  geom_segment(
    aes(x = 12, xend = 12, y = 0, yend = max(bt1_y, tp1_y)),
    color = "black", size = 0.2) +
  geom_point(aes(x = 12, y = tp1_y), color = "#670D2F", size = 1.5) +
  geom_point(aes(x = 12, y = bt1_y), color = "#129990", size = 1.5) +
  geom_text(aes(x = 12.5, y = tp1_y, label = tp1), color = "#670D2F", size = 3, hjust = -0.1) +
  geom_text(aes(x = 12.5, y = bt1_y, label = bt1), color = "#129990", size = 3, hjust = -0.1) +
  labs(y = "", x = "", title = "all") +
  theme_publish(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

# Panel C
surv <-
  ggsurvplot(
    fit2,                                      # objeto com a função
    data = wd,                                # objeto criador da função
    fun = "event",                            # função de transformação da curva de sobrevivência
    xlab = "Age (yr)",                        # titulo do eixo x
    risk.table = FALSE,                       # tabela de risco
    ggtheme = theme_publish(base_size = 10),   # tema
    risk.table.y.text = FALSE,                # usar legenda de linha
    censor.size = 2.2,                        # tamanho do censor
    size = 0.6)                               # tamanho da linha

p2 <-
  surv$plot +
  scale_color_manual(values = c("90th" = "#670D2F", "else" = "#c4c4c470", "10th" = "#129990")) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2), limits = c(0, 0.20), labels = function(y) sprintf("%.2f", y)) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 12, 15, 20, 30)) +
  geom_segment(
    aes(x = 12, xend = 12, y = 0, yend = max(bt2_y, tp2_y)),
    color = "black", size = 0.2) +
  geom_point(aes(x = 12, y = tp2_y), color = "#670D2F", size = 1.5) +
  geom_point(aes(x = 12, y = bt2_y), color = "#129990", size = 1.5) +
  geom_text(aes(x = 12.5, y = tp2_y, label = tp2), color = "#670D2F", size = 3, hjust = -0.1) +
  geom_text(aes(x = 12.5, y = bt2_y, label = bt2), color = "#129990", size = 3, hjust = -0.1) +
  labs(y = "Cumulative event probability", x = "", title = "females") +
  theme_publish(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

# Panel E
surv <-
  ggsurvplot(
    fit3,                                      # objeto com a função
    data = wd,                                # objeto criador da função
    fun = "event",                            # função de transformação da curva de sobrevivência
    xlab = "Age (yr)",                        # titulo do eixo x
    risk.table = FALSE,                       # tabela de risco
    ggtheme = theme_publish(base_size = 7),   # tema
    risk.table.y.text = FALSE,                # usar legenda de linha
    censor.size = 2.2,                        # tamanho do censor
    size = 0.6)                               # tamanho da linha

p3 <-
  surv$plot +
  scale_color_manual(values = c("90th" = "#670D2F", "else" = "#c4c4c470", "10th" = "#129990")) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 12, 15, 20, 30)) +
  geom_segment(
    aes(x = 12, xend = 12, y = 0, yend = max(bt3_y, tp3_y)),
    color = "black", size = 0.2) +
  geom_point(aes(x = 12, y = tp3_y), color = "#670D2F", size = 1.5) +
  geom_point(aes(x = 12, y = bt3_y), color = "#129990", size = 1.5) +
  geom_text(aes(x = 12.5, y = tp3_y, label = tp3), color = "#670D2F", size = 3, hjust = -0.1) +
  geom_text(aes(x = 12.5, y = bt3_y, label = bt3), color = "#129990", size = 3, hjust = -0.1) +
  labs(y = "", x = "Age (yr)", title = "Males") +
  theme_publish(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

# ------------------------------
# Final plot

final <- p1 / p2 / p3

ggsave(
  "Fig6_panelACE.png",
  final,
  device = "png",
  width = 15,
  height = 20,
  units = "cm",
  dpi = 400,
  bg = "white")
