pacman::p_load(dplyr, data.table, ggplot2, envalysis, tidyr, patchwork, ggsurvfit, survminer, png, survival)

source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + gender + site, data = wd)
fit <- survfit(cox)

# Get event prob at 12 yr
strata_names <- names(fit$strata)
strata_sizes <- fit$strata
strata_ends <- cumsum(strata_sizes)
strata_starts <- c(1, head(strata_ends + 1, -1))
res <- lapply(seq_along(strata_names), function(i) {
  times  <- fit$time[strata_starts[i]:strata_ends[i]]
  surv   <- fit$surv[strata_starts[i]:strata_ends[i]]
  event_prob <- 1 - surv
  est <- approx(times, event_prob, xout = 12, method = "linear", rule = 2, f = 0)$y
  data.frame(strata = strata_names[i], time = 12, event_prob_percent = est * 100)})

do.call(rbind, res)

# Event probability
# all data
surv <-
  ggsurvplot(
    fit,                                     # objeto com a função
    data = wd,                     # objeto criador da função
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
  scale_y_continuous(n.breaks = 8, limits = c(0, 0.18), labels = function(y) sprintf("%.2f", y)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 5, 10, 12, 15, 20, 25)) +
  geom_segment(aes(x = 12, xend = 12, y = 0, yend = 0.092), color = "black", linetype = "solid", size = 0.2) +
  geom_point(aes(x = 12, y = 0.09), color = "#670D2F", size = 1) +
  geom_point(aes(x = 12, y = 0.059), color = "#129990", size = 1) +
  geom_text(aes(x = 9, y = 0.10, label = "9.27%"), color = "#670D2F", size = 2.5, hjust = -0.1) +
  geom_text(aes(x = 12.5, y = 0.05, label = "5.78%"), color = "#129990", size = 2.5, hjust = -0.1) +
  labs(y = "Cumulative event probability", x = "Age (yr)") +
  theme_publish(base_size = 7) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(linetype = "dashed", color = "#c1c1c1", size = 0.2))

ggsave(
  "Fig6_panelA.png", p1, device = "png",
  width = 100, height = 60, units = c("mm"),
  dpi = 300, bg = "white")
