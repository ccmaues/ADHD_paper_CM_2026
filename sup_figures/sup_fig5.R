# AUROC and AUCPR over time
pacman::p_load(dplyr, data.table, ggplot2, ggthemr, envalysis, tidyr, nsROC, PRROC)

data <- readRDS("/media/santorolab/C207-3566/cass_BHRC_28042025_ARTICLE.RDS")
database <- data$proband_data %>% # latest version
	mutate(
		W0 = ifelse(W0 == 2, 1, 0),
		W1 = ifelse(W1 == 2, 1, 0),
		W2 = ifelse(W2 == 2, 1, 0))
females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# PCA
all_pcs <-
  data$PCA_all_samples %>%
  inner_join(., select(database, IID, PRS, gender), by = "IID")

fem_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(females, IID, PRS), by = "IID")

man_pcs <-
  data$PCA_by_sex %>%
  inner_join(., select(males, IID, PRS), by = "IID")

# PRS correction
shapiro.test(all_pcs$PRS)
new_PRS <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4 + gender, family = "gaussian", data = all_pcs))
database <- cbind(select(database, -PRS), PRS = new_PRS)

shapiro.test(fem_pcs$PRS)
new_PRS_fem <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = fem_pcs))
females <- cbind(select(females, -PRS), PRS = new_PRS_fem)

shapiro.test(man_pcs$PRS)
new_PRS_man <- residuals(glm(PRS ~ PC1 + PC2 + PC3 + PC4, family = "gaussian", data = man_pcs))
males <- cbind(select(males, -PRS), PRS = new_PRS_man)

# AUROC (FPR and TPR)
all_W0_AUROC <- gROC(database$PRS, database$W0, pvac.auc = TRUE, side = "auto")
all_W1_AUROC <- gROC(database$PRS, database$W1, pvac.auc = TRUE, side = "auto")
all_W2_AUROC <- gROC(database$PRS, database$W2, pvac.auc = TRUE, side = "auto")
    
fem_W0_AUROC <- gROC(females$PRS, females$W0, pvac.auc = TRUE, side = "auto")
fem_W1_AUROC <- gROC(females$PRS, females$W1, pvac.auc = TRUE, side = "auto")
fem_W2_AUROC <- gROC(females$PRS, females$W2, pvac.auc = TRUE, side = "auto")

man_W0_AUROC <- gROC(males$PRS, males$W0, pvac.auc = TRUE, side = "auto")
man_W1_AUROC <- gROC(males$PRS, males$W1, pvac.auc = TRUE, side = "auto")
man_W2_AUROC <- gROC(males$PRS, males$W2, pvac.auc = TRUE, side = "auto")

# AUCPR (Recal and Precision)
all_W0_AUCPR <- pr.curve(scores.class0 = database$PRS, weights.class0 = database$W0, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
all_W1_AUCPR <- pr.curve(scores.class0 = database$PRS, weights.class0 = database$W1, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
all_W2_AUCPR <- pr.curve(scores.class0 = database$PRS, weights.class0 = database$W2, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)

fem_W0_AUCPR <- pr.curve(scores.class0 = females$PRS, weights.class0 = females$W0, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
fem_W1_AUCPR <- pr.curve(scores.class0 = females$PRS, weights.class0 = females$W1, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
fem_W2_AUCPR <- pr.curve(scores.class0 = females$PRS, weights.class0 = females$W2, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)

man_W0_AUCPR <- pr.curve(scores.class0 = males$PRS, weights.class0 = males$W0, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
man_W1_AUCPR <- pr.curve(scores.class0 = males$PRS, weights.class0 = males$W1, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
man_W2_AUCPR <- pr.curve(scores.class0 = males$PRS, weights.class0 = males$W2, curve = TRUE, sorted = FALSE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)

ap1 <-
  rbind(
		data.frame(FPR = all_W0_AUROC$points.coordinates[, "FPR"], TPR = all_W0_AUROC$points.coordinates[, "TPR"], wave = "W0"),
    data.frame(FPR = all_W1_AUROC$points.coordinates[, "FPR"], TPR = all_W1_AUROC$points.coordinates[, "TPR"], wave = "W1"),
    data.frame(FPR = all_W2_AUROC$points.coordinates[, "FPR"], TPR = all_W2_AUROC$points.coordinates[, "TPR"], wave = "W2"))

ap2 <-
  rbind(
		data.frame(Recall = all_W0_AUCPR$curve[, 1], Precision = all_W0_AUCPR$curve[, 2], wave = "W0"),
		data.frame(Recall = all_W1_AUCPR$curve[, 1], Precision = all_W1_AUCPR$curve[, 2], wave = "W1"),
		data.frame(Recall = all_W2_AUCPR$curve[, 1], Precision = all_W2_AUCPR$curve[, 2], wave = "W2"))

ap3 <-
  rbind(
		data.frame(FPR = man_W0_AUROC$points.coordinates[, "FPR"], TPR = man_W0_AUROC$points.coordinates[, "TPR"], wave = "W0"),
    data.frame(FPR = man_W1_AUROC$points.coordinates[, "FPR"], TPR = man_W1_AUROC$points.coordinates[, "TPR"], wave = "W1"),
    data.frame(FPR = man_W2_AUROC$points.coordinates[, "FPR"], TPR = man_W2_AUROC$points.coordinates[, "TPR"], wave = "W2"))

ap4 <-
  rbind(
		data.frame(Recall = man_W0_AUCPR$curve[, 1], Precision = man_W0_AUCPR$curve[, 2], wave = "W0"),
		data.frame(Recall = man_W1_AUCPR$curve[, 1], Precision = man_W1_AUCPR$curve[, 2], wave = "W1"),
		data.frame(Recall = man_W2_AUCPR$curve[, 1], Precision = man_W2_AUCPR$curve[, 2], wave = "W2"))

ap5 <-
  rbind(
		data.frame(FPR = fem_W0_AUROC$points.coordinates[, "FPR"], TPR = fem_W0_AUROC$points.coordinates[, "TPR"], wave = "W0"),
    data.frame(FPR = fem_W1_AUROC$points.coordinates[, "FPR"], TPR = fem_W1_AUROC$points.coordinates[, "TPR"], wave = "W1"),
    data.frame(FPR = fem_W2_AUROC$points.coordinates[, "FPR"], TPR = fem_W2_AUROC$points.coordinates[, "TPR"], wave = "W2"))

ap6 <-
  rbind(
		data.frame(Recall = fem_W0_AUCPR$curve[, 1], Precision = fem_W0_AUCPR$curve[, 2], wave = "W0"),
		data.frame(Recall = fem_W1_AUCPR$curve[, 1], Precision = fem_W1_AUCPR$curve[, 2], wave = "W1"),
		data.frame(Recall = fem_W2_AUCPR$curve[, 1], Precision = fem_W2_AUCPR$curve[, 2], wave = "W2"))

ggthemr("fresh")

p1 <-
	ggplot(ap1, aes(x = FPR, y = TPR, color = wave)) +
	geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#ff0000") +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "False-Positive Rate", y = "True-Positive Rate") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

p2 <-
	ggplot(ap2, aes(x = Recall, y = Precision, color = wave)) +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "Recall", y = "Precision") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

p3 <-
	ggplot(ap3, aes(x = FPR, y = TPR, color = wave)) +
	geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#ff0000") +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "False-Positive Rate", y = "True-Positive Rate") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

p4 <-
	ggplot(ap4, aes(x = Recall, y = Precision, color = wave)) +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "Recall", y = "Precision") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

p5 <-
	ggplot(ap5, aes(x = FPR, y = TPR, color = wave)) +
	geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#ff0000") +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "False-Positive Rate", y = "True-Positive Rate") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

p6 <-
	ggplot(ap6, aes(x = Recall, y = Precision, color = wave)) +
  geom_path(linewidth = 1, aes(linetype = wave)) +
  scale_linetype_manual(values = c("W0" = "dotted", "W1" = "dashed", "W2" = "solid")) +
  # scale_color_manual(values = c("#000000", "#353535", "#757575"), labels = c("W0", "W1", "W2")) +
  labs(x = "Recall", y = "Precision") +
  theme_publish() +
  theme(
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom")

final <- (p1 + p2) / (p3 + p4) / (p5 + p6) + plot_annotation(tag_levels = 'A')
#scale_color_manual(values = c("#65acc2", "#233b43", "#e84646"), labels = c("W0", "W1", "W2")) +

ggsave("fig5_sup.png",final, device = "png", width = 200, height = 300, units = "mm", dpi = 300, bg = "white")