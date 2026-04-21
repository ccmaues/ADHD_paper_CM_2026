# PCA
pacman::p_load(data.table, ggplot2, dplyr)
system2("plink --bfile XXX --pca 10 --out XXX")# REMEMBER DO FUCKING PRUNE
