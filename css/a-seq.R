source("css_func.R")
data.folder <- "../../data/css/raw/a-seq/"

ctrl <- import.bed(paste0(data.folder, "exp_CS_siRNA_Ctrl.bed"),
	genome = "hg19", asRangedData = FALSE)
values(ctrl)$cond <- "ctrl"
m1 <- import.bed(paste0(data.folder, "exp_CS_siRNA_CFIm68.bed"),
	genome = "hg19", asRangedData = FALSE)
values(m1)$cond <- "CFIm68"
m2 <- import.bed(paste0(data.folder, "exp_CS_siRNA_CstF64.bed"),
	genome = "hg19", asRangedData = FALSE)
values(m2)$cond <- "CstF64"

css <- c(ctrl, m1, m2)

data <- "a-seq"
score.threshold <- 20 # removes 94% of all css
clust.threshold <- 20 # cluster 65% of all css
css.num <- 2 # consider only genes with more than 2 css

source("common_processing.R")

saveRDS(css, paste0(out.folder, "a-seq.RDS"))
