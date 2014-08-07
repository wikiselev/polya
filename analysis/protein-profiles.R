source("func.R")
args <- commandArgs(TRUE)
arg <- args[1]
data <- readRDS(paste0(getwd(), "/", "../../data/map/", arg, "/", arg, ".RDS"))
p.data <- prot_data_make_up(data)

p <- prot_plot1(p.data$p.data[gene.id == 1])
t <- p.data$ann.data[gene.id == 1, list(cond, prot, gene.num)]
out.plot <- paste0(out.folder, arg, "-", "ALL", "-", "1", ".pdf")
out.ann <- paste0(out.folder, arg, "-", "ALL", "-", "1", ".txt")
pdf(out.plot, w = 14, h = 7)
print(p)
dev.off()
write.table(t, file = out.ann, row.names = F, quote = F, sep = "\t")

for(us in c("us_un", "un_us", "in_in")) {
	for(con in unique(p.data$p.data[,cond])) {
		p <-
		prot_plot2(p.data$p.data[cond == con & gene.id == 2 & usage.let == us])
		t <- p.data$ann.data[cond == con & gene.id == 2 & usage.let == us,
						list(prot, gene.num)]
		out.plot <- paste0(out.folder, arg, "-", con, "-", us, ".pdf")
		out.ann <- paste0(out.folder, arg, "-", con, "-", us, ".txt")
		pdf(out.plot, w = 14, h = 7)
		print(p)
		dev.off()
		write.table(t, file = out.ann, row.names = F, quote = F, sep = "\t")
	}
}

for(us in c("us_un", "un_us", "in_in")) {
	for(c.id in 1:2) {
		p <-
		prot_plot1(p.data$p.data[gene.id == 2 & css.id == c.id & usage.let == us])
		t <- p.data$ann.data[gene.id == 2 & css.id == c.id & usage.let == us,
						list(cond, prot, gene.num)]
		out.plot <- paste0(out.folder, arg, "-", "ALL", "-", c.id, "-", us, ".pdf")
		out.ann <- paste0(out.folder, arg, "-", "ALL", "-", c.id, "-", us, ".txt")
		pdf(out.plot, w = 14, h = 7)
		print(p)
		dev.off()
		write.table(t, file = out.ann, row.names = F, quote = F, sep = "\t")
	}
}




# ---
# ## Protein heatmaps - PTB

# ```{r, echo=FALSE}

# if (arg == "a-seq") {
# 	heatmap.cond <- "ctrl"
# }

# if (arg == "atlas") {
# 	heatmap.cond <- "brain"
# }

# test <- data[prot == "PTB" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# setkey(test, gene_id)
# colsum_test <- test[,list(rowsum = sum(score)),by=key(test)]
# test <- test[colsum_test]
# test[,usage:="in"]
# test[usage.pos >= 0.75,usage:="us"]
# test[usage.pos <= 0.25,usage:="un"]
# test[score >= 10, score:=as.integer(10)]
# test <- test[score > 1]
# p <- ggplot(as.data.frame(test), aes(coord, reorder(gene_id, rowsum), fill = log10(score))) + geom_tile() + scale_fill_gradient(low = "green", high = "black") + facet_wrap(gene.id + css.id ~ usage, ncol = 4, scale = "free_y") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TDP43

# ```{r, echo=FALSE}
# test <- data[prot == "TDP43" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TDP43_Hela

# ```{r, echo=FALSE}
# test <- data[prot == "TDP43_Hela" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TDP43_ES

# ```{r, echo=FALSE}
# test <- data[prot == "TDP43_ES" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TDP43_SHSY5Y

# ```{r, echo=FALSE}
# test <- data[prot == "TDP43_SHSY5Y" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TIA1

# ```{r, echo=FALSE}
# test <- data[prot == "TIA1" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TIAL1

# ```{r, echo=FALSE}
# test <- data[prot == "TIAL1" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - TIAL1_Hela

# ```{r, echo=FALSE}
# test <- data[prot == "TIAL1_Hela" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
# ---
# ## Protein heatmaps - U2AF65

# ```{r, echo=FALSE}
# test <- data[prot == "U2AF65" & cond == heatmap.cond]
# ```
# ```{r prot_heatmap, echo=FALSE, include=FALSE}
# ```
# ```{r, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE}
# p
# ```
