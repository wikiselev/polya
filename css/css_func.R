# load my general functions
source(paste0(getwd(), "/", "../resources/ann.R"))
out.folder <- paste0(getwd(), "/", "../../data/css/processed/")
report.folder <- "report/"

ens_overlap <- function(d, ens) {
	overlaps <- as.data.frame(findOverlaps(d, ens))
	d <- d[overlaps[ , 1]]
	values(d)$gene_id <- values(ens[overlaps[ , 2]])$gene_id
	# d <- d[ , 2:4]
	return(d)
}

# a function for processing raw css A-seq coordinates
css_processing <- function(css, clust.threshold) {
	css <- css[values(css)$score >= score.threshold]
	ensembl.genes <-
	readRDS("../resources/hsapiens_prot_cod_trs_biomaRt_granges_28Nov2013.RDS")
	css <- ens_overlap(css, ensembl.genes)
	# css.clust - table with all clustered css (see clust.threshold above)
	css.clust <- reduce(css, min.gapwidth = clust.threshold)
	css.clust$clust_id <- c(1:length(css.clust))
	# css.clust <- ens_overlap(css.clust, ensembl.genes)
	# link css and css.clust by clust_id
	overlaps <- as.data.frame(findOverlaps(css, css.clust))
	css$clust_id <- css.clust[overlaps[ , 2]]$clust_id
	css.clust <- css.clust[overlaps[ , 2]]
	css.clust$gene_id <- css$gene_id

	css.clust <- data.table(unique(as.data.frame(sort(css.clust))))
	css <- data.table(as.data.frame(css))

	# print(dim(css))
	# in each cluster in each gene in each condition select css with max score
	setkeyv(css, c("gene_id", "clust_id", "cond"))
	css <-
	unique(css[, list(seqnames, start, end, width, strand, max.score = max(score), 
		sum.score = sum(score)), by = key(css)])

	# from each gene, each condition remove css which are used in less than 10% cases
	setkeyv(css, c("gene_id", "cond"))
	test <- css[, list(clust_id, seqnames, start, end, width, strand, max.score,
		sum.score, cond.usage = sum.score/sum(sum.score)), by = key(css)]
	test <- test[cond.usage > 0.1]

	# index genes and css
	test <- test[, list(clust_id, cond, seqnames, start, end, width, strand,
		max.score, sum.score, gene.ind = length(unique(clust_id))), by = "gene_id"]

	test <- test[order(clust_id)]
	test <- test[, list(clust_id, cond, seqnames, start, end, width, strand,
		max.score, sum.score, gene.ind, css.ind = c(0, diff(clust_id))), by = "gene_id"]
	test[css.ind > 1, css.ind := 1]
	test <- test[, list(clust_id, cond, seqnames, start, end, width, strand,
		max.score, sum.score, gene.ind, css.ind = cumsum(css.ind) + 1), by = "gene_id"]
	test[strand == "+", css.ind := -(css.ind - max(css.ind) - 1), by = "gene_id"]

	# calculate the usage
	setkeyv(test, c("gene_id", "cond"))
	test <- test[, list(clust_id, seqnames, start, end, width, strand, max.score,
		sum.score, gene.ind, css.ind,
		usage = sum.score/sum(sum.score)), by = key(test)]

	# remove genes with distance between css < 300 nt from css.clust
	x <- css.clust[,list(dist = start[-1] - end[-length(end)]),
		by = "gene_id"]
	genes.to.exclude <- x[dist < 300, gene_id]
	test <- test[!(gene_id %in% genes.to.exclude)]
	return(test)
}

css_annotation <- function(d, ind) {
	test <- d[gene.ind <= ind]
	test[,usage_let := "us"]
	test[gene.ind == 2, usage_let := "in_in"]
	test[gene.ind == 2 & css.ind == 1 & usage >= 0.75, usage_let := "us_un"]
	test[gene.ind == 2 & css.ind == 1 & usage <= 0.25, usage_let := "un_us"]
	test[gene.ind == 2 & css.ind == 2 & usage >= 0.75, usage_let := "un_us"]
	test[gene.ind == 2 & css.ind == 2 & usage <= 0.25, usage_let := "us_un"]
	return(test)
}

create_GRange <- function(test) {
	return(
	GRanges(
		seqnames = test[ , seqnames],
		ranges = IRanges(
				 	test[ , start],
				 	test[ , end]
				 ),
		strand = test[ , strand],
		cond = test[ , cond],
		gene.ind = test[ , gene.ind],
		css.ind = test[ , css.ind],
		score = test[ , sum.score],
		usage.pos = test[ , usage],
		usage.let = test[ , usage_let],
		gene_id = test[ , gene_id]
		# clust_id = test[ , clust_id]
		)
	)
}

plot_usage_density <- function(d) {
	p <- ggplot(as.data.frame(d[gene.ind <= 5]), aes(usage)) +
		geom_density(aes(color = cond)) +
		facet_wrap(~ gene.ind + css.ind, scale = "free_y")
	if(data == "a-seq") {
		pdf(file = paste0(report.folder, "a-seq-usage-density.pdf"), 
			width = 10, height = 8)
	}
	if(data == "atlas") {
		pdf(file = paste0(report.folder, "atlas-usage-density.pdf"), 
			width = 10, height = 8)
	}
	print(p)
	dev.off()
}

plot_usage_heatmaps <- function(test) {
	setkeyv(test, c("cond", "usage_let"))
	plot1 <- test[gene.ind == 1,list(freq = length(usage)), by = key(test)]
	plot1[,scaled.freq := freq/sum(freq)]

	plot2 <- test[gene.ind == 2,list(freq = length(usage)), by = key(test)]
	plot2 <- plot2[,list(usage_let, freq, scaled.freq = freq/sum(freq)), by = "cond"]

	p <- ggplot(as.data.frame(plot1), aes(cond, usage_let, fill = scaled.freq)) +
	  geom_tile() +
	  scale_fill_gradient(low = "blue",  high = "yellow", na.value = "blue", name = "Norm. by column") +
	  theme_bw() +
	  geom_text(aes(label = as.character(freq)), hjust=0, show_guide = FALSE, size = 4) +
	  labs(x = "", y = "") +
	  theme(panel.background = element_rect(fill = "blue"), panel.grid.major = element_blank(),
	       panel.grid.minor = element_blank())
	if(data == "a-seq") {
		pdf(file = paste0(report.folder, "a-seq-usage-heatmap1.pdf"), 
			width = 6, height = 3)
	}
	if(data == "atlas") {
		pdf(file = paste0(report.folder, "atlas-usage-heatmap1.pdf"), 
			width = 10, height = 3)
	}
	print(p)
	dev.off()

	p <- ggplot(as.data.frame(plot2), aes(cond, usage_let, fill = scaled.freq)) +
	  geom_tile() +
	  scale_fill_gradient(low = "blue",  high = "yellow", na.value = "blue", name = "Norm. by column") +
	  theme_bw() +
	  geom_text(aes(label = as.character(freq)), hjust=0, show_guide = FALSE, size = 4) +
	  labs(x = "", y = "") +
	  theme(panel.background = element_rect(fill = "blue"), panel.grid.major = element_blank(),
	       panel.grid.minor = element_blank())
	if(data == "a-seq") {
		pdf(file = paste0(report.folder, "a-seq-usage-heatmap2.pdf"), 
			width = 6, height = 6)
	}
	if(data == "atlas") {
		pdf(file = paste0(report.folder, "atlas-usage-heatmap2.pdf"), 
			width = 10, height = 8)
	}
	print(p)
	dev.off()
}
