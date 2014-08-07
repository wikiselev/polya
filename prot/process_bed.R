
# arguments from the command line
args <- commandArgs(TRUE)

# define input, output
data.folder1 <- paste0(getwd(), "/", "../../data/prot/raw/zav/reads/")
data.folder2 <- paste0(getwd(), "/", "../../data/prot/raw/zav/dens/")
data.folder3 <- paste0(getwd(), "/", "../../data/prot/raw/kathi/")
data.folder4 <- paste0(getwd(), "/", "../../data/prot/raw/hnrnps/")
data.folder5 <- paste0(getwd(), "/", "../../data/prot/raw/other/")

out.folder <- paste0(getwd(), "/", "../../data/prot/processed/")

# sorce all my libraries and functions
source(paste0(getwd(), "/", "../resources/ann.R"))

# PAR-Clip protein files from Michaela Zavolan
protein_list1 <- c("CFIm59", "CstF64tau", "CPSF160", "CFIm68", "CstF64",
	"CFIm25", "Fip1", "CPSF73", "CPSF100", "CPSF30")

if (args[1] %in% protein_list1) {
	# import raw reads
	data <- as(import.bed(paste0(data.folder1, args[1], ".bed")), "GRanges")
	# since the reads overlap I want to find a coverage of reads
	# coverage function does not take into account strand, thus first split by
	# strand
	data.pos <- data[strand(data) == "+"]
	data.neg <- data[strand(data) == "-"]
	cov.pos <- coverage(data.pos)
	cov.neg <- coverage(data.neg)
	cov.pos <- as(cov.pos, "GRanges")
	cov.neg <- as(cov.neg, "GRanges")
	seqlengths(cov.pos) <- seqlengths(data)
	seqlengths(cov.neg) <- seqlengths(data)
	strand(cov.pos) = "+"
	strand(cov.neg) = "-"
	cov <- c(cov.pos, cov.neg)
	# remove all intermediate regions which are not covered
	cov <- cov[values(cov)$score >= 1]
	# compute the total.score of the reads (taking into account a score of each
	# read)
	overlaps <- as.data.table(as.data.frame(findOverlaps(cov, data)))
	overlaps[,score:=values(data[overlaps[,subjectHits]])$score]
	overlaps <- overlaps[,list(total.score = sum(score)),by="queryHits"]
	values(cov)$total.score <- overlaps[,total.score]
	# import peak regions information
	peaks <- as(import.bed(paste0(data.folder2, args[1], ".bed")), "GRanges")
	# filter out reads that do not overlap with peak regions
	overlaps <- as.data.frame(findOverlaps(cov, peaks))
	cov.filt <- cov[overlaps[ , 1]]
	# instead of using the whole length of regions I define them as a coordinate of
	# their middle
	tmp <- cov.filt
	start(cov.filt) <- start(tmp) + round((end(tmp) -
		start(tmp))/2)
	end(cov.filt) <- start(tmp) + round((end(tmp) -
		start(tmp))/2)
	saveRDS(cov.filt, paste0(out.folder, args[1], ".RDS"))
}

# iClip protein files from Kathi
protein_list2 <- c("PTB", "TIA1", "TIAL1", "TDP43", "U2AF65", "hnRNPC",
	"TIAL1_Hela", "TDP43_Hela", "TDP43_SHSY5Y", "TDP43_ES")

if (args[1] %in% protein_list2) {
	data <- as(import.bed(paste0(data.folder3, args[1], ".bed")), "GRanges")
	strand(data) <- ifelse(as.numeric(values(data)$name) > 0, "+", "-")
	colnames(values(data))[1] <- "score"
	values(data)$score <- abs(as.numeric(values(data)$score))
	saveRDS(data, paste0(out.folder, args[1], ".RDS"))
}

# # other hnrnps from Yeo_2012 paper
# protein_list3 <- c("hnRNPA1", "hnRNPA2B1", "hnRNPF", "hnRNPH", "hnRNPM",
# 	"hnRNPU")

# if (args[1] %in% protein_list3) {
# 	data <- as(import.bed(paste0(data.folder4, args[1], ".bed")), "GRanges")
# 	starts <- start(values(data)$thick) + round(width(values(data)$thick)/2)
# 	ends <- start(values(data)$thick) + round(width(values(data)$thick)/2)
# 	start(data) <- starts
# 	end(data) <- ends
# 	data <- data[ , 2]
# 	saveRDS(data, paste0(out.folder, args[1]))
# }

if (args[1] == "HuR_Mukherjee") {
	data <- as(import.bed(paste0(data.folder5, args[1], ".bed")), "GRanges")
	data <- data[ , 2]
	starts <- start(data) + round(width(data)/2)
	ends <- start(data) + round(width(data)/2)
	start(data) <- starts
	end(data) <- ends
	saveRDS(data, paste0(out.folder, args[1], ".RDS"))
}
