# TODO: none
#
# Author: Vladimir Kiselev
###############################################################################

# define input, output and annotation folders
annotation.folder <- paste(getwd(), "/", "../../../../../projects/annotations/", 
							sep="") 
data.folder <- paste(getwd(), "/", "../../data/mine/Martin_2012/CLIP/raw_from_clipz_server/reads/", sep="")
data.folder1 <- paste(getwd(), "/", "../../data/mine/Martin_2012/CLIP/", sep="")
out.folder <- paste(getwd(), "/", "../../data/new_with_scores/", sep = "")

# sorce all my libraries and functions
source(paste(annotation.folder, "annotations.r", sep = ""))

files <- list.files(path = data.folder)
files <- files[grepl(".RData", files)]
files <- files[!grepl("_rep_B.RData", files)]

# CFIm59 and CstF64tau are too big compared to others -- I analyze them using
# separate scripts (to be able to allocate more memory to a process)
files <- subset(files, files != "CFIm59.RData" & files != "CstF64tau.RData")

for(file in files) {
	load(paste(data.folder, file, sep = ""))

	# obtain coverage of uniquely mapped reads
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
	cov <- cov[values(cov)$score >= 1]

	# import peak coordinates from read density .wig files (I assume the are
	# the largest ones in the whole coverage)
	peaks <- import.bed(paste(data.folder1, strsplit(file, "\\.")[[1]][1], "_region.bed", sep = ""))
	peaks <- as(peaks, "GRanges")

	# find overlaps between raw coverage and processed read density
	cov.combined <- reduce(cov)
	overlaps <- as.data.frame(findOverlaps(cov.combined, peaks))
	cov.filt <- cov.combined[overlaps[ , 1]]

	# find region of the coverage with the highest read number (in one peak)
	overlaps1 <- as.data.frame(findOverlaps(cov, cov.filt))
	overlaps1$score <- values(cov[overlaps1[ , 1]])$score
	DT <- data.table(overlaps1)
	DT.MAX <- DT[,max(score),by=subjectHits]
	setnames(DT.MAX, "V1", "score")
	DT.MERGE1 <- merge(DT, DT.MAX, by=c("subjectHits"))
	DT.MERGE1 <- DT.MERGE1[score.x==score.y]
	# here are regions of the coverage with the highest read number (in one peak)
	cov.final <- cov[DT.MERGE1[,queryHits]]

	# define the center of the peak region -- this is the final peak coordinate!
	cov.final.center <- cov.final
	start(cov.final.center) <- start(cov.final) + round((end(cov.final) -
		start(cov.final))/2)
	end(cov.final.center) <- start(cov.final) + round((end(cov.final) -
		start(cov.final))/2)

	# I want to put the total read number in the peak are inside the peak
	# coordinate
	DT.SUM <- DT[,sum(score),by=subjectHits]
	DT.MERGE2 <- merge(DT.SUM, DT.MAX, by=c("subjectHits"))
	# setnames(DT.SUM, "V1", "score")
	DT.MERGE3 <- merge(DT, DT.MERGE2, by=c("subjectHits"))
	DT.TOTAL <- DT.MERGE3[score.x==score.y]
	values(cov.final.center)$total_peak_score <- DT.TOTAL[,V1]

	total.read.num <- length(data)

	values(cov.final.center)$score_norm <-
		values(cov.final.center)$score/total.read.num
	values(cov.final.center)$total_peak_score_norm <-
		values(cov.final.center)$total_peak_score/total.read.num

	data <- cov.final.center

	print(file)
	print(range(values(cov.final.center)$score))
	print(range(values(cov.final.center)$total_peak_score))
	print(total.read.num)
	print(range(values(cov.final.center)$score_norm))
	print(range(values(cov.final.center)$total_peak_score_norm))
	print(range(width(cov.final.center)))
	save(data, file = paste(out.folder, file, sep = ""))
}

# [1] "CFIm25.RData"
# [1]    1 1157
# [1]      1 115869
# [1] 1948398
# [1] 5.132422e-07 5.938212e-04
# [1] 5.132422e-07 5.946886e-02
# [1]  1 98
# [1] "CFIm68.RData"
# [1]    1 2026
# [1]      1 205789
# [1] 5611430
# [1] 1.782077e-07 3.610488e-04
# [1] 1.782077e-07 3.667318e-02
# [1]   1 114
# [1] "CPSF100.RData"
# [1]    1 2941
# [1]     1 86728
# [1] 116280
# [1] 8.599931e-06 2.529240e-02
# [1] 8.599931e-06 7.458548e-01
# [1]  1 58
# [1] "CPSF160.RData"
# [1]    1 2767
# [1]      1 936055
# [1] 8050193
# [1] 1.242206e-07 3.437185e-04
# [1] 1.242206e-07 1.162773e-01
# [1]   1 114
# [1] "CPSF30.RData"
# [1]    1 1261
# [1]     1 79143
# [1] 350652
# [1] 2.851830e-06 3.596158e-03
# [1] 2.851830e-06 2.257024e-01
# [1]  1 72
# [1] "CPSF73.RData"
# [1]    1 1278
# [1]     1 63442
# [1] 413431
# [1] 2.418783e-06 3.091205e-03
# [1] 2.418783e-06 1.534525e-01
# [1]  1 65
# [1] "CstF64.RData"
# [1]    1 2957
# [1]      1 396740
# [1] 5396142
# [1] 1.853176e-07 5.479841e-04
# [1] 1.853176e-07 7.352290e-02
# [1]   1 114
# [1] "Fip1.RData"
# [1]   1 845
# [1]     1 51507
# [1] 1441379
# [1] 6.937801e-07 5.862441e-04
# [1] 6.937801e-07 3.573453e-02
# [1]  1 80

kathi_files <- c("PTB.bed", "TIA1.bed", "TIAL1.bed", "TDP43.bed", "U2AF65.bed",
	"hnRNPC.bed")

# this is for clusters
# kathi_files <- c("PTB_clusters.bed", "TIA1_clusters.bed", "TIAL1_clusters.bed",
# 	"TDP43_clusters.bed", "U2AF65_clusters.bed", "hnRNPC_clusters.bed")

for(file in kathi_files) {
	data <- as(import.bed(paste(getwd(), "/", "../../data/from_Kathi/", file, 
				sep = "")), "GRanges")

	strand(data) <- ifelse(as.numeric(values(data)$name) > 0, "+", "-")
	colnames(values(data))[1] <- "score"
	values(data)$score <- abs(as.numeric(values(data)$score))

	# this is for clusters
	# file <- paste(strsplit(file, "_")[[1]][1], ".RData", sep = "")

	file <- paste(strsplit(file, ".")[[1]][1], ".RData", sep = "")

	print(file)
	print(range(values(data)$score))
	print(sum(values(data)$score))
	print(range(width(data)))

	# this is for clusters
	# starts <- start(data) + round(width(data)/2)
	# ends <- start(data) + round(width(data)/2)
	# start(data) <- starts
	# end(data) <- ends

	save(data, file = paste(out.folder, file, sep = ""))
}

# these are new files from Jurnej - 21/11/2013

kathi_files <- c("TIAL1_Hela.bed", "TDP43_Hela.bed", "TDP43_SHSY5Y.bed", "TDP43_ES.bed")

# this is for clusters
# kathi_files <- c("PTB_clusters.bed", "TIA1_clusters.bed", "TIAL1_clusters.bed",
# 	"TDP43_clusters.bed", "U2AF65_clusters.bed", "hnRNPC_clusters.bed")

for(file in kathi_files) {
	data <- as(import.bed(paste(getwd(), "/", "../../data/from_Kathi/", file, 
				sep = "")), "GRanges")

	strand(data) <- ifelse(as.numeric(values(data)$name) > 0, "+", "-")
	colnames(values(data))[1] <- "score"
	values(data)$score <- abs(as.numeric(values(data)$score))

	# this is for clusters
	# file <- paste(strsplit(file, "_")[[1]][1], ".RData", sep = "")

	file <- paste(strsplit(file, "\\.")[[1]][1], ".RData", sep = "")

	print(file)
	print(range(values(data)$score))
	print(sum(values(data)$score))
	print(range(width(data)))

	# this is for clusters
	# starts <- start(data) + round(width(data)/2)
	# ends <- start(data) + round(width(data)/2)
	# start(data) <- starts
	# end(data) <- ends

	save(data, file = paste(out.folder, file, sep = ""))
}

# [1] "PTB.RData"
# [1]     1 13030
# [1] 2453937
# [1] 4.075084e-07 5.309835e-03
# [1]    1 1474
# [1] "TIA1.RData"
# [1]    1 4787
# [1] 291978
# [1] 3.424916e-06 1.639507e-02
# [1]    1 1123
# [1] "TIAL1.RData"
# [1]     1 24885
# [1] 868022
# [1] 1.152045e-06 2.866863e-02
# [1]    1 2000
# [1] "TDP43.RData"
# [1]     1 14202
# [1] 799892
# [1] 1.250169e-06 1.775490e-02
# [1]    1 1997
# [1] "U2AF65.RData"
# [1]     1 70601
# [1] 5693413
# [1] 1.756416e-07 1.240047e-02
# [1]    1 1383
# [1] "hnRNPC.RData"
# [1]     1 68866
# [1] 3943427
# [1] 2.535865e-07 1.746349e-02
# [1]    1 3294



hnrnp_files <- c("from_Kathi/hnRNPA1_Yeo_2012/hnRNPA1_Yeo_2012_clusts_hg19.bed",
	"from_Kathi/hnRNPA2B1_Yeo_2012/hnRNPA2B1_Yeo_2012_clusts_hg19.bed",
	"from_Kathi/hnRNPF_Yeo_2012/hnRNPF_Yeo_2012_clusts_hg19.bed",
	"from_Kathi/hnRNPH_Burge/hnRNPH_Burge_clusts_hg19.bed",
	"from_Kathi/hnRNPM_Yeo_2012/hnRNPM_Yeo_2012_clusts_hg19.bed",
	"from_Kathi/hnRNPU_Yeo_2012/hnRNPU_Yeo_2012_clusts_hg19.bed")

for(file in hnrnp_files) {
	data <- as(import.bed(paste(getwd(), "/", "../../data/", file, 
				sep = "")), "GRanges")
	file <- strsplit(file, "/")[[1]][2]
	file <- strsplit(file, "_")[[1]][1]
	file <- paste(file, ".RData", sep = "")
	print(file)
	print(range(values(data)$score))
	print(sum(values(data)$score))
	# print(range(values(data)$score)/sum(values(data)$score))
	print(range(width(data)))

	# values(data)$score <- values(data)$score/sum(values(data)$score)

	starts <- start(values(data)$thick) + round(width(values(data)$thick)/2)
	ends <- start(values(data)$thick) + round(width(values(data)$thick)/2)
	start(data) <- starts
	end(data) <- ends

	data <- data[ , 2]

	save(data, file = paste(out.folder, file, sep = ""))
}

# [1] "hnRNPA1_Yeo_2012_clusts_hg19.RData"
# [1]    5 1034
# [1] 24863
# [1] 0.000201102 0.041587902
# [1]   2 999
# [1] "hnRNPA2B1_Yeo_2012_clusts_hg19.RData"
# [1]    5 1168
# [1] 83772
# [1] 5.968581e-05 1.394261e-02
# [1]   1 998
# [1] "hnRNPF_Yeo_2012_clusts_hg19.RData"
# [1]    5 1737
# [1] 152372
# [1] 3.281443e-05 1.139973e-02
# [1]    1 1000
# [1] "hnRNPH_Burge_clusts_hg19.RData"
# [1]   7 843
# [1] 602224
# [1] 1.162358e-05 1.399811e-03
# [1]    1 1000
# [1] "hnRNPM_Yeo_2012_clusts_hg19.RData"
# [1]    6 1552
# [1] 75914
# [1] 0.0000790368 0.0204441868
# [1]    3 1000
# [1] "hnRNPU_Yeo_2012_clusts_hg19.RData"
# [1]    6 1245
# [1] 155705
# [1] 3.853441e-05 7.995890e-03
# [1]    1 1000

# PTB Yeo from Kathi
file <- "from_Kathi/PTB_Yeo_2009/PTB_Yeo_2009_clusts_hg19.RData"
load(paste(getwd(), "/", "../../data/", file, sep = ""))
data <- clusts
data <- data[ , 2]

starts <- start(data) + round(width(data)/2)
ends <- start(data) + round(width(data)/2)
start(data) <- starts
end(data) <- ends

colnames(values(data)) <- "score"

file <- "PTB_Yeo_Kathi.RData"
save(data, file = paste(out.folder, file, sep = ""))


file <- "from_Kathi/PTB_Xue_2009.bed"
data <- as(import.bed(paste(getwd(), "/", "../../data/", file, sep = "")),
	"GRanges")
data <- data[ , 2]

file <- "PTB_Yeo_Publ.RData"

print(file)
print(range(values(data)$score))
print(sum(values(data)$score))
print(range(width(data)))

starts <- start(data) + round(width(data)/2)
ends <- start(data) + round(width(data)/2)
start(data) <- starts
end(data) <- ends

save(data, file = paste(out.folder, file, sep = ""))

# [1] "PTB.RData"
# [1]   2 936
# [1] 453160
# [1] 4.413452e-06 2.065496e-03
# [1]   1 792

file <- "mine/HuR_Mukherjee_2011.bed"

data <- as(import.bed(paste(getwd(), "/", "../../data/", file, sep = "")),
	"GRanges")
data <- data[ , 2]

file <- strsplit(file, "/")[[1]][2]
file <- paste(strsplit(file, "_")[[1]][1], strsplit(file, "_")[[1]][2], sep = "_")
file <- paste(file, ".RData", sep = "")

print(file)
print(range(values(data)$score))
print(sum(values(data)$score))
# print(range(values(data)$score)/sum(values(data)$score))
print(range(width(data)))

# values(data)$score <- values(data)$score/sum(values(data)$score)

starts <- start(data) + round(width(data)/2)
ends <- start(data) + round(width(data)/2)
start(data) <- starts
end(data) <- ends

save(data, file = paste(out.folder, file, sep = ""))

# [1] "HuR_Mukherjee.RData"
# [1]  1.00000 13.48344
# [1] 534131.8
# [1] 1.872197e-06 2.524365e-05
# [1]   9 111

file <- "mine/HuR_Lebedeva_2011.bed"

data <- as(import.bed(paste(getwd(), "/", "../../data/", file, 
			sep = "")), "GRanges")
file <- strsplit(file, "/")[[1]][2]
file <- paste(strsplit(file, "_")[[1]][1], strsplit(file, "_")[[1]][2], sep = "_")
file <- paste(file, ".RData", sep = "")

print(file)
print(range(values(data)$score))
print(sum(values(data)$score))
# print(range(values(data)$score)/sum(values(data)$score))
print(range(width(data)))

starts <- start(values(data)$thick) + round(width(values(data)$thick)/2)
ends <- start(values(data)$thick) + round(width(values(data)$thick)/2)
start(data) <- starts
end(data) <- ends

# values(data)$score <- values(data)$score/sum(values(data)$score)

data <- data[ , 2]

save(data, file = paste(out.folder, file, sep = ""))

# [1] "HuR_Lebedeva.RData"
# [1]    1 1000
# [1] 6616036
# [1] 1.511479e-07 1.511479e-04
# [1]  15 324

# file <- "from_Kathi/PTB_Yeo_2009/PTB_Yeo_2009_clusts_hg19.bed" -- this one can't
# be imported -- need to check what the problem

