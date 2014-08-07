# arguments from the command line
args <- commandArgs(TRUE)
# define job id from job array
iscen <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

# define input, output and annotation folders
source(paste0(getwd(), "/", "../resources/ann.R"))
prot.folder <- paste0(getwd(), "/", "../../data/prot/processed/")
out.folder <- paste0(getwd(), "/", "../../data/map/", args[1], "/")
css.folder <- paste0(getwd(), "/", "../../data/css/processed/")

# pick up protein files according to iscen
prot.files <- list.files(path = prot.folder)
prot.names <- sapply(strsplit(prot.files, "\\."), "[[", 1)
prot.file <- prot.files[iscen]
prot.name <- prot.names[iscen]

prot.data <- readRDS(paste0(prot.folder, prot.file))
css.data <- readRDS(paste0(css.folder, args[1], ".RDS"))

# define transcript flanking region around transcript cleavage sites 
flank.length <- 150
transcripts <- flank(css.data, -flank.length, both = TRUE, start = FALSE)

# load protein data file
seqlevels(prot.data, force = TRUE) <- seqlevels(transcripts)

overlaps <- as.data.frame(findOverlaps(prot.data, transcripts))
trs <- as.data.frame(transcripts[overlaps[ , 2]])
mapped <- as.data.frame(prot.data[overlaps[ , 1]])
mapped$start <- (mapped$start - (trs$start + flank.length))*ifelse(mapped$strand == "+", 1, -1)
mapped$end <- (mapped$end - (trs$start + flank.length))*ifelse(mapped$strand == "+", 1, -1)
mapped <- mapped[ , c(2, 6)]
mapped <- cbind(mapped, trs$cond)
mapped <- cbind(mapped, trs$gene.ind)
mapped <- cbind(mapped, trs$css.ind)
mapped <- cbind(mapped, trs$score)
mapped <- cbind(mapped, trs$usage.pos)
mapped <- cbind(mapped, trs$usage.let)
mapped <- cbind(mapped, trs$gene_id)
mapped$prot <- rep(prot.name, length(mapped[ ,1]))
colnames(mapped) <- c("coord", "score", "cond", "gene.id", "css.id", "usage.score",
	"usage.pos", "usage.let", "gene_id", "prot")

saveRDS(mapped, paste0(out.folder, prot.file))
