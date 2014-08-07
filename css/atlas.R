source("css_func.R")
data.folder <- "../../data/css/raw/atlas/"

# import data files
brain <- import.bed(paste(data.folder,
	"GSM747470_human_brain.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
testis <- import.bed(paste(data.folder,
	"GSM747479_human_testis.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
muscle <- import.bed(paste(data.folder,
	"GSM747477_human_muscle.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
liver <- import.bed(paste(data.folder,
	"GSM747472_human_liver.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
maqcUHR1 <- import.bed(paste(data.folder,
	"GSM747475_human_maqc-UHR1.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
maqcUHR2 <- import.bed(paste(data.folder,
	"GSM747476_human_maqc-UHR2.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
maqcbrain1 <- import.bed(paste(data.folder,
	"GSM747473_human_maqc-brain1.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)
maqcbrain2 <- import.bed(paste(data.folder,
	"GSM747474_human_maqc-brain2.sites.clustered.hg19.bed", sep=""), genome = "hg19",
	asRangedData = FALSE)

values(brain)$cond <- "brain"
values(testis)$cond <- "testis"
values(muscle)$cond <- "muscle"
values(liver)$cond <- "liver"
values(maqcUHR1)$cond <- "maqcUHR1"
values(maqcUHR2)$cond <- "maqcUHR2"
values(maqcbrain1)$cond <- "maqcbrain1"
values(maqcbrain2)$cond <- "maqcbrain2"

css <- c(brain, testis, muscle, liver, maqcUHR1, maqcUHR2, maqcbrain1, maqcbrain2)

data <- "atlas"
score.threshold <- 0
clust.threshold <- 10
css.num <- 2  # consider only genes with more than 3 css

source("common_processing.R")

saveRDS(css, paste0(out.folder, "atlas.RDS"))
