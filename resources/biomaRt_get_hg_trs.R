# This script retrieves annotations from ensembl biomart database. Here I am 
# mostly interested in transcripts annotation. However, one could also retrieve
# what is needed by changing the list of attributes "attr". The full list of
# attributes is available by calling listAttributes(mart)
#
# Author: Vladimir Kiselev

source(paste0(getwd(), "/", "../resources/ann.R"))

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# retrieving annotations from ensembl
annotations <- data.frame()

# biomaRt attributes
attr <- c("ensembl_gene_id", 
		  "ensembl_transcript_id",
		  "chromosome_name",
		  "start_position",
		  "end_position",
		  "strand",
		  "transcript_start",
		  "transcript_end",
		  "gene_biotype",
		  "transcript_biotype")

annotations <- getBM(attributes = attr, mart = mart)

# change strand to +/-
annotations$strand <- ifelse(annotations$strand > 0, "+", "-")

# keep only annotations on conventional chromosomes
annotations <- subset(annotations, chromosome_name %in% c(1:22, "X", "Y"))

date <- format(Sys.time(), "%d%b%Y")

saveRDS(annotations, paste0("hsapiens_trs_biomaRt_table_", date, ".RDS"))

# create a GRanges object for protein coding genes only
gene_ann <- annotations[ , c(1, 3:6, 9)]
props <- c("protein_coding")
genes <- GRanges(
	seqnames = gene_ann$chromosome_name[gene_ann$gene_biotype %in% props],
	ranges = IRanges(
			 	gene_ann$start_position[gene_ann$gene_biotype %in% props],
			 	gene_ann$end_position[gene_ann$gene_biotype %in% props]
			 ),
	strand = gene_ann$strand[gene_ann$gene_biotype %in% props],
	gene_id = gene_ann$ensembl_gene_id[gene_ann$gene_biotype %in% props]
)
# put some make up on the object
genes <- unique(genes)
seqlevels(genes) <- paste("chr", seqlevels(genes), sep="")
seqlengths(genes) <- seqlengths(Hsapiens)[names(seqlengths(genes))]

saveRDS(genes, paste0("hsapiens_prot_cod_trs_biomaRt_granges_", date, ".RDS"))
