#!/bin/bash

root="../../data/prot/raw/zav"

# obtain raw reads data from the files downloaded from CLIPZ server (they were
# shared with me by Andreas Gruber):
# http://www.clipz.unibas.ch/publicFiles/
# this url is still working on 27/11/2013
echo "perl clipz2bed.pl $root/reads/sample_250/mapped_sequences $root/reads/sample_250/genome_mappings > $root/reads/CFIm25.bed" | bsub -J CFIm25_read
echo "perl clipz2bed.pl $root/reads/sample_244/mapped_sequences $root/reads/sample_244/genome_mappings > $root/reads/CstF64tau.bed" | bsub -J CstF64tau_read
echo "perl clipz2bed.pl $root/reads/sample_240/mapped_sequences $root/reads/sample_240/genome_mappings > $root/reads/CstF64.bed" | bsub -J CstF64_read
echo "perl clipz2bed.pl $root/reads/sample_248/mapped_sequences $root/reads/sample_248/genome_mappings > $root/reads/Fip1.bed" | bsub -J Fip1_read
echo "perl clipz2bed.pl $root/reads/sample_247/mapped_sequences $root/reads/sample_247/genome_mappings > $root/reads/CPSF160.bed" | bsub -J CPSF160_read
echo "perl clipz2bed.pl $root/reads/sample_899/mapped_sequences $root/reads/sample_899/genome_mappings > $root/reads/CPSF100.bed" | bsub -J CPSF100_read
echo "perl clipz2bed.pl $root/reads/sample_246/mapped_sequences $root/reads/sample_246/genome_mappings > $root/reads/CPSF73.bed" | bsub -J CPSF73_read
echo "perl clipz2bed.pl $root/reads/sample_245/mapped_sequences $root/reads/sample_245/genome_mappings > $root/reads/CPSF30.bed" | bsub -J CPSF30_read
echo "perl clipz2bed.pl $root/reads/sample_239/mapped_sequences $root/reads/sample_239/genome_mappings > $root/reads/CFIm68.bed" | bsub -J CFIm68_read
echo "perl clipz2bed.pl $root/reads/sample_252/mapped_sequences $root/reads/sample_252/genome_mappings > $root/reads/CFIm59.bed" | bsub -J CFIm59_read

# process the WIG files downloaded from NCBI GEO database:
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37401
# wig2bed.py is simply creates BED files from WIG files while losing the density
# information.
# instead of the density I am using raw reads obtained using clipz2bed.pl before
bsub -J CPSF160_pos_wig2bed python wig2bed.py $root/dens/CPSF160_+_density.wig $root/dens/CPSF160_+_region.bed
bsub -J CPSF160_neg_wig2bed python wig2bed.py $root/dens/CPSF160_-_density.wig $root/dens/CPSF160_-_region.bed
bsub -J Fip1_pos_wig2bed python wig2bed.py $root/dens/Fip1_+_density.wig $root/dens/Fip1_+_region.bed
bsub -J Fip1_neg_wig2bed python wig2bed.py $root/dens/Fip1_-_density.wig $root/dens/Fip1_-_region.bed
bsub -J CstF64_neg_wig2bed python wig2bed.py $root/dens/CstF-64_-_density.wig $root/dens/CstF64_-_region.bed
bsub -J CstF64_pos_wig2bed python wig2bed.py $root/dens/CstF-64_+_density.wig $root/dens/CstF64_+_region.bed
bsub -J CstF64tau_pos_wig2bed python wig2bed.py $root/dens/CstF-64tau_+_density.wig $root/dens/CstF64tau_+_region.bed
bsub -J CstF64tau_neg_wig2bed python wig2bed.py $root/dens/CstF-64tau_-_density.wig $root/dens/CstF64tau_-_region.bed
bsub -J CPSF100_pos_wig2bed python wig2bed.py $root/dens/CPSF-100_+_density.wig $root/dens/CPSF100_+_region.bed
bsub -J CPSF100_neg_wig2bed python wig2bed.py $root/dens/CPSF-100_-_density.wig $root/dens/CPSF100_-_region.bed
bsub -J CPSF73_pos_wig2bed python wig2bed.py $root/dens/CPSF-73_+_density.wig $root/dens/CPSF73_+_region.bed
bsub -J CPSF73_neg_wig2bed python wig2bed.py $root/dens/CPSF-73_-_density.wig $root/dens/CPSF73_-_region.bed
bsub -J CPSF30_pos_wig2bed python wig2bed.py $root/dens/CPSF-30_+_density.wig $root/dens/CPSF30_+_region.bed
bsub -J CPSF30_neg_wig2bed python wig2bed.py $root/dens/CPSF-30_-_density.wig $root/dens/CPSF30_-_region.bed
bsub -J CFIm68_pos_wig2bed python wig2bed.py $root/dens/CFIm68_+_density.wig $root/dens/CFIm68_+_region.bed
bsub -J CFIm68_neg_wig2bed python wig2bed.py $root/dens/CFIm68_-_density.wig $root/dens/CFIm68_-_region.bed
bsub -J CFIm59_pos_wig2bed python wig2bed.py $root/dens/CFIm59_+_density.wig $root/dens/CFIm59_+_region.bed
bsub -J CFIm59_neg_wig2bed python wig2bed.py $root/dens/CFIm59_-_density.wig $root/dens/CFIm59_-_region.bed
bsub -J CFIm25_pos_wig2bed python wig2bed.py $root/dens/CFIm25_+_density.wig $root/dens/CFIm25_+_region.bed
bsub -J CFIm25_neg_wig2bed python wig2bed.py $root/dens/CFIm25_-_density.wig $root/dens/CFIm25_-_region.bed

# concatenate region + and - files obtained in the previous step
echo "cat $root/dens/CPSF160_+_region.bed $root/dens/CPSF160_-_region.bed > $root/dens/CPSF160.bed" | bsub -w "ended("CPSF160_pos_wig2bed") && ended("CPSF160_neg_wig2bed")" -J CPSF160_reg_cat
echo "cat $root/dens/Fip1_+_region.bed $root/dens/Fip1_-_region.bed > $root/dens/Fip1.bed" | bsub -w "ended("Fip1_pos_wig2bed") && ended("Fip1_neg_wig2bed")" -J Fip1_reg_cat
echo "cat $root/dens/CstF64_+_region.bed $root/dens/CstF64_-_region.bed > $root/dens/CstF64.bed" | bsub -w "ended("CstF64_pos_wig2bed") && ended("CstF64_neg_wig2bed")" -J CstF64_reg_cat
echo "cat $root/dens/CstF64tau_+_region.bed $root/dens/CstF64tau_-_region.bed > $root/dens/CstF64tau.bed" | bsub -w "ended("CstF64tau_pos_wig2bed") && ended("CstF64tau_neg_wig2bed")" -J CstF64tau_reg_cat
echo "cat $root/dens/CPSF100_+_region.bed $root/dens/CPSF100_-_region.bed > $root/dens/CPSF100.bed" | bsub -w "ended("CPSF100_pos_wig2bed") && ended("CPSF100_neg_wig2bed")" -J CPSF100_reg_cat
echo "cat $root/dens/CPSF73_+_region.bed $root/dens/CPSF73_-_region.bed > $root/dens/CPSF73.bed" | bsub -w "ended("CPSF73_pos_wig2bed") && ended("CPSF73_neg_wig2bed")" -J CPSF73_reg_cat
echo "cat $root/dens/CPSF30_+_region.bed $root/dens/CPSF30_-_region.bed > $root/dens/CPSF30.bed" | bsub -w "ended("CPSF30_pos_wig2bed") && ended("CPSF30_neg_wig2bed")" -J CPSF30_reg_cat
echo "cat $root/dens/CFIm68_+_region.bed $root/dens/CFIm68_-_region.bed > $root/dens/CFIm68.bed" | bsub -w "ended("CFIm68_pos_wig2bed") && ended("CFIm68_neg_wig2bed")" -J CFIm68_reg_cat
echo "cat $root/dens/CFIm59_+_region.bed $root/dens/CFIm59_-_region.bed > $root/dens/CFIm59.bed" | bsub -w "ended("CFIm59_pos_wig2bed") && ended("CFIm59_neg_wig2bed")" -J CFIm59_reg_cat
echo "cat $root/dens/CFIm25_+_region.bed $root/dens/CFIm25_-_region.bed > $root/dens/CFIm25.bed" | bsub -w "ended("CFIm25_pos_wig2bed") && ended("CFIm25_neg_wig2bed")" -J CFIm25_reg_cat

# delete region + and - files leaving only concatenated files
bsub -w "ended("CPSF160_reg_cat")" rm $root/dens/CPSF160_+_region.bed $root/dens/CPSF160_-_region.bed
bsub -w "ended("Fip1_reg_cat")" rm $root/dens/Fip1_+_region.bed $root/dens/Fip1_-_region.bed
bsub -w "ended("CstF64_reg_cat")" rm $root/dens/CstF64_+_region.bed $root/dens/CstF64_-_region.bed
bsub -w "ended("CstF64tau_reg_cat")" rm $root/dens/CstF64tau_+_region.bed $root/dens/CstF64tau_-_region.bed
bsub -w "ended("CPSF100_reg_cat")" rm $root/dens/CPSF100_+_region.bed $root/dens/CPSF100_-_region.bed
bsub -w "ended("CPSF73_reg_cat")" rm $root/dens/CPSF73_+_region.bed $root/dens/CPSF73_-_region.bed
bsub -w "ended("CPSF30_reg_cat")" rm $root/dens/CPSF30_+_region.bed $root/dens/CPSF30_-_region.bed
bsub -w "ended("CFIm68_reg_cat")" rm $root/dens/CFIm68_+_region.bed $root/dens/CFIm68_-_region.bed
bsub -w "ended("CFIm59_reg_cat")" rm $root/dens/CFIm59_+_region.bed $root/dens/CFIm59_-_region.bed
bsub -w "ended("CFIm25_reg_cat")" rm $root/dens/CFIm25_+_region.bed $root/dens/CFIm25_-_region.bed
