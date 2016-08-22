############################################################
###      Effects of NHGRI GWAS SNPs on branchpoints      ###
############################################################

###### download files ######

# download GTF and .fa files
system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz")
system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa.gz")
system("gunzip gencode.v24.annotation.gtf.gz")
system("gunzip GRCh38.p5.genome.fa.gz")

# download NHGRI GWAS file
# date: 22/08/16
system("wget https://www.ebi.ac.uk/gwas/api/search/downloads/full -O gwas_catalog_v1.0.1-associations_e85-r2016-08-14.tsv")

###### run branchpointer ######

library(biomaRt)
library(branchpointer)

options(stringsAsFactors = FALSE)
options(scipen = 999)

# location of the bedtools binary file
# to find location type "which bedtools" in the command line
bedtools <-  "/Applications/apps/bedtools2/bin/bedtools"

# format exon annotation
exons <- readExonAnnotation("gencode.v24.annotation.gtf")

# read in GWAS annotation
NHGRI_GWAS <- read.delim("gwas_catalog_v1.0.1-associations_e85-r2016-08-14.tsv",
                         stringsAsFactors = FALSE)

# convert 'intron_variant' rs ids to query format
mart_snp <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
query <- snpToQuery(NHGRI_GWAS$SNPS[NHGRI_GWAS$CONTEXT == "intron_variant"],
                    mart = mart_snp)

# get location of SNPs
query <- getQueryLoc(query,query_type = "SNP",
                     exons = exons, filter = TRUE)

# get all features for branchpoint prediction
query_attributes <- getBranchpointSequence(query,
                                        query_type = "SNP",
                                        genome = "GRCh38.p5.genome.fa",
                                        bedtools_location = bedtools)
# predict branchpoints
branchpoint_predictions <- predictBranchpoints(query_attributes)

# filter for SNPs that create or delete branchpoints
snp_stats <- predictionsToStats(branchpoint_predictions, query)
snp_stats <- snp_stats[snp_stats$created_n > 0 | snp_stats$deleted_n > 0,]

snp_stats$U2_diff <- abs(snp_stats$max_U2_REF - snp_stats$max_U2_ALT)
snp_stats <- arrange(snp_stats, plyr::desc(U2_diff))

# plot rs17000647
pdf("rs17000647_branchpoints.pdf", width = 8, height = 6)
plotBranchpointWindow(snp_stats$id[1], branchpoint_predictions,
                      query_attributes, plot_mutated = T,
                      plot_structure = T, exons = exons)
dev.off()

