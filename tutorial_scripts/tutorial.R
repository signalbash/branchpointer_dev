############################################################
###           branchpointer example workflows            ###
############################################################

###### download files ######

# download GTF and .fa files
system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz")
system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa.gz")
system("gunzip gencode.v24.annotation.gtf.gz")
system("gunzip GRCh38.p5.genome.fa.gz")

###### setup ######

library(biomaRt)
library(branchpointer)

options(stringsAsFactors = FALSE)
options(scipen = 999)

# location of the bedtools binary file
# to find location type "which bedtools" in the command line
bedtools <-  "/Applications/apps/bedtools2/bin/bedtools"

exons <- readExonAnnotation("gencode.v24.annotation.gtf")

###### branchpoints in user-defined windows ######

# read file
query_intron <- readQueryFile(system.file("extdata", "intron_example.txt",
                                          package="branchpointer")
                              ,query_type = "region")
# get location
query_intron <- getQueryLoc(query_intron,query_type="region",exons = exons)
# get attributes
query_attributes_intron <- getBranchpointSequence(query_intron,
                                                  query_type = "region",
                                                  genome = "GRCh38.p5.genome.fa",
                                                  bedtools_location = bedtools)
# predict branchpoints
branchpoint_predictions_intron <- predictBranchpoints(query_attributes_intron)
# plot window (without structure)
plotBranchpointWindow(query_intron$id[1], branchpoint_predictions_intron,
                      query_attributes_intron, plot_structure = FALSE)

###### branchpoints in specified genes/transcripts/exons ######

#exon
query_intron <- makeRegions("ENSE00003541068", "exon_id", exons)
#transcript
query_intron <- makeRegions("ENST00000357654", "transcript_id", exons)
#gene
query_intron <- makeRegions("ENSG00000139618", "gene_id", exons)

# get location
query_intron <- getQueryLoc(query_intron,query_type="region",exons = exons)
# get attributes
query_attributes_intron <- getBranchpointSequence(query_intron,
                                                  query_type = "region",
                                                  genome = "GRCh38.p5.genome.fa",
                                                  bedtools_location = bedtools)
# predict branchpoints
branchpoint_predictions_intron <- predictBranchpoints(query_attributes_intron)
# plot branchpoints in intron 1
plotBranchpointWindow(query_intron$id[1], branchpoint_predictions_intron,
                      query_attributes_intron, exons = exons)

###### branchpoints in user-defined SNPs ######
# read file
query_snp <- readQueryFile(system.file("extdata", "SNP_example.txt",
                                       package = "branchpointer"),
                           query_type = "SNP")
# get location
query_snp <- getQueryLoc(query_snp,query_type="SNP",exons = exons, filter = FALSE)
# get attributes
query_attributes_snp <- getBranchpointSequence(query_snp,
                                               query_type = "SNP",
                                               genome = "GRCh38.p5.genome.fa",
                                               bedtools_location = bedtools)
# predict branchpoints
branchpoint_predictions_snp <- predictBranchpoints(query_attributes_snp)
# evaluate SNP effects
snp_stats <- predictionsToStats(branchpoint_predictions_snp, query_snp)
# plot branchpoints in reference and alternative sequences (without structure plot)
plotBranchpointWindow(snp_stats$id[2], branchpoint_predictions_snp,
                      query_attributes_snp, plot_structure = FALSE,
                      plot_mutated = TRUE)

###### branchpoints in SNPs with RefSNP IDs ######

# create mart object
mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp",host="www.ensembl.org")
# make query data.frame
query_snp <- snpToQuery(c("rs17000647","rs5031002","rs998731"), mart_snp = mart)

# get location
query_snp <- getQueryLoc(query_snp,query_type="SNP",exons = exons, filter = FALSE)
# get attributes
query_attributes_snp <- getBranchpointSequence(query_snp,
                                               query_type = "SNP",
                                               genome = "GRCh38.p5.genome.fa",
                                               bedtools_location = bedtools)
# predict branchpoints
branchpoint_predictions_snp <- predictBranchpoints(query_attributes_snp)
# evaluate SNP effects
snp_stats <- predictionsToStats(branchpoint_predictions_snp, query_snp)
# plot branchpoints in reference and alternative sequences
plotBranchpointWindow(snp_stats$id[1], branchpoint_predictions_snp,
                      query_attributes_snp, exons = exons,
                      plot_mutated = TRUE)
