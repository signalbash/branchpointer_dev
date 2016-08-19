#' Convert exon annotation data.frame to intron locations
#'
#' Converts exon annotation to intron locations overlapping the branchpoint region
#' for exculsion of non-branchpoint region SNPs
#' Returns a character vector of chromosome_location
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param max_dist Maximum distance from the 3' exon to create the branchpoint region.
#' @return vector of chromosome names and intronic locations
#' @export
#' @examples
#' exons <- readExonAnnotation("gencode.v24.exons.txt")
#' introns <- exonsToIntrons(exons, 50)
#' @author Beth Signal

exonsToIntrons <- function(exons, max_dist = 50){

  #split exon annotation by strand
  exons_pos <- exons[exons$strand=="+",]
  exons_neg <- exons[exons$strand=="-",]

  #make vectors of all co-ordinates within the branchpoint window
  intron_locs_p <- unlist(lapply(exons_pos$start, function(x){(x-max_dist):(x-1)}))
  intron_chroms_p <- rep(exons_pos$chromosome, each=max_dist)
  intron_locs_n <- unlist(lapply(exons_neg$end, function(x){(x+1):(x+max_dist)}))
  intron_chroms_n <- rep(exons_neg$chromosome, each=max_dist)

  introns <- data.frame(chromosome=c(intron_chroms_p,intron_chroms_n),
                        chr_start=c(intron_locs_p,intron_locs_n))

  #format as chrom_1000000
  chrom_and_start <- paste(introns$chromosome, introns$chr_start,sep="_")
  chrom_and_start <- gsub(" ","",chrom_and_start)

  return(chrom_and_start)
}
