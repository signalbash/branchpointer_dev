#' Convert exon annotation data.frame to intron locations
#'
#' Converts exon annotation to intron locations overlapping the branchpoint region
#' for exculsion of non-branchpoint region SNPs
#' Returns a character vector of chromosome_location
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist Maximum distance from the 3' exon to create the branchpoint region.
#' @return vector of chromosome names and intronic locations
#' @export
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt", package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#' introns <- exonsToIntrons(exons, 50)
#' @author Beth Signal

exonsToIntrons <- function(exons, maxDist = 50){

  #split exon annotation by strand
  exons.pos <- exons[exons$strand=="+",]
  exons.neg <- exons[exons$strand=="-",]

  #make vectors of all co-ordinates within the branchpoint window
  intronLocs.p <- unlist(lapply(exons.pos$start, function(x){(x-maxDist):(x-1)}))
  intronChroms.p <- rep(exons.pos$chromosome, each=maxDist)
  intronLocs.n <- unlist(lapply(exons.neg$end, function(x){(x+1):(x+maxDist)}))
  intronChroms.n <- rep(exons.neg$chromosome, each=maxDist)

  introns <- data.frame(chromosome=c(intronChroms.p,intronChroms.n),
                        chr_start=c(intronLocs.p,intronLocs.n))

  #format as chrom_1000000
  chromAndStart <- paste(introns$chromosome, introns$chr_start,sep="_")
  chromAnd_Start <- gsub(" ","",chromAndStart)

  return(chromAndStart)
}
