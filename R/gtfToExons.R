#' Convert GTF file to exon location file
#'
#' Converts a GTF annotation to exon locations
#' @param gtf file containing the gtf annotation.
#' @return exon annotation data.frame
#' @export
#' @import stringr
#' @import rtracklayer
#' @importFrom utils write.table
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges strand
#' @examples
#' small_gtf <- system.file("extdata","gencode.v24.annotation.small.gtf",
#' package = "branchpointer")
#' exons <- gtfToExons(small_gtf)
#' @author Beth Signal

gtfToExons <- function(gtf){

  gtf.rtrack <- rtracklayer::import(gtf)

  #keep exon annotations only
  gtf.rtrack <- gtf.rtrack[gtf.rtrack@elementMetadata$type=="exon"]

  #change "biotype" to "type" in GTFs from Ensembl
  if(any(grepl("biotype", names(gtf.rtrack@elementMetadata)))){
    names(gtf.rtrack@elementMetadata) <- gsub("biotype",
                                              "type",
                                              names(gtf.rtrack@elementMetadata))
  }

  if(all(!grepl("exon", names(gtf.rtrack@elementMetadata)))){
    gtf.rtrack@elementMetadata$exon_number <- NA
    gtf.rtrack@elementMetadata$exon_id <- gtf.rtrack@elementMetadata$transcript_id

    #make exon names
    n <- 1
    while(any(is.na(gtf.rtrack@elementMetadata$exon_number))){
      transcript_ids <- unique(gtf.rtrack@elementMetadata$transcript_id)
      m <- match(transcript_ids,
                 gtf.rtrack@elementMetadata$transcript_id[which(
                   is.na(gtf.rtrack@elementMetadata$exon_number))])
      gtf.rtrack@elementMetadata$exon_number[which(
        is.na(gtf.rtrack@elementMetadata$exon_number))[m]] <- n
      n <- n+1
    }
  }

  exons <- gtf.rtrack
  exons@elementMetadata <- exons@elementMetadata[,c("gene_id","gene_type",
                                                    "transcript_id","transcript_type",
                                                    "exon_id","exon_number")]
  
  return(exons)
}
