#' Make branchpoint window regions
#'
#' Genrate branchpoint window regions corresponding to annotated exon(s) within a
#' queried gene, transcript or exon id
#' @param id identifier for the query gene/transcript/exon id
#' @param idType type of id to match in the exon annotation file (\code{"gene_id"},
#' \code{"transcript_id"}, or \code{"exon_id"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return Granges with formatted query
#' @export
#' @import GenomicRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' windowquery <- makeRegions("ENSG00000139618", "gene_id", exons)
#' windowquery <- makeRegions("ENST00000357654", "transcript_id", exons)
#' windowquery <- makeRegions("ENSE00003518965", "exon_id", exons)
#' @author Beth Signal

makeRegions <- function(id, idType, exons) {

  validTypes <- c("gene_id", "transcript_id","exon_id")

  #missing or invalid types
  noType <- missing(idType)

  if(!noType){
    noType <-  noType | !(idType %in% validTypes)
  }


  if (!noType) {
    x <- which(colnames(exons@elementMetadata) == idType)
    y <- grep(id, exons@elementMetadata[,x])
  }else{
    y <- vector()
  }

  #go through possible columns if no matches found
  if (length(y) == 0 | noType) {

    idType <- validTypes[1]
    x <- which(colnames(exons@elementMetadata) == idType)
    y <- y <- grep(id, exons@elementMetadata[,x])

    if (length(y) == 0) {

      idType <- validTypes[2]
      x <- which(colnames(exons@elementMetadata) == idType)
      y <- y <- grep(id, exons@elementMetadata[,x])

    }

    if (length(y) == 0) {

      idType <- validTypes[3]
      x <- which(colnames(exons@elementMetadata) == idType)
      y <- y <- grep(id, exons@elementMetadata[,x])

    }

    if (length(y) == 0) {

      stop(paste0("cannot find ", id," in the exon annotation"))

    }

  }

  #use a subset of the exon annotation for faster processing
  if (idType != "gene_id") {
    gene_id <- exons@elementMetadata$gene_id[y[1]]
    y2 <- which(!is.na(match(exons@elementMetadata$gene_id,gene_id)))
  }else{
    y2 <- y
  }

  exons.subset <- exons[y]
  
  #by definition first exons shouldn' have branchpoints
  keep <- which(exons.subset@elementMetadata$exon_number > 1)

  if (as.logical(exons.subset@strand[1] == "+")) {
    windowStarts <- (exons.subset@ranges@start - 50)[keep]
  }else{
    windowStarts <- (exons.subset@ranges@start + exons.subset@ranges@width -1 + 10)[keep]
  }

  window <- exons.subset[keep]
  window@ranges@start <- as.integer(windowStarts)
  window@ranges@width[1:length(windowStarts)] <- as.integer(41)

  return(getQueryLoc(window,queryType = "region",exons = exons[y2]))

}
