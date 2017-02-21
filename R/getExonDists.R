#' Get the closest 3' and 5' exons
#'
#' Finds the closest annotated exons from a genomic co-ordinate.
#' Returns the distance to the 3' exon, distance to the 5' exon,
#' ids of the 3' and 5' exon,
#' and if the exons are from the same parent gene
#' @param queryLine line from a query data.frame (or a character vector)
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and starnd at position 4.
#' @param exonAnnotation data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param queryType type of query. "SNP" or "region"
#' @return vector with distance to the closest 3' and 5' exons,
#' whether these exons are part of the same gene (i.e is the location intronic),
#' and the identifiers for the 3' and 5' exons.
#' @export
#' @import GenomicRanges
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query_line <- readQueryFile(query_snp,query_type = "SNP")[1,]
#' exons_dists <- getExonDists(query_line, exons, query_type = "SNP")
#' @author Beth Signal

getExonDists <- function(queryLine, exonAnnotation, queryType){
  
  # if using regions, region start will be used for negative strand queries
  # & region end for positive strand queries
  
  #faster when exonAnnotation is filtered
  exonAnnotation.subset <- exonAnnotation[(exonAnnotation@seqnames == 
                                           as.character(queryLine@seqnames) & 
                                           exonAnnotation@strand == 
                                           queryLine@strand)]
  
  queryLine2 <- queryLine
  if(as.logical(queryLine2@strand == '-')){
    queryLine2@ranges@width  = as.integer(1)
  }else{
    end <- queryLine2@ranges@width - 1 + queryLine2@ranges@start
    queryLine2@ranges@width  = as.integer(1)
    queryLine2@ranges@start = as.integer(end)
  }
  
  # follow to 5'
  f <- follow(queryLine2, exonAnnotation.subset)
  gene5 <- exonAnnotation.subset[f]@elementMetadata$gene_id
  exon5 <- exonAnnotation.subset[f]@elementMetadata$exon_id
  to5prime <- distance(queryLine2,exonAnnotation.subset[f]) + 1
  
  # preceed to 3'
  p <- precede(queryLine2, exonAnnotation.subset)
  gene3 <- exonAnnotation.subset[p]@elementMetadata$gene_id
  exon3 <- exonAnnotation.subset[p]@elementMetadata$exon_id
  to3prime <- distance(queryLine2,exonAnnotation.subset[p]) + 1
  
  sameGene <- gene3 == gene5
  
  queryLine@elementMetadata$to_3prime <- to3prime
  queryLine@elementMetadata$to_5prime <- to5prime
  queryLine@elementMetadata$same_gene <- sameGene
  queryLine@elementMetadata$exon_3prime <- exon3
  queryLine@elementMetadata$exon_5prime <- exon5
  
  return(queryLine)
}
