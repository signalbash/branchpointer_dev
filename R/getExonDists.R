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
  if(queryType=="region"){
    query.start <- as.numeric(queryLine[3])
    query.end <- as.numeric(queryLine[4])
    query.chrom <- as.character(queryLine[2])
    query.strand <- as.character(queryLine[5])
  }else if(queryType=="SNP"){
    query.start <- as.numeric(queryLine[3])
    query.end <- as.numeric(queryLine[3])
    query.chrom <- as.character(queryLine[2])
    query.strand <- as.character(queryLine[4])
  }

  #match strand and chromosome
  exonAnnotation.subset <- exonAnnotation[which(exonAnnotation$strand == query.strand &
                                                    exonAnnotation$chromosome == query.chrom),]
  if(query.strand == "+"){

    #get distance to from query point to all exon starts (3' exon distance)
    #need to remove exon # == 1 for 3' exons
    to3prime <- exonAnnotation.subset$start[which(exonAnnotation.subset$exon_number > 1)] - query.end

    #get gene and exon id of nearest exon
    ind3prime <- which.min(to3prime[to3prime > 0])

    if(length(ind3prime) == 0){

      gene3=NA
      exon3=NA

    } else {

      gene3 <- exonAnnotation.subset$gene_id[which(exonAnnotation.subset$exon_number > 1)][which(to3prime >0)[ind3prime]]
      exon3 <- exonAnnotation.subset$exon_id[which(exonAnnotation.subset$exon_number > 1)][which(to3prime >0)[ind3prime]]

    }

    #get minimum distance
    if(all(to3prime <= 0)){
      to3prime <- -1
    }else{
      to3prime <- min(to3prime[to3prime > 0])
    }

    #get distance to from query point to all exon ends (5' exon distance)
    to5prime <- query.end - exonAnnotation.subset$end

    #get gene and exon id of nearest exon
    ind5prime <- which.min(to5prime[to5prime > 0])

    if(length(ind5prime) == 0){

      gene5=NA
      exon5=NA

    } else {

      gene5 <- exonAnnotation.subset$gene_id[which(to5prime >0)[ind5prime]]
      exon5 <- exonAnnotation.subset$exon_id[which(to5prime >0)[ind5prime]]

    }
    #get minimum distance
    if(all(to5prime <= 0)){
      to5prime <- -1
    }else{
      to5prime <- min(to5prime[to5prime > 0])
    }

  }else{
    #get distance to from query point to all exon ends (3' exon distance)
    to3prime <- query.start - exonAnnotation.subset$end[which(exonAnnotation.subset$exon_number > 1)]

    #get gene and exon id of nearest exon
    ind3prime <- which.min(to3prime[to3prime > 0])

    if(length(ind3prime) == 0){

      gene3=NA
      exon3=NA

    } else {

      gene3 <- exonAnnotation.subset$gene_id[which(exonAnnotation.subset$exon_number > 1)][which(to3prime >0)[ind3prime]]
      exon3 <- exonAnnotation.subset$exon_id[which(exonAnnotation.subset$exon_number > 1)][which(to3prime >0)[ind3prime]]

    }

    #get minimum distance
    if(all(to3prime <= 0)){
      to3prime <- -1
    }else{
      to3prime <- min(to3prime[to3prime > 0])
    }
    #get distance to from query point to all exon starts (5' exon distance)
    to5prime <- exonAnnotation.subset$start - query.start

    #get gene and exon id of nearest exon
    ind5prime <- which.min(to5prime[to5prime > 0])

    if(length(ind5prime) == 0){

      gene5=NA
      exon5=NA

    } else {

      gene5 <- exonAnnotation.subset$gene_id[which(to5prime >0)[ind5prime]]
      exon5 <- exonAnnotation.subset$exon_id[which(to5prime >0)[ind5prime]]

    }
    #get minimum distance
    if(all(to5prime <= 0)){
      to5prime <- -1
    }else{
      to5prime <- min(to5prime[to5prime > 0])
    }
  }

  #check if the 3' and 5' exons are from the same gene
  #(i.e. if the co-ordinate is intronic)
  sameGene <- gene3==gene5
  return(c(to3prime,to5prime, sameGene, exon3, exon5))
}
