#' Get the closest 3' and 5' exons
#'
#' Finds the closest annotated exons from a genomic co-ordinate.
#' Returns the distance to the 3' exon, distance to the 5' exon,
#' ids of the 3' and 5' exon,
#' and if the exons are from the same parent gene
#' @param query_line line from a query data.frame (or a character vector)
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and starnd at position 4.
#' @param exon_annotation data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param query_type type of query. "SNP" or "region"
#' @return vector with distance to the closest 3' and 5' exons,
#' whether these exons are part of the same gene (i.e is the location intronic),
#' and the identifiers for the 3' and 5' exons.
#' @export
#' @examples
#' exons_dists <- getExonDists(query_line, exons)
#' @author Beth Signal

getExonDists <- function(query_line, exon_annotation, query_type){

  # if using regions, region start will be used for negative strand queries
  # & region end for positive strand queries
  if(query_type=="region"){
    query_start <- as.numeric(query_line[3])
    query_end <- as.numeric(query_line[4])
    query_chrom <- as.character(query_line[2])
    query_strand <- as.character(query_line[5])
  }else if(query_type=="SNP"){
    query_start <- as.numeric(query_line[3])
    query_end <- as.numeric(query_line[3])
    query_chrom <- as.character(query_line[2])
    query_strand <- as.character(query_line[4])
  }

  #match strand and chromosome
  exon_annotation.subset <- exon_annotation[which(exon_annotation$strand == query_strand &
                                                    exon_annotation$chromosome == query_chrom),]
  if(query_strand == "+"){

    #get distance to from query point to all exon starts (3' exon distance)
    #need to remove exon # == 1 for 3' exons
    to_3prime <- exon_annotation.subset$start[which(exon_annotation.subset$exon_number > 1)] - query_end

    #get gene and exon id of nearest exon
    ind_3prime <- which.min(to_3prime[to_3prime > 0])

    if(length(ind_3prime) == 0){

      gene_3=NA
      exon_3=NA

    } else {

      gene_3 <- exon_annotation.subset$gene_id[which(exon_annotation.subset$exon_number > 1)][which(to_3prime >0)[ind_3prime]]
      exon_3 <- exon_annotation.subset$exon_id[which(exon_annotation.subset$exon_number > 1)][which(to_3prime >0)[ind_3prime]]

    }

    #get minimum distance
    if(all(to_3prime <= 0)){
      to_3prime <- -1
    }else{
      to_3prime <- min(to_3prime[to_3prime > 0])
    }

    #get distance to from query point to all exon ends (5' exon distance)
    to_5prime <- query_end - exon_annotation.subset$end

    #get gene and exon id of nearest exon
    ind_5prime <- which.min(to_5prime[to_5prime > 0])

    if(length(ind_5prime) == 0){

      gene_5=NA
      exon_5=NA

    } else {

      gene_5 <- exon_annotation.subset$gene_id[which(to_5prime >0)[ind_5prime]]
      exon_5 <- exon_annotation.subset$exon_id[which(to_5prime >0)[ind_5prime]]

    }
    #get minimum distance
    if(all(to_5prime <= 0)){
      to_5prime <- -1
    }else{
      to_5prime <- min(to_5prime[to_5prime > 0])
    }

  }else{
    #get distance to from query point to all exon ends (3' exon distance)
    to_3prime <- query_start - exon_annotation.subset$end[which(exon_annotation.subset$exon_number > 1)]

    #get gene and exon id of nearest exon
    ind_3prime <- which.min(to_3prime[to_3prime > 0])

    if(length(ind_3prime) == 0){

      gene_3=NA
      exon_3=NA

    } else {

      gene_3 <- exon_annotation.subset$gene_id[which(exon_annotation.subset$exon_number > 1)][which(to_3prime >0)[ind_3prime]]
      exon_3 <- exon_annotation.subset$exon_id[which(exon_annotation.subset$exon_number > 1)][which(to_3prime >0)[ind_3prime]]

    }

    #get minimum distance
    if(all(to_3prime <= 0)){
      to_3prime <- -1
    }else{
      to_3prime <- min(to_3prime[to_3prime > 0])
    }
    #get distance to from query point to all exon starts (5' exon distance)
    to_5prime <- exon_annotation.subset$start - query_start

    #get gene and exon id of nearest exon
    ind_5prime <- which.min(to_5prime[to_5prime > 0])

    if(length(ind_5prime) == 0){

      gene_5=NA
      exon_5=NA

    } else {

      gene_5 <- exon_annotation.subset$gene_id[which(to_5prime >0)[ind_5prime]]
      exon_5 <- exon_annotation.subset$exon_id[which(to_5prime >0)[ind_5prime]]

    }
    #get minimum distance
    if(all(to_5prime <= 0)){
      to_5prime <- -1
    }else{
      to_5prime <- min(to_5prime[to_5prime > 0])
    }
  }

  #check if the 3' and 5' exons are from the same gene
  #(i.e. if the co-ordinate is intronic)
  same_gene <- gene_3==gene_5
  return(c(to_3prime,to_5prime, same_gene, exon_3, exon_5))
}
