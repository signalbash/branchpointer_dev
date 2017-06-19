#' Convert exon annotation GRanges to intron locations
#'
#' Converts exon annotation to intron locations overlapping the branchpoint region
#' for exculsion of non-branchpoint region SNPs
#' Returns a character vector of chromosome locations
#' @param exons GRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist Maximum distance from the 3' exon to create the branchpoint region.
#' @return GRanges containing intron window co-ordinates
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal
exonsToIntrons <- function(exons, maxDist = 50){
  
  introns <- exons
  
  posInd <- which(as.logical(strand(exons) == "+"))
  start(introns[posInd]) <- start(introns[posInd]) - (maxDist + 1)
  end(introns[posInd]) <- start(introns[posInd]) + maxDist
  
  negInd <- which(as.logical(strand(exons) == "-"))
  end(introns[negInd]) <- end(introns[negInd]) + (maxDist + 1)
  start(introns[negInd]) <- end(introns[negInd]) - maxDist
  
  return(introns)
}

#' Get the closest 3' and 5' exons
#'
#' Finds the closest annotated exons from a genomic co-ordinate.
#' Returns the distance to the 3' exon, distance to the 5' exon,
#' ids of the 3' and 5' exon,
#' and if the exons are from the same parent gene
#' @param query GenomicRangesquery
#' @param exons GenomicRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param queryType type of query. "SNP" or "region"
#' @return GenomicRanges with distance to the closest 3' and 5' exons,
#' whether these exons are part of the same gene (i.e is the location intronic),
#' and the identifiers for the 3' and 5' exons.
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal

getExonDists <- function(query, exons, queryType){
  
  # if using regions, region start will be used for negative strand queries
  # & region end for positive strand queries

  negInd <- which(as.logical(strand(query) == '-'))
  posInd <- which(as.logical(strand(query) == '+'))
  
  queryPoint <- query
  
  width(ranges(queryPoint))[negInd] <- 1
  end <- end(ranges(queryPoint))[posInd]
  start(ranges(queryPoint[posInd])) <- end
  width(ranges(queryPoint[posInd])) <- 1
  
  f <- GenomicRanges::follow(queryPoint,exons)
  gene5 <- exons$gene_id[f]
  exon5 <- exons$exon_id[f]
  
  keep <- which(!is.na(f))
  f <- f[keep]
  to5prime <- rep(NA, length(queryPoint))
  to5prime[keep] <- GenomicRanges::distance(queryPoint[keep],exons[f]) + 1
  
  
  p <- GenomicRanges::precede(queryPoint, exons)
  gene3 <- exons$gene_id[p]
  exon3 <- exons$exon_id[p]
  
  keep <- which(!is.na(p))
  p <- p[keep]
  to3prime <- rep(NA, length(queryPoint))
  to3prime[keep] <- GenomicRanges::distance(queryPoint[keep],exons[p]) + 1
  
  sameGene <- gene3 == gene5
  
  query$to_3prime <- to3prime
  query$to_5prime <- to5prime
  query$same_gene <- sameGene
  query$exon_3prime <- exon3
  query$exon_5prime <- exon5
  
  return(query)
}

#' Find the closest 3' and 5' exons to a branchpointer query
#'
#' Finds the closest annotated exons from genomic co-ordinates in a branchpointer query GRanges
#' @param query branchpointer query GenomicRanges
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and strand at position 4.
#' @param queryType type of query file (\code{"SNP"} or \code{"region"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons.
#' @return GenomicRanges with the query and its location relative to the 3' and 5' exons
#' @import GenomicRanges
#' @author Beth Signal

getQueryLoc <- function(query, queryType,maxDist=50, filter=TRUE, exons){

  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){

    stop("please specify queryType as \"region\" or \"SNP\"")

  }

  if(missing(exons)){

    stop("please specify exon annotation object")

  }

  if(queryType=="SNP"){

      if(filter){
        message("filtering for SNPs in branchpoint windows")
        introns <- exonsToIntrons(exons, maxDist)
        o <- GenomicRanges::findOverlaps(query, introns)
        keep <- unique(o@from)
        query <- query[keep]

      }
  }
  
  nearestExons <- getExonDists(query, exons, queryType)

  #remove any points not in same gene
  rm <- which(nearestExons$same_gene == FALSE |
    is.na(nearestExons$exon_3prime) | 
    is.na(nearestExons$exon_5prime))

  if(length(rm) > 0){
    nearestExons <- nearestExons[-rm]
  }


  if(queryType=="region"){
    #if introns regions are closer to the 3'SS than the branchpoint region
    # move the region to cover the window of the nearest exon
    near3 <- which(nearestExons$to_3prime < 18)
    if(length(near3) > 0){
      addDistance <- 18 - nearestExons$to_3prime[near3]
      negStrand <- which(as.logical(strand(nearestExons)[near3] == "-"))
      posStrand <- which(as.logical(strand(nearestExons)[near3] == "+"))
      
      if(length(negStrand) > 0){
        start(ranges(nearestExons))[near3][negStrand] <- 
          start(ranges(nearestExons))[near3][negStrand] + addDistance[negStrand]
      }
      if(length(posStrand) > 0){
          end(ranges(nearestExons))[near3][posStrand] <- 
              end(ranges(nearestExons))[near3][posStrand] - addDistance[posStrand]
      }
      
    }

    #if introns regions are further from the 3'SS than the branchpoint region
    notNear3 <- which(nearestExons$to_3prime > 44)
    if(length(notNear3) > 0){
      nearestExons <- nearestExons[-notNear3]
    }

    nearestExons <- getExonDists(nearestExons,exons,queryType)
    
    #adjust query location to only cover the 27nt window
    
    #negStrand -- move end
    negStrand <- which(as.logical(strand(nearestExons) == "-"))
    if(length(negStrand) > 0){
        end(ranges(nearestExons[negStrand])) <- 
            start(ranges(nearestExons[negStrand])) - 
            nearestExons$to_3prime[negStrand] + 44
    }
    #posStrand -- move start
    posStrand <- which(as.logical(strand(nearestExons) == "+"))
    if(length(posStrand) > 0){
        start(ranges(nearestExons[posStrand])) <- 
            end(ranges(nearestExons[posStrand])) +
            nearestExons$to_3prime[posStrand] - 44
    }

    }else if(queryType=="SNP"){

      #if a 5' or 3' exon can't be found remove from analysis
      rm <- which(nearestExons$to_3prime==-1 | nearestExons$to_5prime==-1)
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
      
      rm <- which(nearestExons$to_3prime > maxDist)
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
  
      #check both exons come from the same parent gene 
      #(i.e query is in intronic region)
      rm <- which(!(nearestExons$same_gene))
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
  }
  return(nearestExons)
}
