#' Convert exon annotation GRanges to intron locations
#'
#' Converts exon annotation to intron locations overlapping the branchpoint region
#' for exculsion of non-branchpoint region SNPs
#' Returns a character vector of chromosome locations
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist Maximum distance from the 3' exon to create the branchpoint region.
#' @return vector of chromosome names and intronic locations
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal
exonsToIntrons <- function(exons, maxDist = 50){
  
  #split exon annotation by strand
  exons.pos <- exons[exons@strand=="+",]
  exons.neg <- exons[exons@strand=="-",]
  
  #make vectors of all co-ordinates within the branchpoint window
  intronLocs.p <- unlist(lapply(exons.pos@ranges@start, function(x){(x-maxDist):(x-1)}))
  intronChroms.p <- rep(exons.pos@seqnames, each=maxDist)
  intronLocs.n <- unlist(lapply(exons.neg@ranges@start + exons.neg@ranges@width - 1, function(x){(x+1):(x+maxDist)}))
  intronChroms.n <- rep(exons.neg@seqnames, each=maxDist)
  
  introns <- data.frame(chromosome=c(intronChroms.p,intronChroms.n),
                        chr_start=c(intronLocs.p,intronLocs.n))
  
  #format as chrom_1000000
  chromAndStart <- paste(introns$chromosome, introns$chr_start,sep="_")
  chromAndStart <- gsub(" ","",chromAndStart)
  
  return(chromAndStart)
}

#' Get the closest 3' and 5' exons
#'
#' Finds the closest annotated exons from a genomic co-ordinate.
#' Returns the distance to the 3' exon, distance to the 5' exon,
#' ids of the 3' and 5' exon,
#' and if the exons are from the same parent gene
#' @param queryLine line from a query GenomicRanges
#' @param exons GenomicRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param queryType type of query. "SNP" or "region"
#' @return GenomicRanges with distance to the closest 3' and 5' exons,
#' whether these exons are part of the same gene (i.e is the location intronic),
#' and the identifiers for the 3' and 5' exons.
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal

getExonDists <- function(queryLine, exons, queryType){
  
  # if using regions, region start will be used for negative strand queries
  # & region end for positive strand queries
  
  #faster when exonAnnotation is filtered
  exons.subset <- exons[(exons@seqnames == 
                           as.character(queryLine@seqnames) & 
                           exons@strand == queryLine@strand)]
  
  queryLine2 <- queryLine
  if(as.logical(queryLine2@strand == '-')){
    queryLine2@ranges@width  = as.integer(1)
  }else{
    end <- queryLine2@ranges@width - 1 + queryLine2@ranges@start
    queryLine2@ranges@width  = as.integer(1)
    queryLine2@ranges@start = as.integer(end)
  }
  
  # follow to 5'
  f <- GenomicRanges::follow(queryLine2, exons.subset)
  if(length(f) > 0){
    gene5 <- exons.subset[f]@elementMetadata$gene_id
    exon5 <- exons.subset[f]@elementMetadata$exon_id
    to5prime <- GenomicRanges::distance(queryLine2,exons.subset[f]) + 1
  }else{
    gene5 <- NA
    exon5 <- NA
    to5prime <- NA
  }
  
  # preceed to 3'
  p <- GenomicRanges::precede(queryLine2, exons.subset)
  if(length(p) > 0){
    gene3 <- exons.subset[p]@elementMetadata$gene_id
    exon3 <- exons.subset[p]@elementMetadata$exon_id
    to3prime <- GenomicRanges::distance(queryLine2,exons.subset[p]) + 1
  }else{
    gene3 <- NA
    exon3 <- NA
    to3prime <- NA
  }
  
  sameGene <- gene3 == gene5
  
  queryLine@elementMetadata$to_3prime <- to3prime
  queryLine@elementMetadata$to_5prime <- to5prime
  queryLine@elementMetadata$same_gene <- sameGene
  queryLine@elementMetadata$exon_3prime <- exon3
  queryLine@elementMetadata$exon_5prime <- exon5
  
  return(queryLine)
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
#' @param useParallel use parallelisation to speed up code?
#' (reccomended for > 500 query entries)
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @return GenomicRanges with the query and its location relative to the 3' and 5' exons
#' @export
#' @import parallel
#' @import GenomicRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#'
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,queryType = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' @author Beth Signal

getQueryLoc <- function(query, queryType,maxDist=50, filter=TRUE, exons,
                        useParallel=FALSE, cores=1){

  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){

    stop("please specify queryType as \"region\" or \"SNP\"")

  }

  if(missing(exons)){

    stop("please specify exon annotation object")

  }

  if(useParallel){

    maxCores <- parallel::detectCores()

    if(maxCores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", maxCores,")"))
      message(paste0("using all available cores"))
      cores <- maxCores

    }

  }

  if(queryType=="SNP"){

      if(filter){
        message("filtering for SNPs in branchpoint windows")
        introns <- exonsToIntrons(exons, maxDist)
        query.2 <- paste(query@seqnames, query@ranges@start, sep="_")
        query.2 <- gsub(" ","", query.2)
        x <- match(query.2, introns)
        rm <- which(is.na(x))
        if(length(rm)>0){
          query <- query[-rm,]
        }
      }
    }

  if(useParallel){
    cluster <- makeCluster(cores)
    nearestExons <- parApply(cluster,query,1,getExonDists,exons, queryType)
    stopCluster(cluster)
  }else{
    nearestExons <- lapply(query,getExonDists,exons,queryType)
  }

  nearestExons <- do.call("c", nearestExons)
  
  #nearestExons$to_3prime <- as.numeric(nearestExons$to_3prime)
  #nearestExons$to_5prime <- as.numeric(nearestExons$to_5prime)

  #remove any points not in same gene
  rm <- which(nearestExons@elementMetadata$same_gene == FALSE |
    is.na(nearestExons@elementMetadata$exon_3prime) | 
    is.na(nearestExons@elementMetadata$exon_5prime))

  if(length(rm) > 0){
    nearestExons <- nearestExons[-rm]
  }


  if(queryType=="region"){
    #if introns regions are closer to the 3'SS than the branchpoint region
    # move the region to cover the window of the nearest exon
    near3 <- which(nearestExons$to_3prime < 18)
    if(length(near3) > 0){
      addDistance <- 18 - nearestExons$to_3prime[near3]
      negStrand <- which(as.logical(nearestExons@strand[near3] == "-"))
      
      if(length(negStrand) > 0){
        nearestExons@ranges@start[near3][negStrand] <- nearestExons@ranges@start[near3][negStrand] + as.integer(addDistance[negStrand])
      }
      nearestExons@ranges@width[near3] <- nearestExons@ranges@width[near3] - as.integer(addDistance)
    }

    # remove instances where the query region is too close to the exon
    #startToEnd <- query$chrom_start <= query$chrom_end
    #if(any(!startToEnd)){
    #  rm <- which(startTo_End == FALSE)
    #  query <- query[-rm,]
    #  nearestExons <- nearestExons[-rm,]
    #}

    #if introns regions are further from the 3'SS than the branchpoint region
    notNear3 <- which(nearestExons@elementMetadata$to_3prime > 44)
    if(length(notNear3) > 0){
      nearestExons <- nearestExons[-notNear3]
    }

    #get exon distances for re-aligned windows
    if(useParallel){
      cluster <- makeCluster(cores)
      nearestExons <- parApply(cluster,query,1,getExonDists,exons,queryType)
      stopCluster(cluster)
    }else{
      nearestExons <- lapply(nearestExons,getExonDists,exons,queryType)
    }
    nearestExons <- do.call("c", nearestExons)

    #adjust query location to only cover the 27nt window
    move <- 18 - nearestExons@elementMetadata$to_3prime
    
    nearestExons <- nearestExons

    #negStrand -- move start
    negStrand <- which(as.logical(nearestExons@strand == "-"))
    nearestExons@ranges@start[negStrand] <- 
      nearestExons@ranges@start[negStrand] + as.integer(move)[negStrand]
    nearestExons@ranges@width <- nearestExons@ranges@width - as.integer(move)
    
    #posStrand
    posStrand <- which(as.logical(nearestExons@strand == "+"))
    end <- (nearestExons@ranges@start + nearestExons@ranges@width - 1)[posStrand]
    nearestExons@ranges@start[posStrand] <- as.integer(end - 26)
    
    #setting multiple GRanges values to a single value requires indexing
    nearestExons@ranges@width[c(posStrand, negStrand)] <- as.integer(27)
    
    }else if(queryType=="SNP"){

      #if a 5' or 3' exon can't be found remove from analysis
      rm <- which(nearestExons@elementMetadata$to_3prime==-1 | nearestExons@elementMetadata$to_5prime==-1)
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
      
      rm <- which(nearestExons@elementMetadata$to_3prime > maxDist)
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
  
      #check both exons come from the same parent gene (i.e query is in intronic region)
      rm <- which(!(nearestExons@elementMetadata$same_gene))
      if(length(rm) > 0){
        nearestExons <- nearestExons[-rm,]
      }
  }
  return(nearestExons)
}
