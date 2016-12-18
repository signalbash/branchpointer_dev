#' Find the closest 3' and 5' exons to a branchpointer query
#'
#' Finds the closest annotated exons from genomic co-ordinates in a branchpointer query data.frame.

#' @param query branchpointer query data.frame
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and strand at position 4.
#' @param query_type type of query file (\code{"SNP"} or \code{"region"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param max_dist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons.
#' @param use_parallel use parallelisation to speed up code?
#' (reccomended for > 500 query entries)
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @return data.frame with the query and its location relative to the 3' and 5' exons
#' @export
#' @import parallel
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#'
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,query_type = "SNP")
#' query <- getQueryLoc(query,query_type = "SNP",exons = exons, filter = FALSE)
#' @author Beth Signal

getQueryLoc <- function(query, query_type,max_dist=50, filter=TRUE, exons,
                        use_parallel=FALSE, cores=1){

  if(missing(query_type) | !(query_type %in% c("SNP", "region"))){

    stop("please specify query_type as \"region\" or \"SNP\"")

  }

  if(missing(exons)){

    stop("please specify exon annotation data.frame")

  }

  if(use_parallel){

    max_cores <- parallel::detectCores()

    if(max_cores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", max_cores,")"))
      message(paste0("using all available cores"))
      cores=max_cores

    }

  }

  if(query_type=="SNP"){

      if(filter){
        message("filtering for SNPs in branchpoint windows")
        introns=exonsToIntrons(exons, max_dist)
        query_2=paste(query$chromosome, query$chrom_start, sep="_")
        query_2=gsub(" ","", query_2)
        x=match(query_2, introns)
        rm=which(is.na(x))
        if(length(rm)>0){
          query <- query[-rm,]
        }
      }
    }

  if(use_parallel){
    cluster <- makeCluster(cores)
    nearest_exons <- parApply(cluster,query,1,getExonDists,exons, query_type)
    stopCluster(cluster)
  }else{
    nearest_exons <- apply(query,1,getExonDists,exons,query_type)
  }

  nearest_exons <- matrix(unlist(nearest_exons), ncol=5, byrow = TRUE)
  nearest_exons <- as.data.frame((nearest_exons), stringsAsFactors=FALSE)
  rownames(nearest_exons) <- query$id
  colnames(nearest_exons) <- c("to_3prime","to_5prime","same_gene","exon_3prime","exon_5prime")
  nearest_exons$to_3prime <- as.numeric(nearest_exons$to_3prime)
  nearest_exons$to_5prime <- as.numeric(nearest_exons$to_5prime)

  #remove any points not in same gene
  rm <- which((!is.na(nearest_exons$same_gene) & nearest_exons$same_gene == FALSE) |
    is.na(nearest_exons$exon_3prime) | is.na(nearest_exons$exon_5prime))

  if(length(rm) > 0){

    query <- query[-rm,]
    nearest_exons <- nearest_exons[-rm,]

  }


  if(query_type=="region"){
    #if introns regions are closer to the 3'SS than the branchpoint region
    # move the region to cover the window of the nearest exon
    near_3 <- which(nearest_exons$to_3prime < 18)
    if(length(near_3) > 0){
      add_distance <- 18 - nearest_exons$to_3prime[near_3]
      #positive strand
      query$chrom_end[near_3][(query$strand[near_3] == "+")] <-
        query$chrom_end[near_3][(query$strand[near_3] == "+")] -
        add_distance[(query$strand[near_3] == "+")]
      #negative strand
      query$chrom_start[near_3][(query$strand[near_3] == "-")] <-
        query$chrom_start[near_3][(query$strand[near_3] == "-")] +
        add_distance[(query$strand[near_3] == "-")]
    }

    # remove instances where the query region is to close to the exon
    start_to_end <- query$chrom_start <= query$chrom_end
    if(any(!start_to_end)){
      rm <- which(start_to_end == FALSE)
      query <- query[-rm,]
      nearest_exons <- nearest_exons[-rm,]
    }

    #if introns regions are further from the 3'SS than the branchpoint region
    not_near_3 <- which(nearest_exons$to_3prime > 44)
    if(length(not_near_3) > 0){
      query <- query[-not_near_3,]
      nearest_exons <- nearest_exons[-not_near_3,]
    }

    #get exon distances for re-aligned windows
    if(use_parallel){
      cluster <- makeCluster(cores)
      nearest_exons <- parApply(cluster,query,1,getExonDists,exons,query_type)
      stopCluster(cluster)
    }else{
      nearest_exons <- apply(query,1,getExonDists,exons,query_type)
    }

    nearest_exons <- as.data.frame(t(nearest_exons),stringsAsFactors = FALSE)
    rownames(nearest_exons) <- query$id
    colnames(nearest_exons) <- c("to_3prime","to_5prime","same_gene","exon_3prime","exon_5prime")
    nearest_exons$to_3prime <- as.numeric(nearest_exons$to_3prime)
    nearest_exons$to_5prime <- as.numeric(nearest_exons$to_5prime)


    #adjust query location to only cover the 27nt window
    move <- 18 - nearest_exons$to_3prime

    query_region_bp <- query
    query_region_bp$chrom_end[query_region_bp$strand == "+"] <-
      query_region_bp$chrom_end[query_region_bp$strand == "+"] -
      move[query_region_bp$strand == "+"]
    query_region_bp$chrom_start[query_region_bp$strand == "-"] <-
      query_region_bp$chrom_start[query_region_bp$strand == "-"] +
      move[query_region_bp$strand == "-"]

    new_start <- query_region_bp$chrom_end[query_region_bp$strand == "+"] -26
    replace_start <- query_region_bp$chrom_start[query_region_bp$strand == "+"] < new_start
    query_region_bp$chrom_start[query_region_bp$strand == "+"][replace_start] <- new_start[replace_start]

    new_start <- query_region_bp$chrom_start[query_region_bp$strand == "-"] +26
    replace_start <- query_region_bp$chrom_end[query_region_bp$strand == "-"] > new_start
    query_region_bp$chrom_end[query_region_bp$strand == "-"][replace_start] <- new_start[replace_start]

    queryLoc=cbind(query_region_bp, nearest_exons, stringsAsFactors=FALSE)
  }else if(query_type=="SNP"){

    #if a 5' or 3' exon can't be found remove from analysis
    if(any(nearest_exons$to_3prime==-1 | nearest_exons$to_5prime==-1)){
      rm <- which((nearest_exons$to_3prime==-1 | nearest_exons$to_5prime==-1))
      query <- query[-rm,]
      nearest_exons <- nearest_exons[-rm,]
    }

    if(any(nearest_exons$to_3prime > max_dist)){
      rm <- which(nearest_exons$to_3prime >max_dist)
      query <- query[-rm,]
      nearest_exons <- nearest_exons[-rm,]
    }

    #check both exons come from the same parent gene (i.e query is in intronic region)
    #fix 1/0 TRUE/FALSE encoding
    if(any(nearest_exons$same_gene==0)){
      rm <- which(nearest_exons$same_gene==0)
      query <- query[-rm,]
      nearest_exons <- nearest_exons[-rm,]
    }
    queryLoc <- cbind(query, nearest_exons)
  }
  return(queryLoc)
}
