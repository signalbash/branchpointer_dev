#' Read a query file
#'
#' Reads and formats a query file, generated from snpToQuery, or manually generated.
#' Converts unstranded SNPs to two entries for each strand.
#' Checks for duplicate names and replaces if found.
#' @param queryFile tab delimited file containing query information.
#' For intronic regions should be in the format: region id, chromosome name, region start, region id, strand.
#' For SNP variants should be in the format: SNP id, chromosome name, SNP position, strand, reference allele (A/T/C/G), alternative allele (A/T/C/G)
#' @param queryType type of query file (\code{"SNP"} or \code{"region"})
#' @param exons GRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons.
#' @param useParallel use parallelisation to speed up code?
#' (reccomended for > 500 query entries)
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @return Formatted query GRanges
#' @export
#' @import GenomicRanges
#' @importFrom data.table fread
#' @importClassesFrom data.table data.table
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' 
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,queryType = "SNP", exons)
#'
#' queryIntron <- system.file("extdata","intron_example.txt", package = "branchpointer")
#' query <- readQueryFile(queryIntron,queryType = "region", exons)
#' @author Beth Signal

readQueryFile <- function(queryFile, queryType,exons, maxDist=50, filter=TRUE, 
                          useParallel=FALSE, cores=1){

  if(missing(queryType)){
    queryTest <- data.table::fread(
      queryFile,header = TRUE,
      nrows=1)

    if(ncol(queryTest) == 5){
      queryType <- "region"
      message("first line of file has 5 columns")
      message("using query_type = \"region\"")
    }else if(ncol(queryTest) == 6){
      queryType <- "SNP"
      message("first line of file has 6 columns")
      message("using query_type = \"SNP\"")
    }else{
      stop("please specify query_type and provide correctly formatted file")
    }
  }

  if (!(queryType %in% c("SNP", "region"))) {
    message("please specify query_type as \"region\" or \"SNP\"")
  }else{
    if (queryType == "SNP") {
      query <- data.table::fread(
        queryFile,header = TRUE,
        colClasses = c(
          "character", "character", "numeric",
          "character", "character", "character"
        ), data.table=FALSE)
      colnames(query)[1:6] <-
        c("id", "chromosome","chrom_start","strand","ref_allele", "alt_allele")

      #check single nts
      query$ref_allele <- toupper(query$ref_allele)
      query$alt_allele <- toupper(query$alt_allele)

      #check for unstranded queries & replace with positive & negative
      unstranded <- which(query$strand != "+" | query$strand != "-")
      if (length(unstranded) > 0) {
        query.pos <- query[unstranded,]
        query.pos$id <- paste0(query.pos$id, "_pos")
        query.pos$strand <- "+"
        query.neg <- query[unstranded,]
        query.neg$id <- paste0(query.neg$id, "_neg")
        query.neg$strand <- "-"
        query <- rbind(query[-unstranded,], query.pos,query.neg)
        }
      
      queryGRanges <- GRanges(seqnames=S4Vectors::Rle(query$chromosome),
                              ranges=IRanges::IRanges(start=query$chrom_start, width=1),
                              strand=query$strand,
                              id=query$id,
                              ref_allele=query$ref_allele,
                              alt_allele=query$alt_allele)
    }else if (queryType == "region") {
      query <- data.table::fread(
        queryFile, header = TRUE,colClasses = c(
          "character","character","numeric","numeric","character"
        ), data.table=FALSE)
      colnames(query) <-
        c("id","chromosome","chrom_start","chrom_end","strand")
      
      queryGRanges <- GRanges(seqnames=S4Vectors::Rle(query$chromosome),
                              ranges=IRanges::IRanges(start=query$chrom_start, 
                                             end=query$chrom_end),
                              strand=query$strand,
                              id=query$id)
    }

    #check for duplicated query ids
    if (any(duplicated(query$id))) {
      message(paste0(length(which(
        duplicated(query$id)
      ))," query ids are not unique"))
      message("Check output for new names or rename")
      query$id = make.names(query$id, unique = TRUE)
    }
    
    #find 3'/5'exons
    if(length(queryGRanges) > 0){
      queryGRanges.loc <- getQueryLoc(queryGRanges, queryType, maxDist = maxDist, filter = filter,
                                    exons = exons, useParallel = useParallel, cores = cores)
      return(queryGRanges.loc)
    }
    
  }
}
