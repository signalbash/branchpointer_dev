#' Read a query file
#'
#' Reads and formats a query file, generated from snpToQuery, or manually generated.
#' Converts unstranded SNPs to two entries for each strand.
#' Checks for duplicate names and replaces if found.
#' @param queryFile tab delimited file containing query information.
#' For intronic regions should be in the format: region id, chromosome name, region start, region id, strand.
#' For SNP variants should be in the format: SNP id, chromosome name, SNP position, strand, reference allele (A/T/C/G), alternative allele (A/T/C/G)
#' @param queryType type of query file (\code{"SNP"} or \code{"region"})
#' @return data.frame with formatted query
#' @export
#' @import data.table
#' @examples
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,query_type = "SNP")
#'
#' query_snp <- system.file("extdata","intron_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,query_type = "region")
#' @author Beth Signal

readQueryFile <- function(queryFile, queryType){

  if(missing(queryType)){

    queryTest <- fread(
      queryFile,header = TRUE,
      nrows=1
      )

    if(ncol(queryTest) == 5){

      query_type <- "region"
      message("first line of file has 5 columns")
      message("using query_type = \"region\"")

    }else if(ncol(queryTest) == 6){

      query_type <- "SNP"
      message("first line of file has 6 columns")
      message("using query_type = \"SNP\"")

    }else{

      stop("please specify query_type and provide correctly formatted file")

    }

  }

  if (!(query_type %in% c("SNP", "region"))) {

    message("please specify query_type as \"region\" or \"SNP\"")

  }else{

    if (query_type == "SNP") {

      query <- fread(
        queryFile,header = TRUE,
        colClasses = c(
          "character", "character", "numeric",
          "character", "character", "character"
        )
      )
      query <- as.data.frame(query)
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

    }else if (queryType == "region") {

      query <- fread(
        queryFile, header = TRUE,colClasses = c(
          "character","character","numeric","numeric","character"
        )
      )
      query <- as.data.frame(query)
      colnames(query) <-
        c("id","chromosome","chrom_start","chrom_end","strand")

    }

    #check for duplicated query ids
    if (any(duplicated(query$id))) {

      message(paste0(length(which(
        duplicated(query$id)
      ))," query ids are not unique"))
      message("Check output for new names or rename")
      query$id = make.names(query$id, unique = TRUE)

    }

    return(query)

  }

}
