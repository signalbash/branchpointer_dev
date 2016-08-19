#' Gets a branchpointer formatted query from refsnp ids
#'
#' Searches Biomart for refsnp ids, and pulls genomic location and sequence identity information
#' Reformats alleles so each query has only one alternative allele
#' @param rs_id Vector of refsnp ids
#' @param mart_snp biomaRt mart object specifying the BioMart database and dataset to be used
#' @return formatted SNP query data.frame
#' @export
#' @import biomaRt
#' @importFrom stringr str_split
#' @examples
#' mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="www.ensembl.org")
#' query <- snpToQuery("rs17000647", mart)
#' @author Beth Signal

snpToQuery <- function(rs_id, mart_snp) {
  snp_info <- biomaRt::getBM(attributes = c("refsnp_id",'refsnp_source', "chr_name",
                                   "chrom_start", "allele"),
                    filters = "snp_filter", values = rs_id, mart = mart_snp)

  #make sure each SNP has only 1 ref and 1 alternate allele
  multi_alleles <- which(nchar(snp_info$allele) != 3)

  if (length(multi_alleles) > 0) {
    snp_info_remade <- snp_info[-multi_alleles,]
    snp_info <- snp_info[multi_alleles,]

    for (i in seq_along(snp_info$refsnp_id)) {
      nts <- unlist(stringr::str_split(snp_info$allele[i],"/"))
      ref <- nts[1]
      alt <- nts[-1]
      alleles <- paste0(ref,"/",alt)

      remade <- snp_info[c(rep(i, length(alleles))),]
      remade$allele <- alleles
      snp_info_remade <- rbind(snp_info_remade, remade)
    }
    snp_info <- snp_info_remade
  }

  snp_info_query = data.frame(
    id = snp_info$refsnp_id, chromosome = paste0("chr",snp_info$chr_name),
    chrom_start = snp_info$chrom_start, strand = 2,
    ref_allele = str_sub(snp_info$allele,1,1),
    alt_allele = str_sub(snp_info$allele,3,3)
  )


  #check for unstranded queries & replace with positive & negative
  unstranded <- which(snp_info_query$strand != "+" | snp_info_query$strand !="-")
  if(length(unstranded)>0){
    query_pos <- snp_info_query[unstranded,]
    query_pos$id <- paste0(query_pos$id, "_pos")
    query_pos$strand <- "+"
    query_neg <- snp_info_query[unstranded,]
    query_neg$id <- paste0(query_neg$id, "_neg")
    query_neg$strand <- "-"
    snp_info_query <- rbind(snp_info_query[-unstranded,], query_pos,query_neg)
  }

  #check for duplicated query ids
  if(any(duplicated(snp_info_query$id))){
    message(paste0(length(which(duplicated(snp_info_query$id)))," query ids are not unique"))
    message("Check output for new names or rename")
    snp_info_query$id <- make.names(snp_info_query$id, unique=TRUE)
  }

  return(snp_info_query)

}
