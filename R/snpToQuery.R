#' Gets a branchpointer formatted query from refsnp ids
#'
#' Searches Biomart for refsnp ids, and pulls genomic location and sequence identity information
#' Reformats alleles so each query has only one alternative allele
#' @param refSNP Vector of refsnp ids
#' @param mart.snp biomaRt mart object specifying the BioMart database and dataset to be used
#' @return formatted SNP query data.frame
#' @export
#' @import biomaRt
#' @importFrom stringr str_split
#' @examples
#' mart <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="www.ensembl.org")
#' query <- snpToQuery("rs17000647", mart)
#' @author Beth Signal

snpToQuery <- function(refSNP, mart.snp) {
  snpInfo <- biomaRt::getBM(attributes = c("refsnp_id",'refsnp_source', "chr_name",
                                   "chrom_start", "allele"),
                    filters = "snp_filter", values = refSNP, mart = mart.snp)

  #make sure each SNP has only 1 ref and 1 alternate allele
  multiAlleles <- which(nchar(snpInfo$allele) != 3)

  if (length(multiAlleles) > 0) {
    snpInfo.remade <- snpInfo[-multiAlleles,]
    snpInfo <- snpInfo[multiAlleles,]

    for (i in seq_along(snpInfo$refsnp_id)) {
      nts <- unlist(stringr::str_split(snpInfo$allele[i],"/"))
      ref <- nts[1]
      alt <- nts[-1]
      alleles <- paste0(ref,"/",alt)

      remade <- snpInfo[c(rep(i, length(alleles))),]
      remade$allele <- alleles
      snpInfo.remade <- rbind(snpInfo.remade, remade)
    }
    snpInfo <- snpInfo.remade
  }

  snpInfoQuery = data.frame(
    id = snp_info$refsnp_id, chromosome = paste0("chr",snpInfo$chr_name),
    chrom_start = snpInfo$chrom_start, strand = 2,
    ref_allele = str_sub(snpInfo$allele,1,1),
    alt_allele = str_sub(snpInfo$allele,3,3)
  )


  #check for unstranded queries & replace with positive & negative
  unstranded <- which(snpInfoQuery$strand != "+" | snpInfoQuery$strand !="-")
  if(length(unstranded)>0){
    query.pos <- snpInfoQuery[unstranded,]
    query.pos$id <- paste0(query.pos$id, ".pos")
    query.pos$strand <- "+"
    query.neg <- snpInfoQuery[unstranded,]
    query.neg$id <- paste0(query.neg$id, ".neg")
    query.neg$strand <- "-"
    snpInfoQuery <- rbind(snpInfoQuery[-unstranded,], query.pos, query.neg)
  }

  #check for duplicated query ids
  if(any(duplicated(snpInfoQuery$id))){
    message(paste0(length(which(duplicated(snpInfoQuery$id)))," query ids are not unique"))
    message("Check output for new names or rename")
    snpInfoQuery$id <- make.names(snpInfoQuery$id, unique=TRUE)
  }

  return(snpInfoQuery)

}
