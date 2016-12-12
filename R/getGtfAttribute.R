#' extract a field from the GTF annotation
#'
#' extract a field from the GTF annotation (usually column 9)
#' @param gtf gtf annotation data.frame
#' @param attribute name of the attribute to exact the value for
#' @return vector with attribute value, or NA if attribute name is not found.
#' @export
#' @import stringr
#' @examples
#' gtf <- as.data.frame(data.table::fread("~/Downloads/gencode.v24.annotation.gtf"))
#' gtf$gene_id <- getGtfAttribute(gtf, "gene_id")
#' @author Beth Signal

getGtfAttribute <-  function(gtf, attribute){

  values <- unlist(stringr::str_split(gtf$V9, ";"))

  values <- grep(attribute, values, value = TRUE)

  values <- gsub(paste0(attribute," "), "", values)

  values <- gsub(" ", "", values)

  values <- gsub('"', "", values)

  if(length(values) == length(gtf[,9])){
    return(values)
  }else{
    return(NA)
  }

}
