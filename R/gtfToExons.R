#' Convert GTF file to exon location file
#'
#' Converts a GTF annotation to exon locations
#' @param gtf file containing the gtf annotation.
#' @param output_file name of the file to write the exon annotation to.
#' If not specified, will write to the same location as the gtf file.
#' @return exon annotation data.frame
#' @export
#' @import data.table
#' @import stringr
#' @importFrom utils write.table
#' @examples
#' exons <- gtfToExons("gencode.v24.annotation.gtf", output_file="gencode.v24.exons.txt")
#' @author Beth Signal

gtfToExons <- function(gtf, output_file=NA){

  gtf <- as.data.frame(data.table::fread(gtf))

  #keep exon annotations only
  gtf <- gtf[gtf$V3 == "exon",]

  if(grepl("biotype", gtf$V9[1])){
    type_name <- "biotype"
  }else{
    type_name <- "type"
  }

  gtf$gene_id <- getGtfAttribute(gtf, "gene_id")
  gtf$gene_biotype <- getGtfAttribute(gtf, paste0("gene_", type_name))
  gtf$transcript_id <- getGtfAttribute(gtf, "transcript_id")
  gtf$transcript_biotype <- getGtfAttribute(gtf, paste0("transcript_", type_name))

  if(any(grepl("exon", gtf$V9))){
    gtf$exon_number <- getGtfAttribute(gtf, "exon_number")
    gtf$exon_id <- getGtfAttribute(gtf, "exon_id")
  }else{
    gtf$exon_number <- NA
    gtf$exon_id <- gtf$transcript_id

    #make exon names
    n <- 1
    while(any(is.na(gtf$exon_number))){
      transcript_ids <- unique(gtf$transcript_id)
      m <- match(transcript_ids, gtf$transcript_id[which(is.na(gtf$exon_number))])
      gtf$exon_number[which(is.na(gtf$exon_number))[m]] <- n
      n <- n+1
    }
  }

  if(!is.na(output_file)){
  utils::write.table(gtf[,c(1,3,4,5,7,10:15)], file=output_file,
              quote=FALSE, col.names =FALSE, row.names = FALSE, sep ="\t")
  }

  return(gtf[,c(1,3,4,5,7,10:15)])
}
