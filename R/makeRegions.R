#' Make branchpoint window regions
#'
#' Genrate branchpoint window regions corresponding to annotated exon(s) within a
#' queried gene, transcript or exon id
#' @param id identifier for the query gene/transcript/exon id
#' @param id_type type of id to match in the exon annotation file (\code{"gene_id"},
#' \code{"transcript_id"}, or \code{"exon_id"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return data.frame with frormatted query
#' @export
#' @examples
#' windowquery <- makeRegions("ENSG00000139618", "gene_id", exons)
#' windowquery <- makeRegions("ENST00000357654", "transcript_id", exons)
#' windowquery <- makeRegions("ENSE00003518965", "exon_id", exons)
#' @author Beth Signal

makeRegions <- function(id, id_type, exons) {

  valid_types <- c("gene_id", "transcript_id","exon_id")

  #missing or invalid types
  no_type <- missing(id_type)

  if(!no_type){

    no_type <-  no_type | !(id_type %in% valid_types)

  }


  if (!no_type) {

    x <- which(colnames(exons) == id_type)
    y <- grep(id, exons[,x])

  }else{

    y <- vector()

  }

  #go through possible columns if no matches found
  if (length(y) == 0 | no_type) {

    id_type <- valid_types[1]
    x <- which(colnames(exons) == id_type)
    y <- grep(id, exons[,x])

    if (length(y) == 0) {

      id_type <- valid_types[2]
      x <- which(colnames(exons) == id_type)
      y <- grep(id, exons[,x])

    }

    if (length(y) == 0) {

      id_type <- valid_types[3]
      x <- which(colnames(exons) == id_type)
      y <- grep(id, exons[,x])

    }

    if (length(y) == 0) {

      stop(paste0("cannot find ", id," in the exon annotation"))

    }

  }

  #use a subset of the exon annotation for faster processing
  if (id_type != "gene_id") {

    gene_id <- exons$gene_id[y[1]]
    y2 <- which(!is.na(match(exons$gene_id,gene_id)))

  }else{

    y2 <- y

  }

  exons_subset <- exons[y,]

  #by definition first exons shouldn' have branchpoints
  keep <- which(exons_subset$exon_number > 1)

  if (exons_subset$strand[1] == "+") {

    window_starts <- (exons_subset$start - 50)[keep]
    window_ends <- (exons_subset$start - 10)[keep]

  }else{

    window_starts <- (exons_subset$end + 10)[keep]
    window_ends <- (exons_subset$end + 50)[keep]

  }

  window_df <- data.frame(
    id = exons_subset$exon_id[keep],
    chromosome = exons_subset$chromosome[keep],
    chrom_start = window_starts,
    chrom_end = window_ends,
    strand = exons_subset$strand[keep]
  )

  window_df <- window_df[!duplicated(window_df$id),]

  exons_subset <- exons[y2,]

  return(getQueryLoc(window_df,query_type = "region",exons = exons_subset))

}
