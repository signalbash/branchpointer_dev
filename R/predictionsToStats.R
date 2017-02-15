#' Convert SNP branchpoint predictions across the branchpoint window to an intronic summary
#'
#' Takes predictions of branchpoint probabilities from a reference
#' and alternative SNP and summarises the effect(s) of the SNP.
#' @param predictions site-wide branchpoint proability predictions
#' produced from predictBranchpoints()
#' @param query query data.frame containing all snp ids to be summarised
#' @param probabilityCutoff Value to be used as the cutoff for
#' discriminating branchpoint sites from non-branchpoint sites
#' (default = \code{0.5})
#' @param probabilityChange Minimum probability score change
#' required to call a branchpoint site as deleted or created by
#' a SNP (default = \code{0.2})
#' @return data.frame with summarised branchpoint changes
#' occuring within the intron due to a SNP.
#' @export
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,query_type = "SNP")
#' query <- getQueryLoc(query,query_type = "SNP",exons = exons, filter = FALSE)
#' query_attributes <- getBranchpointSequence(query,
#' query_type = "SNP",
#' useBSgenome = TRUE,
#' BSgenome = genome)
#' branchpoint_predictions <- predictBranchpoints(query_attributes)
#' snp_stats <- predictionsToStats(branchpoint_predictions, query)
#' @author Beth Signal
#'
predictionsToStats <- function(predictions, query,
                               probabilityCutoff = 0.5,
                               probabilityChange = 0.2){

  snpIDs <- query$id

  branchpointSNPstats <- data.frame(id = snpIDs,
                                      BP_num_REF = NA,
                                      BP_num_ALT = NA,
                                      deleted_n = NA,
                                      created_n = NA,
                                      dist_to_exon = NA,
                                      dist_to_BP_REF = NA,
                                      dist_to_BP_ALT = NA,
                                      max_prob_REF = NA,
                                      max_prob_ALT = NA,
                                      max_U2_REF = NA,
                                      max_U2_ALT = NA)

    for(z in seq(along = snpIDs)){
      #match snp id to predictions
      branchpoint.mutated <- predictions[which(!is.na(match(predictions$id, snp_ids[z]))),]

      #match snp id to query attributes
      info <- query[which(!is.na(match(query$id, snpIDs[z]))),]

      #find location of the SNP relative to the annotated 3'exon
      if(min(abs(branchpoint.mutated$end - info$chrom_start)) == 0){
        branchpointSNPstats$dist_to_exon[z] <- branchpoint.mutated$distance[which.min(abs(branchpoint.mutated$end - info$chrom_start))]
      }else{
        if(info$strand == "+"){
          exonStart <- branchpoint.mutated$end[1] + branchpoint.mutated$distance[1]
          branchpointSNPstats$dist_to_exon[z] <- exonStart - info$chrom_start
        }else{
          exonStart <- branchpoint_mutated$end[1] - branchpoint.mutated$distance[1]
          branchpointSNPstats$dist_to_exon[z] <- info$chrom_start - exonStart
        }
      }

      # find location of the SNP relative to the predicted BPs
      # in reference and alternative sequences
      diffs.ref <- branchpointSNPstats$dist_to_exon[z] -
        branchpoint.mutated$distance[branchpoint.mutated$allele_status=="REF"][which(
        branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status=="REF"] >= probabilityCutoff)]

      branchpointSNPstats$BP_num_REF[z] <- length(diffs.ref)

      if(branchpointSNPstats$BP_num_REF[z] > 0){
        branchpointSNPstats$dist_to_BP_REF[z] <- diffs.ref[which.min(abs(diffs.ref))]
      }else{
        branchpointSNPstats$dist_to_BP_REF[z] <- NA
      }

      diffs.alt <- branchpointSNPstats$dist_to_exon[z]-branchpoint.mutated$distance[branchpoint.mutated$allele_status=="ALT"][which(
        branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status=="ALT"] >= probabilityCutoff)]

      branchpointSNPstats$BP_num_ALT[z] <- length(diffs.alt)

      if(branchpointSNPstats$BP_num_ALT[z] >0){
        branchpointSNPstats$dist_to_BP_ALT[z] <- diffs.alt[which.min(abs(diffs.alt))]
      }else{
        branchpointSNPstats$dist_to_BP_ALT[z] <- NA
      }

      # maximum branchpoint probability score
      BP.ref <- which(branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status == "REF"] > probabilityCutoff)
      branchpointSNPstats$max_prob_REF[z] <- max(branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status == "REF"])

      # maximum U2 binding energy of all branchpoint sites
      # above the probability cutoff
      if(branchpointSNPstats$BP_num_REF[z] > 0){
        branchpointSNPstats$max_U2_REF[z] <- max(
          branchpoint.mutated$U2_binding_energy[branchpoint.mutated$allele_status == "REF"][BP.ref])
      }else{
        branchpointSNPstats$max_U2_REF[z] <- NA
      }

      # maximum branchpoint probability score
      BP.alt <- which(branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status=="ALT"] > probabilityCutoff)
      branchpointSNPstats$max_prob_ALT[z] <- max(branchpoint.mutated$branchpoint_prob[branchpoint.mutated$allele_status=="ALT"])

      # maximum U2 binding energy of all branchpoint sites
      # above the probability cutoff
      if(branchpointSNPstats$BP_num_ALT[z] > 0){
        branchpointSNPstats$max_U2_ALT[z] <- max(branchpoint.mutated$U2_binding_energy[branchpoint.mutated$allele_status=="ALT"][BP.alt] )
      }else{
        branchpointSNPstats$max_U2_ALT[z] <- NA
      }

      # all sites (ref & alt) called as branchpoints
      dists <- unique(branchpoint.mutated$distance[branchpoint.mutated$branchpoint_prob > probabilityCutoff])

      #probabiltiy scores for ref/alt at all BP sites
      prob.ref <- branchpoint.mutated[branchpoint.mutated$allele_status=="REF" & branchpoint.mutated$distance %in% dists,]$branchpoint_prob
      prob.alt <- branchpoint.mutated[branchpoint.mutated$allele_status=="ALT" & branchpoint.mutated$distance %in% dists,]$branchpoint_prob

      #change must be of sufficient magnitude to be called as created or deleted
      branchpointSNPstats$deleted_n[z] = length(which((prob.ref - prob.alt) > probabilityChange &
                                                          (prob.ref < probabilityCutoff | prob.alt < probabilityCutoff)))
      branchpointSNPstats$created_n[z] = length(which((prob.ref - prob.alt) < (probabilityChange * -1) &
                                                          (prob.ref < probabilityCutoff | prob.alt < probabilityCutoff)))
    }

  m <- match(snpIDs, query$id)
  n <- which(colnames(query) %in% c("id","chromosome","chrom_start","strand",
                                    "ref_allele","alt_allele") & !duplicated(colnames(query)))

  branchpointSNPstats <- cbind(query[m,n],
                                 branchpointSNPstats[,-1])

  return(branchpointSNPstats)
}
