#' Convert SNP branchpoint predictions across the branchpoint window to an intronic summary
#'
#' Takes predictions of branchpoint probabilities from a reference
#' and alternative SNP and summarises the effect(s) of the SNP.
#' @param predictions site-wide branchpoint proability predictions
#' produced from predictBranchpoints()
#' @param query query data.frame containing all snp ids to be summarised
#' @param BP_probability_cutoff Value to be used as the cutoff for
#' discriminating branchpoint sites from non-branchpoint sites
#' (default = \code{0.5})
#' @param BP_probability_change Minimum probability score change
#' required to call a branchpoint site as deleted or created by
#' a SNP (default = \code{0.2})
#' @return data.frame with summarised branchpoint changes
#' occuring within the intron due to a SNP.
#' @export
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,query_type = "SNP")
#' query <- getQueryLoc(query,query_type = "SNP",exons = exons, filter = FALSE)
#' query_attributes <- getBranchpointSequence(query,
#' query_type = "SNP",
#' genome = "~/Downloads/GRCh38.p5.genome.fa",
#' bedtools_location = "/Applications/apps/bedtools2/bin/bedtools")
#' branchpoint_predictions <- predictBranchpoints(query_attributes)
#' snp_stats <- predictionsToStats(branchpoint_predictions, query)
#' @author Beth Signal
#'
predictionsToStats <- function(predictions, query,
                               BP_probability_cutoff = 0.5,
                               BP_probability_change = 0.2){

  snp_ids <- query$id

  branchpoint_snp_stats <- data.frame(id = snp_ids,
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

    for(z in seq(along = snp_ids)){
      #match snp id to predictions
      branchpoint_mutated <- predictions[which(!is.na(match(predictions$id, snp_ids[z]))),]

      #match snp id to query attributes
      info <- query[which(!is.na(match(query$id, snp_ids[z]))),]

      #find location of the SNP relative to the annotated 3'exon
      if(min(abs(branchpoint_mutated$end - info$chrom_start)) == 0){
        branchpoint_snp_stats$dist_to_exon[z] <- branchpoint_mutated$distance[which.min(abs(branchpoint_mutated$end - info$chrom_start))]
      }else{
        if(info$strand == "+"){
          exon_start <- branchpoint_mutated$end[1] + branchpoint_mutated$distance[1]
          branchpoint_snp_stats$dist_to_exon[z] <- exon_start - info$chrom_start
        }else{
          exon_start <- branchpoint_mutated$end[1] - branchpoint_mutated$distance[1]
          branchpoint_snp_stats$dist_to_exon[z] <- info$chrom_start - exon_start
        }
      }

      # find location of the SNP relative to the predicted BPs
      # in reference and alternative sequences
      diffsREF <- branchpoint_snp_stats$dist_to_exon[z] -
        branchpoint_mutated$distance[branchpoint_mutated$allele_status=="REF"][which(
        branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="REF"] >= BP_probability_cutoff)]

      branchpoint_snp_stats$BP_num_REF[z] <- length(diffsREF)

      if(branchpoint_snp_stats$BP_num_REF[z] > 0){
        branchpoint_snp_stats$dist_to_BP_REF[z] <- diffsREF[which.min(abs(diffsREF))]
      }else{
        branchpoint_snp_stats$dist_to_BP_REF[z] <- NA
      }

      diffsALT <- branchpoint_snp_stats$dist_to_exon[z]-branchpoint_mutated$distance[branchpoint_mutated$allele_status=="ALT"][which(
        branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"] >= BP_probability_cutoff)]

      branchpoint_snp_stats$BP_num_ALT[z] <- length(diffsALT)

      if(branchpoint_snp_stats$BP_num_ALT[z] >0){
        branchpoint_snp_stats$dist_to_BP_ALT[z] <- diffsALT[which.min(abs(diffsALT))]
      }else{
        branchpoint_snp_stats$dist_to_BP_ALT[z] <- NA
      }

      # maximum branchpoint probability score
      BP_Norm <- which(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status == "REF"] > BP_probability_cutoff)
      branchpoint_snp_stats$max_prob_REF[z] <- max(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status == "REF"])

      # maximum U2 binding energy of all branchpoint sites
      # above the probability cutoff
      if(branchpoint_snp_stats$BP_num_REF[z] > 0){
        branchpoint_snp_stats$max_U2_REF[z] <- max(
          branchpoint_mutated$U2_binding_energy[branchpoint_mutated$allele_status == "REF"][BP_Norm])
      }else{
        branchpoint_snp_stats$max_U2_REF[z] <- NA
      }

      # maximum branchpoint probability score
      BP_Mut <- which(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"] > BP_probability_cutoff)
      branchpoint_snp_stats$max_prob_ALT[z] <- max(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"])

      # maximum U2 binding energy of all branchpoint sites
      # above the probability cutoff
      if(branchpoint_snp_stats$BP_num_ALT[z] > 0){
        branchpoint_snp_stats$max_U2_ALT[z] <- max(branchpoint_mutated$U2_binding_energy[branchpoint_mutated$allele_status=="ALT"][BP_Mut] )
      }else{
        branchpoint_snp_stats$max_U2_ALT[z] <- NA
      }

      # all sites (ref & alt) called as branchpoints
      dists <- unique(branchpoint_mutated$distance[branchpoint_mutated$branchpoint_prob > BP_probability_cutoff])

      #probabiltiy scores for ref/alt at all BP sites
      ref_p <- branchpoint_mutated[branchpoint_mutated$allele_status=="REF" & branchpoint_mutated$distance %in% dists,]$branchpoint_prob
      alt_p <- branchpoint_mutated[branchpoint_mutated$allele_status=="ALT" & branchpoint_mutated$distance %in% dists,]$branchpoint_prob

      #change must be of sufficient magnitude to be called as created or deleted
      branchpoint_snp_stats$deleted_n[z] = length(which((ref_p - alt_p) > BP_probability_change &
                                                          (ref_p < BP_probability_cutoff | alt_p < BP_probability_cutoff)))
      branchpoint_snp_stats$created_n[z] = length(which((ref_p - alt_p) < (BP_probability_change * -1) &
                                                          (ref_p < BP_probability_cutoff | alt_p < BP_probability_cutoff)))
    }

  m <- match(snp_ids, query$id)
  n <- which(colnames(query) %in% c("id","chromosome","chrom_start","strand",
                                    "ref_allele","alt_allele") & !duplicated(colnames(query)))

  branchpoint_snp_stats <- cbind(query[m,n],
                                 branchpoint_snp_stats[,-1])

  return(branchpoint_snp_stats)
}
