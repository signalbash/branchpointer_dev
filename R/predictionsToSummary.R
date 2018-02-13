#' Convert SNP branchpoint predictions across the branchpoint window to an intronic summary
#'
#' Takes predictions of branchpoint probabilities from a reference
#' and alternative SNP and summarises the effect(s) of the SNP.
#' @param predictions site-wide branchpoint proability predictions
#' produced from predictBranchpoints()
#' @param query query GRanges containing all SNP ids to be summarised
#' @param probabilityCutoff Value to be used as the cutoff for
#' discriminating branchpoint sites from non-branchpoint sites
#' (default = \code{0.52})
#' @param probabilityChange Minimum probability score change
#' required to call a branchpoint site as deleted or created by
#' a SNP (default = \code{0.15})
#' @return GRanges with summarised branchpoint changes
#' occuring within the intron due to a SNP.
#' @export
#' @import GenomicRanges
#' @importFrom stats aggregate
#' @examples
#' smallExons <- system.file("extdata","gencode.v26.annotation.small.gtf", package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' 
#' querySNPFile <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' querySNP <- readQueryFile(querySNPFile,queryType = "SNP",exons = exons, filter = FALSE)
#' predictionsSNP <- predictBranchpoints(querySNP,queryType = "SNP",BSgenome = g)
#' 
#' summarySNP <- predictionsToSummary(querySNP,predictionsSNP)
#' @author Beth Signal
#'
predictionsToSummary <- function(query,
                                 predictions,
                                 probabilityCutoff = 0.52,
                                 probabilityChange = 0.15){
  
  mcols(query)$BP_num_REF = 0
  mcols(query)$BP_num_ALT = 0
  mcols(query)$deleted_n = 0
  mcols(query)$created_n = 0
  mcols(query)$is_in_window = 0
  mcols(query)$dist_to_BP_REF = NA
  mcols(query)$dist_to_BP_ALT = NA
  mcols(query)$max_prob_REF = NA
  mcols(query)$max_prob_ALT = NA
  mcols(query)$max_U2_REF = NA
  mcols(query)$max_U2_ALT = NA
  
  predictions <- as.data.frame(mcols(predictions))

  # Reference sequence summary
  bps.ref <- which(predictions$branchpoint_prob >= probabilityCutoff & 
                     predictions$status == "REF")
  bps.ref <- predictions[bps.ref,]
  bps.ref$distance <- bps.ref$to_3prime - bps.ref$to_3prime_point
  
  tab <- as.data.frame(table(bps.ref$id))
  m <- match(query$id, tab$Var1)
  query$BP_num_REF[which(!is.na(m))] <- tab$Freq[m[which(!is.na(m))]]
  
  bps.ref <- bps.ref[order(abs(bps.ref[,which(colnames(bps.ref) == "distance")])),]
  m <- match(query$id, bps.ref$id)
  query$dist_to_BP_REF[which(!is.na(m))] <- bps.ref$distance[m[which(!is.na(m))]]
  
  bps.ref.all <- predictions[which(predictions$status == "REF"),]
  bps.ref.all <- bps.ref.all[rev(order(bps.ref.all[,which(colnames(bps.ref.all) == "branchpoint_prob")])),]
  
  m <- match(query$id, bps.ref.all$id)
  query$max_prob_REF[which(!is.na(m))] <- bps.ref.all$branchpoint_prob[m[which(!is.na(m))]]
  
  bps.ref <- bps.ref[rev(order(bps.ref[,which(colnames(bps.ref) == "U2_binding_energy")])),]
  
  m <- match(query$id, bps.ref$id)
  query$max_U2_REF[which(!is.na(m))] <- bps.ref$U2_binding_energy[m[which(!is.na(m))]]
  
  # Alternative sequence summary
  bps.alt <- which(predictions$branchpoint_prob >= probabilityCutoff 
                   & predictions$status == "ALT")
  bps.alt <- predictions[bps.alt,]
  bps.alt$distance <- bps.alt$to_3prime - bps.alt$to_3prime_point
  
  tab <- as.data.frame(table(bps.alt$id))
  m <- match(query$id, tab$Var1)
  query$BP_num_ALT[which(!is.na(m))] <- tab$Freq[m[which(!is.na(m))]]
  
  bps.alt <- bps.alt[order(abs(bps.alt[,which(colnames(bps.alt) == "distance")])),]
  
  m <- match(query$id, bps.alt$id)
  query$dist_to_BP_ALT[which(!is.na(m))] <- bps.alt$distance[m[which(!is.na(m))]]
  
  bps.alt.all <- predictions[which(predictions$status == "ALT"),]
  bps.alt.all <- bps.alt.all[rev(order(bps.alt.all[,which(colnames(bps.alt.all) == "branchpoint_prob")])),]
  
  m <- match(query$id, bps.alt.all$id)
  query$max_prob_ALT[which(!is.na(m))] <- bps.alt.all$branchpoint_prob[m[which(!is.na(m))]]
  
  bps.alt <- bps.alt[rev(order(bps.alt[,which(colnames(bps.alt) == "U2_binding_energy")])),]
  
  m <- match(query$id, bps.alt$id)
  query$max_U2_ALT[which(!is.na(m))] <- bps.alt$U2_binding_energy[m[which(!is.na(m))]]
  
  
  #?indel
  shift <- ifelse(query$ref_allele =="-", 0, nchar(query$ref_allele)) - 
      ifelse(query$alt_allele =="-", 0, nchar(query$alt_allele))
  
  #id.loc <- with(predictions, paste0(id,"_", to_3prime_point))
  id.loc <- with(predictions, paste0(id,"_", test_site))
  
  notInAltWindow <- bps.ref$id[which(
      bps.ref$to_3prime_point - bps.ref$extend > 44 | 
          bps.ref$to_3prime_point - bps.ref$extend < 18)]
  
  notInAltWindow <- as.data.frame(table(notInAltWindow))
  m <- match(notInAltWindow$notInAltWindow, query$id)
  query$is_in_window[m] <- notInAltWindow$Freq
  
  
  notInAltWindow <- bps.ref$id[which(
      bps.ref$to_3prime_point - bps.ref$extend > 44 | 
          bps.ref$to_3prime_point - bps.ref$extend < 18)]
  
  notInAltWindow <- as.data.frame(table(notInAltWindow))
  m <- match(notInAltWindow$notInAltWindow, query$id)
  query$is_in_window[m] <- notInAltWindow$Freq
  
  
  notInRefWindow <- bps.alt$id[which(
      bps.alt$to_3prime_point - bps.alt$extend > 44 | 
          bps.alt$to_3prime_point - bps.alt$extend < 18)]
  
  notInRefWindow <- as.data.frame(table(notInRefWindow))
  m <- match(notInRefWindow$notInRefWindow, query$id)
  query$is_in_window[m] <- query$is_in_window[m] + notInRefWindow$Freq
  
  index <- which(id.loc %in% id.loc[predictions$branchpoint_prob >
                                        probabilityCutoff])
  
  # Any sites changed?
  bps.all <- predictions[index,]
  bps.all <- bps.all[,c('id','test_site','to_3prime_point',
                        'branchpoint_prob','status')]
  
  #bps.all$id_site <- paste0(bps.all$id,"_",bps.all$to_3prime_point)
  bps.all$id_site <- paste0(bps.all$id,"_",bps.all$test_site)
  
  bps.ref <- bps.all[bps.all$status=="REF",]
  bps.alt <- bps.all[bps.all$status=="ALT",]
  
  
  bps.ref$prob_ref <- bps.ref$branchpoint_prob
  bps.alt$prob_ref <- bps.ref$branchpoint_prob[match(bps.alt$id_site,
                                                     bps.ref$id_site)]
  
  bps.alt$prob_alt <- bps.alt$branchpoint_prob
  bps.ref$prob_alt <- bps.alt$branchpoint_prob[match(bps.ref$id_site,
                                                     bps.alt$id_site)]
  
  bps.all <- rbind(bps.ref, bps.alt[which(!(bps.alt$id_site %in% 
                                                bps.ref$id_site)),])
  
  deleted <- as.data.frame(table(bps.ref$id[which(
      bps.all$prob_ref > probabilityCutoff &
          ((bps.all$prob_ref - bps.all$prob_alt > probabilityCutoff &
          bps.all$prob_alt < probabilityCutoff) | 
          (is.na(bps.all$prob_alt)))
  )]))
  
  created <- as.data.frame(table(bps.ref$id[which(
      bps.all$prob_alt > probabilityCutoff &
          ((bps.all$prob_alt - bps.all$prob_ref > probabilityCutoff &
                bps.all$prob_ref < probabilityCutoff) | 
               (is.na(bps.all$prob_ref)))
  )]))
  
  m <- match(query$id, deleted$Var1)
  query$deleted_n[which(!is.na(m))] <- deleted$Freq[m[which(!is.na(m))]]
  m <- match(query$id, created$Var1)
  query$created_n[which(!is.na(m))] <- created$Freq[m[which(!is.na(m))]]
  
  return(query) 
}
