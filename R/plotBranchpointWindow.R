#' Plots branchpointer predictions
#'
#' Plots branchpointer predictions
#' @param queryName query id used to identify the SNP or region
#' @param queryAttributes attribute data.frame
#' @param predictions prediction of branchpoint probabilities for each site
#' @param probabilityCutoff probability score cutoff value for displaying U2 binding energy
#' @param plotStructure plot structures for gene and 3' exon containing and skipping isoforms
#' @param plotMutated plot alternative sequence predicitons alongside reference sequence predictions
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return ggplot2 plot with branchpoint features in the specified intronic region
#' @export
#' @import ggplot2
#' @import ggbio
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
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
#' plotBranchpointWindow(query$id[1],  branchpoint_predictions,
#' query_attributes, plot_mutated = TRUE, exons = exons)
#' @author Beth Signal

plotBranchpointWindow <- function(queryName,
                                 predictions,
                                 probabilityCutoff = 0.5,
                                 plotMutated = FALSE,
                                 plotStructure = TRUE,
                                 exons) {
  #find the query
  ind <- which(!is.na(match(
    predictions@elementMetadata$id, queryName
  )))
  predictions.snp <- predictions[which(predictions@elementMetadata$id == queryName)]

  predictions.snp@elementMetadata$to_3prime_point <- predictions.snp@elementMetadata$to_3prime_point * -1
  
  predictions.snp.ref <- as.data.frame(predictions.snp@elementMetadata[predictions.snp@elementMetadata$status == "REF",])
  colnames(predictions.snp.ref) <- gsub("to_3prime_point","distance",colnames(predictions.snp.ref))
  colnames(predictions.snp.ref) <- gsub("seq_pos0","nucleotide",colnames(predictions.snp.ref))
  predictions.snp.ref$nucleotide <- as.character(predictions.snp.ref$nucleotide)
  
  #make plots for reference sequence
  
  #probability score by sequence position
  plot.prob.ref <- ggplot(predictions.snp.ref, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                           col = 'nucleotide',alpha = 'U2_binding_energy')) +
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                       breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                          "0.60","0.70","0.80","0.90","1.00")) +
    scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                      limits = c("A","C","G","T")) +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  #U2 binding energy by sequence position
  plot.U2.ref <- ggplot(predictions.snp.ref[predictions.snp.ref$branchpoint_prob >= probabilityCutoff,],
              aes_string(x = 'distance', y = 'U2_binding_energy')) +
    geom_bar(stat = "identity",width = 1) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    theme(legend.position = "none") +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

  #Sequence identity
  plot.seq.ref <- ggplot(predictions.snp.ref, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
    geom_text(size = 4, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  if (plotMutated == TRUE) {
    
    predictions.snp.alt <- as.data.frame(predictions.snp@elementMetadata[predictions.snp@elementMetadata$status == "ALT",])
    colnames(predictions.snp.alt) <- gsub("to_3prime_point","distance",colnames(predictions.snp.alt))
    colnames(predictions.snp.alt) <- gsub("seq_pos0","nucleotide",colnames(predictions.snp.alt))
    predictions.snp.alt$nucleotide <- as.character(predictions.snp.alt$nucleotide)
    
    #probability score by sequence position
    plot.prob.alt <- ggplot(predictions.snp.alt, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                             col = 'nucleotide',alpha = 'U2_binding_energy')) +
      geom_bar(stat = "identity") +
      scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                         breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                            "0.60","0.70","0.80","0.90","1.00")) +
      scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                        limits = c("A","C","G","T")) +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = c("A","C","G","T")) +
      theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17)) + ggtitle("Alternative")

    plot.prob.ref <- plot.prob.ref + ggtitle("Reference")

    #U2 binding energy by sequence position
    plot.U2.alt <- ggplot(predictions.snp.alt[predictions.snp.alt$branchpoint_prob >= probabilityCutoff,],
                aes_string(x = 'distance', y = 'U2_binding_energy')) +
      geom_bar(stat = "identity",width = 1) +
      scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                         breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
      theme(legend.position = "none") +
      scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                         breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

    #Sequence identity
    plot.seq.alt <- ggplot(predictions.snp.alt, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = c("A","C","G","T")) +
      theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17))

    ######Get reference vs. alt sequence for 0-50 window#####

  #   seq <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
  #                                    grepl("REF",queryAttributes$id))[1],]$seq
  #   pos <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
  #                                    grepl("REF",queryAttributes$id))[1],]$to_3prime
  #   seq.ref <- unlist(stringr::str_split(stringr::str_sub(seq,251 + (pos - 50),250 + pos),""))
  # 
  #   seq <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
  #                                    grepl("ALT",queryAttributes$id))[1],]$seq
  #   pos <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
  #                                    grepl("ALT",queryAttributes$id))[1],]$to_3prime
  #   seq.alt <- unlist(stringr::str_split(stringr::str_sub(seq,251 + (pos - 50),250 + pos),""))
  # 
  #   wholeSeq <- data.frame(nt = c(seq_1, seq_2),
  #                          pos = rep(seq(-50,-1,1),2),
  #                          set = rep(c("REF","ALT"), each =50))
  #   mutPosition <- wholeSeq$pos[which(wholeSeq$nt[wholeSeq$set == "REF"] !=
  #                                   wholeSeq$nt[wholeSeq$set == "ALT"])]
  # 
  #   plot.seqComparison <- ggplot(wholeSeq, aes_string(x = 'pos', y = 'set', col = 'nt',label = 'nt')) +
  #     geom_rect(xmin = mutPosition - 0.5,xmax = mutPosition + 0.5, ymin = 0.5, ymax = 0.5 + 2,
  #               col = NA, fill = "grey90") +
  #     geom_text(size = 4, family = "Courier") +
  #     scale_color_manual(values = mercer_nt_cols) +
  #     theme(legend.position = "none",
  #           axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
  #           panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
  #           axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  #     scale_x_continuous(limits = c(-50,-1)) +
  #     scale_y_discrete(labels = c("Alternative", "Reference")) +
  #     geom_rect(xmin = -44.5,xmax = -17.5, ymin = 0.5, ymax = 0.5 + 2, fill = NA, col ="black")
  # 
  }

  #plot gene/transcript structure
  if (plotStructure == TRUE) {
    
    exonInd <- which(exons@elementMetadata$exon_id == predictions.snp@elementMetadata$exon_3prime[1])
    geneID <- exons@elementMetadata$gene_id[exonInd][1]
    
    exonsForPlot <- exons[exons@elementMetadata$gene_id == geneID]
    exonsForPlot@elementMetadata$tested_exon <- "X"
    
    
    if(as.logical(exonsForPlot@strand[1] == "+")){
      exon5primeLoc <- exons@ranges@start[exonInd][1]
      exonsForPlot@elementMetadata$tested_exon[which(exonsForPlot@ranges@start == exon5primeLoc)] <- "Y"
    }else{
      exon5primeLoc <- exons@ranges@start[exonInd][1] + exons@ranges@width[exonInd][1] - 1
      exonsForPlot@elementMetadata$tested_exon[which((exonsForPlot@ranges@start + exonsForPlot@ranges@width - 1) == exon5primeLoc)] <- "Y"
      
    }
     
    followTranscripts <- unique(exonsForPlot@elementMetadata$transcript_id[which(!is.na(follow(exonsForPlot,exons[exonInd][1])))])
    precedeTranscripts <- unique(exonsForPlot@elementMetadata$transcript_id[which(!is.na(precede(exonsForPlot,exons[exonInd][1])))])
    
    keepTranscripts <- followTranscripts[followTranscripts %in% precedeTranscripts]
    exonsForPlot <- exonsForPlot[exonsForPlot@elementMetadata$transcript_id %in% keepTranscripts]
    
    plot.structure <-  ggplot() + ggbio::geom_alignment(exonsForPlot, aes(fill=tested_exon)) + 
      scale_fill_manual(values = c("black","grey60")) +
      theme(panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title.y = element_blank(), axis.ticks.y = element_blank(),
            legend.position = "none")
    
  }
    
  theme_set(theme_gray())
  if (plotStructure == TRUE & plotMutated == TRUE) {
    ggdraw() +
      draw_plot(plot.structure,0,0.775,1,0.225) + 
      #draw_plot(plot.seq.comparison,0,0.675,1,0.1) +
      draw_plot(plot.prob.ref,0,.275,0.5,.40) + draw_plot(plot.seq.ref,0,.2,0.5,.075) + draw_plot(plot.U2.ref,0,0,0.5,.2) +
      draw_plot(plot.prob.alt,0.5,.275,0.5,.40) + draw_plot(plot.seq.alt,0.5,.2,0.5,.075) + draw_plot(plot.U2.alt,0.5,0,0.5,.2)

  }else if (plotStructure == TRUE & plotMutated == FALSE) {
    ggdraw() +
      draw_plot(plot.structure,0,0.775,1,0.225) +
      draw_plot(plot.prob.ref,0,.325,1,.45) + draw_plot(plot.seq.ref,0,.25,1,.075) + draw_plot(plot.U2.ref,0,0,1,.25)
  }else if (plotStructure == FALSE & plotMutated == FALSE) {
    ggdraw() +
      draw_plot(plot.prob.ref,0,.325,1,.45) + draw_plot(plot.seq.ref,0,.25,1,.075) + draw_plot(plot.U2.ref,0,0,1,.25)
  }else if (plotStructure == FALSE & plotMutated == TRUE) {
    ggdraw() +
      #draw_plot(plot.seq.comparison,0,0.9,1,0.1) +
      draw_plot(plot.prob.ref,0,.4,0.5,.5) + draw_plot(plot.seq.ref,0,.3,0.5,.1) + draw_plot(plot.U2.ref,0,0,0.5,.3) +
      draw_plot(plot.prob.alt,0.5,.4,0.5,.5) + draw_plot(plot.seq.alt,0.5,.3,0.5,.1) + draw_plot(plot.U2.alt,0.5,0,0.5,.3)
  }

}
