#' Plots branchpointer predictions
#'
#' Plots branchpointer predictions
#' @param queryName query id used to identify the SNP or region
#' @param queryAttributes attribute data.frame
#' @param branchpointPredictions prediction of branchpoint probabilities for each site
#' @param probabilityCutoff probability score cutoff value for displaying U2 binding energy
#' @param plotStructure plot structures for gene and 3' exon containing and skipping isoforms
#' @param plotMutated plot alternative sequence predicitons alongside reference sequence predictions
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return ggplot2 plot with branchpoint features in the specified intronic region
#' @export
#' @import ggplot2
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
                                 branchpointPredictions,
                                 queryAttributes,
                                 probabilityCutoff = 0.5,
                                 plotMutated = FALSE,
                                 plotStructure = TRUE,
                                 exons) {
  #find the query
  ind <- which(!is.na(match(
    branchpointPredictions$id, queryName
  )))
  branchpoint.mutated <- branchpointPredictions[ind,]

  branchpoint.mutated$distance <- branchpoint.mutated$distance * -1
  branchpoint.mutated.ref <- branchpoint.mutated[branchpoint.mutated$allele_status ==
                                                  "REF",]

  #make plots for reference sequence

  #probability score by sequence position
  plot.prob.ref <- ggplot(branchpoint.mutated.ref, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                           col = 'nucleotide',alpha = 'U2_binding_energy')) +
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                       breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                          "0.60","0.70","0.80","0.90","1.00")) +
    scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                      limits = levels(branchpoint_mutated_ref$nucleotide)) +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = levels(branchpoint_mutated_ref$nucleotide)) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  #U2 binding energy by sequence position
  plot.U2.ref <- ggplot(branchpoint.mutated.ref[branchpoint.mutated.ref$branchpoint_prob >= probabilityCutoff,],
              aes_string(x = 'distance', y = 'U2_binding_energy')) +
    geom_bar(stat = "identity",width = 1) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    theme(legend.position = "none") +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

  #Sequence identity
  plot.seq.ref <- ggplot(branchpoint.mutated.ref, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
    geom_text(size = 4, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = levels(branchpoint.mutated.ref$nucleotide)) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  if (plotMutated == TRUE) {
    branchpoint.mutated.alt <- branchpoint.mutated[branchpoint.mutated$allele_status == "ALT",]

    #probability score by sequence position
    plot.prob.alt <- ggplot(branchpoint.mutated.alt, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                             col = 'nucleotide',alpha = 'U2_binding_energy')) +
      geom_bar(stat = "identity") +
      scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                         breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                            "0.60","0.70","0.80","0.90","1.00")) +
      scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                        limits = levels(branchpoint.mutated.ref$nucleotide)) +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = levels(branchpoint.mutated.ref$nucleotide)) +
      theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17)) + ggtitle("Alternative")

    plot.prob.ref <- plot.prob.ref + ggtitle("Reference")

    #U2 binding energy by sequence position
    plot.U2.alt <- ggplot(branchpoint.mutated.alt[branchpoint.mutated.alt$branchpoint_prob >= probabilityCutoff,],
                aes_string(x = 'distance', y = 'U2_binding_energy')) +
      geom_bar(stat = "identity",width = 1) +
      scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                         breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
      theme(legend.position = "none") +
      scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                         breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

    #Sequence identity
    plot.seq.alt <- ggplot(branchpoint.mutated.alt, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = levels(branchpoint.mutated.ref$nucleotide)) +
      theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17))


    ######Get reference vs. alt sequence for 0-50 window#####

    seq <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
                                     grepl("REF",queryAttributes$id))[1],]$seq
    pos <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
                                     grepl("REF",queryAttributes$id))[1],]$to_3prime
    seq.ref <- unlist(stringr::str_split(stringr::str_sub(seq,251 + (pos - 50),250 + pos),""))

    seq <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
                                     grepl("ALT",queryAttributes$id))[1],]$seq
    pos <- queryAttributes[which(grepl(queryName, queryAttributes$id) &
                                     grepl("ALT",queryAttributes$id))[1],]$to_3prime
    seq.alt <- unlist(stringr::str_split(stringr::str_sub(seq,251 + (pos - 50),250 + pos),""))

    wholeSeq <- data.frame(nt = c(seq_1, seq_2),
                           pos = rep(seq(-50,-1,1),2),
                           set = rep(c("REF","ALT"), each =50))
    mutPosition <- wholeSeq$pos[which(wholeSeq$nt[wholeSeq$set == "REF"] !=
                                    wholeSeq$nt[wholeSeq$set == "ALT"])]

    plot.seqComparison <- ggplot(wholeSeq, aes_string(x = 'pos', y = 'set', col = 'nt',label = 'nt')) +
      geom_rect(xmin = mutPosition - 0.5,xmax = mutPosition + 0.5, ymin = 0.5, ymax = 0.5 + 2,
                col = NA, fill = "grey90") +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols) +
      theme(legend.position = "none",
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-50,-1)) +
      scale_y_discrete(labels = c("Alternative", "Reference")) +
      geom_rect(xmin = -44.5,xmax = -17.5, ymin = 0.5, ymax = 0.5 + 2, fill = NA, col ="black")

  }

  #plot gene/transcript structure
  if (plotStructure == TRUE) {

    x <- match(exons$gene_id, exons$gene_id[grep(branchpoint.mutated$exon_3prime[1], exons$exon_id)][1])
    transcripts <- exons[which(!is.na(x)),]

    transcripts.gene <- transcripts
    transcripts.gene$transcript_id <- transcripts.gene$gene_id
    transcripts <- rbind(transcripts,transcripts.gene)
    transcripts$length <- transcripts$end - transcripts$start

    if (transcripts$strand[1] == "+") {
      bpExons = which(transcripts$start ==
                         transcripts$start[which(transcripts$exon_id ==
                                                   branchpoint.mutated$exon_3prime[1])][1])

    }else{
      bpExons <- which(transcripts$end ==
                         transcripts$end[which(transcripts$exon_id ==
                                                 branchpoint.mutated$exon_3prime[1])][1])
    }

    transcripts$SNP_exon <- 0
    transcripts$SNP_exon[bp_exons] <- 1
    transcripts$SNP_exon[transcripts$transcript_id == transcripts$gene_id[1]] <- 2

    #make smaller data.frame for plotting connecting lines
    transcriptsLines <- transcripts[!duplicated(transcripts$transcript_id),]

    max.start <- vector()
    min.start <- vector()
    max.end <- vector()
    min.end <- vector()
    for (t in seq_along(transcriptsLines$chromosome)) {
      ind <- grep(transcriptsLines$transcript_id[t], transcripts$transcript_id)
      max.start[t] <- max(transcripts$start[ind])
      min.end[t] <- min(transcripts$end[ind])
      max.end[t] <- max(transcripts$end[ind])
      min.start[t] <- min(transcripts$start[ind])

    }

    #only plot transcripts overlapping the query
    keep <- which((max.end >= max(transcripts$end[transcripts$SNP_exon == 1])) &
                   (min.start <= min(transcripts$start[transcripts$SNP_exon == 1])))

    transcriptsLines$max <- max.start
    transcriptsLines$min <- min.end
    transcriptsLines <- transcriptsLines[keep,]
    transcriptsBP <- unique(transcriptsLines$transcript_id)
    transcripts <- transcripts[transcripts$transcript_id %in% transcriptsBP,]
    transcripts$transcript_id_num <- as.numeric(as.factor(transcripts$transcript_id)) * -1
    transcriptsLines$transcript_id_num <- as.numeric(as.factor(transcriptsLines$transcript_id)) * -1

    if (transcripts$strand[1] == "+") {
      transcriptsLines$transcript_id <- paste0(transcriptsLines$transcript_id, " (+)")
    }else{
      transcriptsLines$transcript_id <- paste0(transcriptsLines$transcript_id, " (-)")
    }

    plot.structure <- ggplot(transcripts, aes_string(xmin = 'start', xmax = 'start + length',
                                              ymin = 'transcript_id_num - 0.4',
                                              ymax = 'transcript_id_num + 0.4', fill = 'factor(SNP_exon)')) +
      geom_segment(data = transcriptsLines, aes_string(x = 'min',xend = 'max', y = 'transcript_id_num',
                                                 yend = 'transcript_id_num')) +
      geom_rect() +
      scale_fill_manual(values = c("black","blue","grey60")) +
      scale_y_continuous(breaks = seq(-1, min(transcriptsLines$transcript_id_num),-1),
                         labels = transcriptsLines$transcript_id[
                           match(seq(-1, min(transcriptsLines$transcript_id_num),-1),
                                 transcriptsLines$transcript_id_num)]) +
      theme(panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

  }
  theme_set(theme_gray())
  if (plotStructure == TRUE & plotMutated == TRUE) {
    ggdraw() +
      draw_plot(plot.structure,0,0.775,1,0.225) + draw_plot(plot.seq.comparison,0,0.675,1,0.1) +
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
      draw_plot(plot.seq.comparison,0,0.9,1,0.1) +
      draw_plot(plot.prob.ref,0,.4,0.5,.5) + draw_plot(plot.seq.ref,0,.3,0.5,.1) + draw_plot(plot.U2.ref,0,0,0.5,.3) +
      draw_plot(plot.prob.alt,0.5,.4,0.5,.5) + draw_plot(plot.seq.alt,0.5,.3,0.5,.1) + draw_plot(plot.U2.alt,0.5,0,0.5,.3)
  }

}
