#' Plots branchpointer predictions
#'
#' Plots branchpointer predictions
#' @param query_name query id used to identify the SNP or region
#' @param query_attributes attribute data.frame
#' @param branchpoint_predictions prediction of branchpoint probabilities for each site
#' @param probability_cutoff probability score cutoff value for displaying U2 binding energy
#' @param plot_structure plot structures for gene and 3' exon containing and skipping isoforms
#' @param plot_mutated plot alternative sequence predicitons alongside reference sequence predictions
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return ggplot2 plot with branchpoint features in the specified intronic region
#' @export
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @examples
#' exons <- readExonAnnotation("gencode.v24.exons.txt")
#' query <- readQueryFile("SNP_example.txt",query_type = "SNP")
#' query <- getQueryLoc(query,query_type="SNP",exons = exons, filter=FALSE)
#' query_attributes <- getBranchpointSequence(query,query_type = "SNP",
#' genome = "GRCh38.p5.genome.fa",
#' bedtools_location="/Applications/bedtools2/bin/bedtools")
#' branchpoint_predictions <- predictBranchpoints(query_attributes)
#' plotBranchpointWindow(query$id[1],  branchpoint_predictions,
#' query_attributes, plot_mutated = TRUE, exons = exons)
#' @author Beth Signal

plotBranchpointWindow <- function(query_name,
                                 branchpoint_predictions,
                                 query_attributes,
                                 probability_cutoff = 0.5,
                                 plot_mutated = FALSE,
                                 plot_structure = TRUE,
                                 exons) {
  #find the query
  ind <- which(!is.na(match(
    branchpoint_predictions$id, query_name
  )))
  branchpoint_mutated <- branchpoint_predictions[ind,]

  branchpoint_mutated$distance <- branchpoint_mutated$distance * -1
  branchpoint_mutated_ref <- branchpoint_mutated[branchpoint_mutated$allele_status ==
                                                  "REF",]

  #make plots for reference sequence

  #probability score by sequence position
  plot_prob_ref <- ggplot(branchpoint_mutated_ref, aes(x = distance, y = branchpoint_prob, fill = nucleotide,
                                           col = nucleotide,alpha = U2_binding_energy)) +
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
  plot_U2_ref <- ggplot(branchpoint_mutated_ref[branchpoint_mutated_ref$branchpoint_prob >= probability_cutoff,],
              aes(x = distance, y = U2_binding_energy)) +
    geom_bar(stat = "identity",width = 1) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    theme(legend.position = "none") +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

  #Sequence identity
  plot_seq_ref <- ggplot(branchpoint_mutated_ref, aes(x = distance, y = 1, col = nucleotide,label = nucleotide)) +
    geom_text(size = 4, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = levels(branchpoint_mutated_ref$nucleotide)) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  if (plot_mutated == TRUE) {
    branchpoint_mutated_alt <- branchpoint_mutated[branchpoint_mutated$allele_status == "ALT",]

    #probability score by sequence position
    plot_prob_alt <- ggplot(branchpoint_mutated_alt, aes(x = distance, y = branchpoint_prob, fill = nucleotide,
                                             col = nucleotide,alpha = U2_binding_energy)) +
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
      scale_x_continuous(limits = c(-45,-17)) + ggtitle("Alternative")

    plot_prob_ref <- plot_prob_ref + ggtitle("Reference")

    #U2 binding energy by sequence position
    plot_U2_alt <- ggplot(branchpoint_mutated_alt[branchpoint_mutated_alt$branchpoint_prob >= probability_cutoff,],
                aes(x = distance, y = U2_binding_energy)) +
      geom_bar(stat = "identity",width = 1) +
      scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                         breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
      theme(legend.position = "none") +
      scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                         breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

    #Sequence identity
    plot_seq_alt <- ggplot(branchpoint_mutated_alt, aes(x = distance, y = 1, col = nucleotide,label = nucleotide)) +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = levels(branchpoint_mutated_ref$nucleotide)) +
      theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17))


    ######Get reference vs. alt sequence for 0-50 window#####

    seq <- query_attributes[which(grepl(query_name, query_attributes$id) &
                                     grepl("REF",query_attributes$id))[1],]$seq
    pos <- query_attributes[which(grepl(query_name, query_attributes$id) &
                                     grepl("REF",query_attributes$id))[1],]$to_3prime
    seq_1 <- unlist(str_split(str_sub(seq,251 + (pos - 50),250 + pos),""))

    seq <- query_attributes[which(grepl(query_name, query_attributes$id) &
                                     grepl("ALT",query_attributes$id))[1],]$seq
    pos <- query_attributes[which(grepl(query_name, query_attributes$id) &
                                     grepl("ALT",query_attributes$id))[1],]$to_3prime
    seq_2 <- unlist(str_split(str_sub(seq,251 + (pos - 50),250 + pos),""))

    whole_seq <- data.frame(nt = c(seq_1, seq_2),
                           pos = rep(seq(-50,-1,1),2),
                           set = rep(c("REF","ALT"), each =50))
    mut_pos <- whole_seq$pos[which(whole_seq$nt[whole_seq$set == "REF"] !=
                                    whole_seq$nt[whole_seq$set == "ALT"])]

    plot_seq_comparison <- ggplot(whole_seq, aes(x = pos, y = set, col = nt,label = nt)) +
      geom_rect(xmin = mut_pos - 0.5,xmax = mut_pos + 0.5, ymin = 0.5, ymax = 0.5 + 2,
                col = NA, fill = "grey90") +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values =mercer_nt_cols) +
      theme(legend.position = "none",
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-50,-1)) +
      scale_y_discrete(labels = c("Alternative", "Reference")) +
      geom_rect(xmin = -44.5,xmax = -17.5, ymin = 0.5, ymax = 0.5 + 2, fill = NA, col ="black")

  }

  #plot gene/transcript structure
  if (plot_structure == TRUE) {

    x <- match(exons$gene_id, exons$gene_id[grep(branchpoint_mutated$exon_3prime[1], exons$exon_id)][1])
    transcripts <- exons[which(!is.na(x)),]

    transcripts_gene <- transcripts
    transcripts_gene$transcript_id <- transcripts_gene$gene_id
    transcripts <- rbind(transcripts,transcripts_gene)
    transcripts$length <- transcripts$end - transcripts$start

    if (transcripts$strand[1] == "+") {
      bp_exons = which(transcripts$start ==
                         transcripts$start[which(transcripts$exon_id ==
                                                   branchpoint_mutated$exon_3prime[1])][1])

    }else{
      bp_exons <- which(transcripts$end ==
                         transcripts$end[which(transcripts$exon_id ==
                                                 branchpoint_mutated$exon_3prime[1])][1])
    }

    transcripts$SNP_exon <- 0
    transcripts$SNP_exon[bp_exons] <- 1
    transcripts$SNP_exon[transcripts$transcript_id == transcripts$gene_id[1]] <- 2

    #make smaller data.frame for plotting connecting lines
    transcripts_lines <- transcripts[!duplicated(transcripts$transcript_id),]

    maxs <- vector()
    mins <- vector()
    maxe <- vector()
    mine <- vector()
    for (t in seq_along(transcripts_lines$chromosome)) {
      ind <- grep(transcripts_lines$transcript_id[t], transcripts$transcript_id)
      maxs[t] <- max(transcripts$start[ind])
      mine[t] <- min(transcripts$end[ind])
      maxe[t] <- max(transcripts$end[ind])
      mins[t] <- min(transcripts$start[ind])

    }

    #only plot transcripts overlapping the query
    keep <- which((maxe >= max(transcripts$end[transcripts$SNP_exon == 1])) &
                   (mins <= min(transcripts$start[transcripts$SNP_exon == 1])))

    transcripts_lines$max <- maxs
    transcripts_lines$min <- mine
    transcripts_lines <- transcripts_lines[keep,]
    BP_transcripts <- unique(transcripts_lines$transcript_id)
    transcripts <- transcripts[transcripts$transcript_id %in% BP_transcripts,]
    transcripts$transcript_id_num <- as.numeric(as.factor(transcripts$transcript_id)) * -1
    transcripts_lines$transcript_id_num <- as.numeric(as.factor(transcripts_lines$transcript_id)) * -1

    if (transcripts$strand[1] == "+") {
      transcripts_lines$transcript_id <- paste0(transcripts_lines$transcript_id, " (+)")
    }else{
      transcripts_lines$transcript_id <- paste0(transcripts_lines$transcript_id, " (-)")
    }

    structure_plot <- ggplot(transcripts, aes(xmin = start, xmax = start + length,
                                             ymin = transcript_id_num - 0.4,
                                             ymax = transcript_id_num + 0.4, fill = factor(SNP_exon))) +
      geom_segment(data = transcripts_lines, aes(x = min,xend = max, y = transcript_id_num,
                                                 yend = transcript_id_num)) +
      geom_rect() +
      scale_fill_manual(values = c("black","blue","grey60")) +
      scale_y_continuous(breaks = seq(-1, min(transcripts_lines$transcript_id_num),-1),
                         labels = transcripts_lines$transcript_id[
                           match(seq(-1, min(transcripts_lines$transcript_id_num),-1),
                                 transcripts_lines$transcript_id_num)]) +
      theme(panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

  }
  theme_set(theme_gray())
  if (plot_structure == TRUE & plot_mutated == TRUE) {
    ggdraw() +
      draw_plot(structure_plot,0,0.775,1,0.225) + draw_plot(plot_seq_comparison,0,0.675,1,0.1) +
      draw_plot(plot_prob_ref,0,.275,0.5,.40) + draw_plot(plot_seq_ref,0,.2,0.5,.075) + draw_plot(plot_U2_ref,0,0,0.5,.2) +
      draw_plot(plot_prob_alt,0.5,.275,0.5,.40) + draw_plot(plot_seq_alt,0.5,.2,0.5,.075) + draw_plot(plot_U2_alt,0.5,0,0.5,.2)

  }else if (plot_structure == TRUE & plot_mutated == FALSE) {
    ggdraw() +
      draw_plot(structure_plot,0,0.775,1,0.225) +
      draw_plot(plot_prob_ref,0,.325,1,.45) + draw_plot(plot_seq_ref,0,.25,1,.075) + draw_plot(plot_U2_ref,0,0,1,.25)
  }else if (plot_structure == FALSE & plot_mutated == FALSE) {
    ggdraw() +
      draw_plot(plot_prob_ref,0,.325,1,.45) + draw_plot(plot_seq_ref,0,.25,1,.075) + draw_plot(plot_U2_ref,0,0,1,.25)
  }else if (plot_structure == FALSE & plot_mutated == TRUE) {
    ggdraw() +
      draw_plot(plot_seq_comparison,0,0.9,1,0.1) +
      draw_plot(plot_prob_ref,0,.4,0.5,.5) + draw_plot(plot_seq_ref,0,.3,0.5,.1) + draw_plot(plot_U2_ref,0,0,0.5,.3) +
      draw_plot(plot_prob_alt,0.5,.4,0.5,.5) + draw_plot(plot_seq_alt,0.5,.3,0.5,.1) + draw_plot(plot_U2_alt,0.5,0,0.5,.3)
  }

}
