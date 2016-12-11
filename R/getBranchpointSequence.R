#' Get branchpoint sequence features
#' Gets intronic sequence covering the branchpoint window and extracts predictive features
#' @param query branchpointer query data.frame
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and strand at position 4.
#' @param rm_chr remove "chr" before chromosome names before writing bed file.
#' Required if genome sequence names do not contain "chr"
#' @param query_type type of branchpointer query. "SNP" or "region".
#' @param genome .fa genome file location
#' @param bedtools_location bedtools binary location (which bedtools)
#' @param unique_id unique string identifier for intermediate .bed and .fa files.
#' @param working_directory directory where intermediate .bed and .fa are located
#' @param use_parallel use parallelisation to speed up code?
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @param useBSgenome Use a BSgenome object for sequence retrieval?  (default = \code{1})
#' Overridden if a genome .fa and a bedtools_location are sepcified.
#' @param BSgenome BSgenome object
#' @return data.frame with all features required to predict branchpoint probability scores
#' @export
#' @import data.table
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom Biostrings complement
#' @examples
#' query_attributes <- getBranchpointSequence(query,
#' query_type = "SNP"
#' genome = "GRCh38.p5.genome.fa",
#' bedtools_location = "/Applications/apps/bedtools2/bin/bedtools")
#'
#' genome <- BSgenome.Hsapiens.UCSC.hg38
#' query_attributes <- getBranchpointSequence(query,
#' query_type = "region"
#' useBSgenome = TRUE,
#' BSgenome = genome)

#' @author Beth Signal

getBranchpointSequence <- function(query, unique_id = "test",
                                   query_type,
                                   working_directory = ".",
                                   genome = NA,
                                   rm_chr = FALSE,
                                   bedtools_location,
                                   use_parallel = FALSE,
                                   cores = 1,
                                   useBSgenome = FALSE,
                                   BSgenome = NULL) {


  if(missing(query_type) | !(query_type %in% c("SNP", "region"))){

    stop("please specify query_type as \"region\" or \"SNP\"")

  }

  if(is.na(genome) & useBSgenome == FALSE){

        stop("please specify a genome .fa file for sequence extraction or set useBSgenome to TRUE and specify a BSgenome object")

  }

  if(is.na(genome) & useBSgenome == TRUE & is.null(BSgenome)){

    stop("please specify a BSgenome object")

  }

  if(is.na(genome) & useBSgenome == FALSE & !missing(bedtools_location)){

    stop("please specify a genome .fa file for sequence extraction")

  }

  if(missing(bedtools_location) & !is.na(genome) & (useBSgenome == FALSE | (useBSgenome == TRUE & is.null(BSgenome)))){

    stop("please specify the bedtools binary location")

  }

  if(!is.na(genome) & !missing(bedtools_location) & useBSgenome == TRUE & !is.null(BSgenome)){

    stop("Both a .fa genome and BSgenome have been specified.\n
         Using the .fa genome...\n
         to use BSgenome, don't specify a .fa file")

  }

  if(use_parallel){

    max_cores <- parallel::detectCores()

    if(max_cores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", max_cores,")"))
      message(paste0("using all available cores"))
      cores <- max_cores

    }

  }

  #make bed format file
  if (query_type == "SNP") {
    bed <- query[,c(2,3,3,1,4)]

  }else if (query_type == "region") {
    bed <- query[,c(2,4,4,1,5)]
    bed[bed[,5] == "-",] <- query[bed[,5] == "-",c(2,3,3,1,5)]
  }

  bed[,2] <- (bed[,3] - 1)
  colnames(bed)[2:3] <- c("start","end")
  bed$score <- 0

  #extend bed file to cover +/- 250 nt from each query point
  if (length(which(bed$strand == "+")) > 0) {
    bed[which(bed$strand == "+"),]$start <-
      bed[which(bed$strand == "+"),]$end - 251 - (44 - query$to_3prime[which(bed$strand ==
                                                                               "+")])
    bed[which(bed$strand == "+"),]$end <-
      bed[which(bed$strand == "+"),]$start + 501 + 27
  }

  if (length(which(bed$strand == "-")) > 0) {
    bed[which(bed$strand == "-"),]$end <-
      bed[which(bed$strand == "-"),]$start + 251 + (44 - query$to_3prime[which(bed$strand ==
                                                                                 "-")])
    bed[which(bed$strand == "-"),]$start <-
      bed[which(bed$strand == "-"),]$end - (501 + 27)
  }

  bed <- bed[,c(1,2,3,4,6,5)]

  if (rm_chr == TRUE) {
    bed[,1] <- gsub("chr","", bed[,1])
  }

  if(useBSgenome){
    bed_seq <- getSeq(BSgenome, bed$chromosome, start=bed$start+1,
                   end=bed$end, strand=bed$strand)
    s <- as.character(bed_seq)
    ids <- rownames(bed)
  }else{
  #convert to .fasta using bedtools
    utils::write.table(
      bed, sep = "\t", file = paste0(working_directory,"/mutation_",unique_id,".bed"),
      row.names = FALSE,col.names = FALSE,quote = FALSE
    )
    cmd <- paste0(
      bedtools_location," getfasta -fi ", genome,
      " -bed ",working_directory,"/mutation_",unique_id,".bed -fo ",
      working_directory,"/mutation_",unique_id,".fa -name -s"
    )
    system(cmd)
    fasta <-
      data.table::fread(paste0(working_directory,"/mutation_",unique_id,".fa"),
                        header = FALSE, stringsAsFactors = FALSE)
    fasta <- as.data.frame(fasta)
    system(paste0("rm -f ",working_directory,"/mutation_",unique_id,"*"))

    s <- fasta[seq(2,dim(fasta)[1],by = 2),1]
    ids <- gsub(">","",fasta[seq(1,dim(fasta)[1],by = 2),1])
    m <- match(query$id,ids)
    query <- query[which(!is.na(m)),]
  }

  ##mutate at SNP location
  if (query_type == "SNP") {
    #location of SNP
    loc <- 44 - query$to_3prime

    ref_nt <- as.character(query$ref_allele)
    alt_nt <- as.character(query$alt_allele)

    #change to compliment if on negative strand
    ref_nt[query$strand == "-"] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(ref_nt[query$strand == "-"])))
    alt_nt[query$strand == "-"] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(alt_nt[query$strand == "-"])))

    #check ref allele
    ref_allele_correct <- substr(s, 251 + (loc),251 + (loc)) == ref_nt

    if (any(!ref_allele_correct)) {
      rm <- which(ref_allele_correct == FALSE)

      if (all(ref_allele_correct[query$strand == "-"] == FALSE) &
          all(ref_allele_correct[query$strand == "+"])) {
        message("reference alleles are incorrect for all negative strand introns")
        message("please input alleles as positive strand sequences")
      }else{
        message("reference alleles do not match sequence for:")
        message(paste(query$id[rm], collapse = "\n"))
      }
      message("removing from analysis")
      query <- query[-rm,]
      s <- s[-rm]
      ref_nt <- ref_nt[-rm]
      alt_nt <- alt_nt[-rm]
      bed <- bed[-rm,]
      loc <- loc[-rm]
    }

    #create mutated sequence
    s_mut <-
      paste0(substr(s, 1,250 + (loc)), alt_nt, substr(s, 252 + (loc),nchar(s)))

    seqs_mut <- vector()
    for (i in 18:44) {
      seqs_mut <- append(seqs_mut, substr(s_mut, (i - 17),(i - 17) + 500))
    }


  }

  #create 501nt sequences with each query point centered at 251
  seqs <- vector()
  #from 44 to 18 (dist.2)
  for (i in 18:44) {
    seqs <- append(seqs, substr(s, (i - 17),(i - 17) + 500))
  }

  if (query_type == "SNP") {
    df <-
      data.frame(id = c(
        paste0(bed$id,"_",rep(44:18,each = length(s)), "_REF"),
        paste0(bed$id,"_",rep(44:18,each = length(s)), "_ALT")
      ),
      seq = c(seqs, seqs_mut) , stringsAsFactors = FALSE)
    reps <- 2

  }else{
    df <-
      data.frame(id = paste0(bed$id,"_",rep(44:18,each = length(s)), "_REF"),
                 seq = seqs, stringsAsFactors = FALSE)
    reps <- 1
  }


  df$to_3prime_vec <- rep(rep(44:18,each = length(s)),reps)
  df$to_5prime_vec <-
    rep(query$to_3prime + query$to_5prime, reps) - df$to_3prime_vec

  if(query_type == "region"){
    df$end <- rep(rep(query$chrom_end, length(18:44)),reps)
  }else{
    df$end <- rep(rep(query$chrom_start, length(18:44)),reps)
  }
  start <- rep(rep(query$chrom_start, length(18:44)),reps)

  df$dist.2 <- rep(rep(query$to_3prime, length(18:44)),reps)
  df$strand <- rep(rep(query$strand, length(18:44)),reps)

  df$end[df$strand == "-"] <- start[df$strand == "-"]

  df$chromosome <- rep(rep(query$chromosome, length(18:44)),reps)
  df$exon_3prime <- rep(rep(query$exon_3prime, length(18:44)),reps)
  df$exon_5prime <- rep(rep(query$exon_5prime, length(18:44)),reps)

  df$end[df$strand == "+"] <- df$end[df$strand == "+"] +
    df$dist.2[df$strand == "+"] - df$to_3prime_vec[df$strand == "+"]
  df$end[df$strand == "-"] <- df$end[df$strand == "-"]  -
    df$dist.2[df$strand == "-"] + df$to_3prime_vec[df$strand == "-"]

  query_attributes <- data.frame(
    exon_id = df$id,
    dist.1 = df$to_5prime_vec,
    dist.2 = df$to_3prime_vec, stringsAsFactors = FALSE
  )


  #get sequence identity at position -5 to +5 relative to testing point
  seq_pos0 <-
    factor(substr(df$seq,251,251), levels = c("A","C","G","T"))
  seq_pos1 <-
    factor(substr(df$seq,252,252), levels = c("A","C","G","T"))
  seq_pos2 <-
    factor(substr(df$seq,253,253), levels = c("A","C","G","T"))
  seq_pos3 <-
    factor(substr(df$seq,254,254), levels = c("A","C","G","T"))
  seq_pos4 <-
    factor(substr(df$seq,255,255), levels = c("A","C","G","T"))
  seq_pos5 <-
    factor(substr(df$seq,256,256), levels = c("A","C","G","T"))
  seq_neg1 <-
    factor(substr(df$seq,250,250), levels = c("A","C","G","T"))
  seq_neg2 <-
    factor(substr(df$seq,249,249), levels = c("A","C","G","T"))
  seq_neg3 <-
    factor(substr(df$seq,248,248), levels = c("A","C","G","T"))
  seq_neg4 <-
    factor(substr(df$seq,247,247), levels = c("A","C","G","T"))
  seq_neg5 <-
    factor(substr(df$seq,246,246), levels = c("A","C","G","T"))

  #find canonical AG splice dinucleotides
  f <- gregexpr("AG",substr(df$seq, 252,501),perl = TRUE)

  if (use_parallel) {
    cluster <- parallel::makeCluster(cores)

    canon_hits <- parallel::parLapply(cluster,f, getCanonical3SS)
    canon_df <- matrix(unlist(canon_hits), ncol = 5, byrow = TRUE)
    canon_df <- as.data.frame(canon_df, stringsAsFactors=FALSE)
    colnames(canon_df) <-
      c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
    query_attributes <- cbind(query_attributes, seq = df$seq)

    pyra_df <-
      parallel::parApply(cluster,query_attributes, 1, getPPT)

    parallel::stopCluster(cluster)
  }else{
    canon_hits <- lapply(f, getCanonical3SS)
    canon_df <- matrix(unlist(canon_hits), ncol = 5, byrow = TRUE)
    canon_df <- as.data.frame(canon_df, stringsAsFactors=FALSE)
    colnames(canon_df) <-
      c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
    query_attributes <- cbind(query_attributes, seq = df$seq)

    pyra_df <- apply(query_attributes, 1, getPPT)
  }

  pyra_df <- as.data.frame(t(pyra_df),stringsAsFactors=FALSE)
  colnames(pyra_df) <- c("ppt_start","ppt_run_length")
  query_attributes$seq <- NULL

  query_attributes <-
    cbind(
      query_attributes, pyra_df, canon_df, seq_neg5,seq_neg4,seq_neg3,seq_neg2,
      seq_neg1,seq_pos0,seq_pos1,seq_pos2,seq_pos3,seq_pos4,seq_pos5,
      row.names = NULL, stringsAsFactors=FALSE
    )

  df <- cbind(df[,c(1,8,5,7,2,9,10)], query_attributes[,c(-c(1))])

  colnames(df)[c(8,9)] <- c("to_5prime","to_3prime")

  return(df)
}
