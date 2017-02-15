#' Get branchpoint sequence features
#' Gets intronic sequence covering the branchpoint window and extracts predictive features
#' @param query branchpointer query data.frame
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and strand at position 4.
#' @param rmChr remove "chr" before chromosome names before writing bed file.
#' Required if genome sequence names do not contain "chr"
#' @param queryType type of branchpointer query. "SNP" or "region".
#' @param genome .fa genome file location
#' @param bedtoolsLocation bedtools binary location (which bedtools)
#' @param uniqueId unique string identifier for intermediate .bed and .fa files.
#' @param workingDirectory directory where intermediate .bed and .fa are located
#' @param useParallel use parallelisation to speed up code?
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @param useBSgenome Use a BSgenome object for sequence retrieval?  (default = \code{1})
#' Overridden if a genome .fa and a bedtools_location are sepcified.
#' @param BSgenome BSgenome object
#' @return data.frame with all features required to predict branchpoint probability scores
#' @export
#' @import data.table
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom Biostrings complement
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @examples
#' small_exons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(small_exons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' query_snp <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(query_snp,queryType = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' query_attributes <- getBranchpointSequence(query,
#' queryType = "SNP",
#' useBSgenome = TRUE,
#' BSgenome = genome)
#'
#' query <- makeRegions("ENSE00003541068.1", "exon_id", exons)
#' query_attributes <- getBranchpointSequence(query,
#' queryType = "region",
#' useBSgenome = TRUE,
#' BSgenome = genome)

#' @author Beth Signal

getBranchpointSequence <- function(query, uniqueId = "test",
                                   queryType,
                                   workingDirectory = ".",
                                   genome = NA,
                                   rmChr = FALSE,
                                   bedtoolsLocation,
                                   useParallel = FALSE,
                                   cores = 1,
                                   useBSgenome = FALSE,
                                   BSgenome = NULL) {


  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){

    stop("please specify queryType as \"region\" or \"SNP\"")

  }

  if(is.na(genome) & useBSgenome == FALSE){

        stop("please specify a genome .fa file for sequence extraction or set useBSgenome to TRUE and specify a BSgenome object")

  }

  if(is.na(genome) & useBSgenome == TRUE & is.null(BSgenome)){

    stop("please specify a BSgenome object")

  }

  if(is.na(genome) & useBSgenome == FALSE & !missing(bedtoolsLocation)){

    stop("please specify a genome .fa file for sequence extraction")

  }

  if(missing(bedtoolsLocation) & !is.na(genome) & (useBSgenome == FALSE | (useBSgenome == TRUE & is.null(BSgenome)))){

    stop("please specify the bedtools binary location")

  }

  if(!is.na(genome) & !missing(bedtoolsLocation) & useBSgenome == TRUE & !is.null(BSgenome)){

    stop("Both a .fa genome and BSgenome have been specified.\n
         Using the .fa genome...\n
         to use BSgenome, don't specify a .fa file")

  }

  if(useParallel){

    maxCores <- parallel::detectCores()

    if(maxCores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", maxCores,")"))
      message(paste0("using all available cores"))
      cores <- maxCores

    }

  }

  #make bed format file
  if (queryType == "SNP") {
    bed <- query[,c(2,3,3,1,4)]

  }else if (queryType == "region") {
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
    bed.seq <- Biostrings::getSeq(BSgenome, bed$chromosome, start=bed$start+1,
                   end=bed$end, strand=bed$strand)
    s <- as.character(bed.seq)
    ids <- rownames(bed)
  }else{
  #convert to .fasta using bedtools
    utils::write.table(
      bed, sep = "\t", file = paste0(workingDirectory,"/mutation_",uniqueId,".bed"),
      row.names = FALSE,col.names = FALSE,quote = FALSE
    )
    cmd <- paste0(
      bedtoolsLocation," getfasta -fi ", genome,
      " -bed ",workingDirectory,"/mutation_",uniqueId,".bed -fo ",
      workingDirectory,"/mutation_",uniqueId,".fa -name -s"
    )
    system(cmd)
    fasta <-
      data.table::fread(paste0(workingDirectory,"/mutation_",uniqueId,".fa"),
                        header = FALSE, stringsAsFactors = FALSE)
    fasta <- as.data.frame(fasta)
    system(paste0("rm -f ",workingDirectory,"/mutation_",uniqueId,"*"))

    s <- fasta[seq(2,dim(fasta)[1],by = 2),1]
    ids <- gsub(">","",fasta[seq(1,dim(fasta)[1],by = 2),1])
    m <- match(query$id,ids)
    query <- query[which(!is.na(m)),]
  }

  ##mutate at SNP location
  if (queryType == "SNP") {
    #location of SNP
    loc <- 44 - query$to_3prime

    ref.nt <- as.character(query$ref_allele)
    alt.nt <- as.character(query$alt_allele)

    #change to compliment if on negative strand
    ref.nt[query$strand == "-"] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(ref.nt[query$strand == "-"])))
    alt.nt[query$strand == "-"] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(alt.nt[query$strand == "-"])))

    #check ref allele
    refAlleleCorrect <- substr(s, 251 + (loc),251 + (loc)) == ref_nt

    if (any(!refAlleleCorrect)) {
      rm <- which(refAlleleCorrect == FALSE)

      if (all(refAlleleCorrect[query$strand == "-"] == FALSE) &
          all(refAlleleCorrect[query$strand == "+"])) {
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
    s.mut <-
      paste0(substr(s, 1,250 + (loc)), alt.nt, substr(s, 252 + (loc),nchar(s)))

    seqs.mut <- vector()
    for (i in 18:44) {
      seqs.mut <- append(seqs.mut, substr(s.mut, (i - 17),(i - 17) + 500))
    }


  }

  #create 501nt sequences with each query point centered at 251
  seqs <- vector()
  #from 44 to 18 (dist.2)
  for (i in 18:44) {
    seqs <- append(seqs, substr(s, (i - 17),(i - 17) + 500))
  }

  if (queryType == "SNP") {
    df <-
      data.frame(id = c(
        paste0(bed$id,"_",rep(44:18,each = length(s)), "_REF"),
        paste0(bed$id,"_",rep(44:18,each = length(s)), "_ALT")
      ),
      seq = c(seqs, seqs.mut) , stringsAsFactors = FALSE)
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

  if(queryType == "region"){
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

  queryAttributes <- data.frame(
    exon_id = df$id,
    dist.1 = df$to_5prime_vec,
    dist.2 = df$to_3prime_vec, stringsAsFactors = FALSE
  )


  #get sequence identity at position -5 to +5 relative to testing point
  seq.pos.0 <-
    factor(substr(df$seq,251,251), levels = c("A","C","G","T"))
  seq.pos.1 <-
    factor(substr(df$seq,252,252), levels = c("A","C","G","T"))
  seq.pos.2 <-
    factor(substr(df$seq,253,253), levels = c("A","C","G","T"))
  seq.pos.3 <-
    factor(substr(df$seq,254,254), levels = c("A","C","G","T"))
  seq.pos.4 <-
    factor(substr(df$seq,255,255), levels = c("A","C","G","T"))
  seq.pos.5 <-
    factor(substr(df$seq,256,256), levels = c("A","C","G","T"))
  seq.neg.1 <-
    factor(substr(df$seq,250,250), levels = c("A","C","G","T"))
  seq.neg.2 <-
    factor(substr(df$seq,249,249), levels = c("A","C","G","T"))
  seq.neg.3 <-
    factor(substr(df$seq,248,248), levels = c("A","C","G","T"))
  seq.neg.4 <-
    factor(substr(df$seq,247,247), levels = c("A","C","G","T"))
  seq.neg.5 <-
    factor(substr(df$seq,246,246), levels = c("A","C","G","T"))

  #find canonical AG splice dinucleotides
  f <- gregexpr("AG",substr(df$seq, 252,501),perl = TRUE)

  if (useParallel) {
    cluster <- parallel::makeCluster(cores)

    canonHits <- parallel::parLapply(cluster,f, getCanonical3SS)
    canon <- matrix(unlist(canonHits), ncol = 5, byrow = TRUE)
    canon <- as.data.frame(canon, stringsAsFactors=FALSE)
    colnames(canon) <-
      c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
    queryAttributes <- cbind(queryAttributes, seq = df$seq)

    pyra <-
      parallel::parApply(cluster,queryAttributes, 1, getPPT)

    parallel::stopCluster(cluster)
  }else{
    canonHits <- lapply(f, getCanonical3SS)
    canon <- matrix(unlist(canonHits), ncol = 5, byrow = TRUE)
    canon <- as.data.frame(canon, stringsAsFactors=FALSE)
    colnames(canon) <-
      c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
    queryAttributes <- cbind(queryAttributes, seq = df$seq)

    pyra <- apply(queryAttributes, 1, getPPT)
  }

  pyra <- as.data.frame(t(pyra),stringsAsFactors=FALSE)
  colnames(pyra) <- c("ppt_start","ppt_run_length")
  queryAttributes$seq <- NULL

  queryAttributes <-
    cbind(
      queryAttributes, pyra, canon, seq.neg.5,seq.neg.4,seq.neg.3,seq.neg.2,
      seq.neg.1,seq.pos.0,seq.pos.1,seq.pos.2,seq.pos.3,seq.pos.4,seq.pos.5,
      row.names = NULL, stringsAsFactors=FALSE
    )

  df <- cbind(df[,c(1,8,5,7,2,9,10)], queryAttributes[,c(-c(1))])

  colnames(df)[c(8,9)] <- c("to_5prime","to_3prime")

  return(df)
}
