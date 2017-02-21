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
#' BSgenome = genome)
#'
#' query <- makeRegions("ENSE00003541068.1", "exon_id", exons)
#' query_attributes <- getBranchpointSequence(query,
#' queryType = "region",
#' BSgenome = genome)

#' @author Beth Signal

getBranchpointSequence <- function(query, uniqueId = "test",
                                   queryType,
                                   workingDirectory = ".",
                                   genome = NA,
                                   bedtoolsLocation,
                                   BSgenome = NULL,
                                   useParallel = FALSE,
                                   cores = 1,
                                   rmChr = FALSE) {
  
  if(is.na(genome) & is.null(BSgenome)){
        stop("please specify a genome .fa file for sequence extraction or specify a BSgenome object")
  }

  if(is.na(genome) | missing(bedtoolsLocation)){
    stop("please specify a genome .fa file for sequence extraction and a bedtools binary location")
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
  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){
    stop("please specify queryType as \"region\" or \"SNP\"")
  }else if (queryType == "SNP") {
    bed <- query
    bed@ranges@start <- as.integer(bed@ranges@start - 1)
  }else if (queryType == "region") {
    bed <- query
    bed@ranges@start[as.logical(bed@strand== "+")] <- as.integer(bed@ranges@start + bed@ranges@width - 1)[as.logical(bed@strand== "+")]
    bed@ranges@start <- bed@ranges@start - as.integer(1)
    bed@ranges@width[1:length(bed@ranges@width)] <- as.integer(2)
  }

  #extend bed file to cover +/- 250 nt from each query point
  bed@ranges@width[1:length(bed@ranges@width)] <- as.integer(529)
  bed@ranges@start[which(as.logical(bed@strand == "+"))] <- as.integer(bed@ranges@start - 250 - (44 - query@elementMetadata$to_3prime))[which(as.logical(bed@strand == "+"))]
  bed@ranges@start[which(as.logical(bed@strand == "-"))] <- as.integer(bed@ranges@start - 277 + (44 - query@elementMetadata$to_3prime))[which(as.logical(bed@strand == "-"))]

  
  if(!is.na(genome) & !missing(bedtoolsLocation)){
    #convert to .fasta using bedtools
    bedTable <- data.frame(bed@seqnames, 
                           bed@ranges@start,
                           bed@ranges@start + bed@ranges@width -1,
                           bed@elementMetadata$id,
                           score=0,
                           bed@strand)    
    
    if (rmChr == TRUE) {
      bedTable[,1] <- gsub("chr","", bedTable[,1])
    }
    
    utils::write.table(
      bedTable, sep = "\t", file = paste0(workingDirectory,"/mutation_",uniqueId,".bed"),
      row.names = FALSE,col.names = FALSE,quote = FALSE)
    cmd <- paste0(
      bedtoolsLocation," getfasta -fi ", genome,
      " -bed ",workingDirectory,"/mutation_",uniqueId,".bed -fo ",
      workingDirectory,"/mutation_",uniqueId,".fa -name -s")
    system(cmd)
    fasta <-
      data.table::fread(paste0(workingDirectory,"/mutation_",uniqueId,".fa"),
                        header = FALSE, stringsAsFactors = FALSE)
    fasta <- as.data.frame(fasta)
    system(paste0("rm -f ",workingDirectory,"/mutation_",uniqueId,"*"))
    
    s <- fasta[seq(2,dim(fasta)[1],by = 2),1]
    query@elementMetadata$seq <- s
    
    }else{
      # need to +1 for BSgenomes sequence retreval 
      # ranges given are bedtools (legacy) 
      bed@ranges@start <- as.integer(bed@ranges@start +1)
      bed@ranges@width <- as.integer(bed@ranges@width -1)
      
      bed.seq <- Biostrings::getSeq(BSgenome, bed)
      query@elementMetadata$seq <- as.character(bed.seq)
    }

  ##mutate at SNP location
  if (queryType == "SNP") {
    #location of SNP
    loc <- 44 - query@elementMetadata$to_3prime
    nt.ref <- as.character(query@elementMetadata$ref_allele)
    nt.alt <- as.character(query@elementMetadata$alt_allele)

    #change to compliment if on negative strand
    nt.ref[which(as.logical(query@strand == "-"))] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(nt.ref[which(as.logical(query@strand == "-"))])))
    nt.alt[which(query$strand == "-")] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(nt.alt[which(as.logical(query@strand == "-"))])))

    #check ref allele
    refAlleleCorrect <- substr(query@elementMetadata$seq, 251 + (loc),251 + (loc)) == nt.ref

    if (any(!refAlleleCorrect)) {
      rm <- which(refAlleleCorrect == FALSE)
      if (all(refAlleleCorrect[which(as.logical(query@strand == "-"))] == FALSE) &
          all(refAlleleCorrect[which(as.logical(query@strand == "+"))])) {
        message("reference alleles are incorrect for all negative strand introns")
        message("please input alleles as positive strand sequences")
      }else{
        message("reference alleles do not match sequence for:")
        message(paste(query@elementMetadata$id[rm], collapse = "\n"))
      }
      message("removing from analysis")
      query <- query[-rm]
    }
  }

  #create 501nt sequences with each query point centered at 251
  seqs <- vector()
  #from 44 to 18 (dist.2)
  for (i in 18:44) {
    seqs <- append(seqs, substr(query@elementMetadata$seq, (i - 17),(i - 17) + 500))
  }

  if (queryType == "SNP") {
    #create mutated sequence
    s.mut <-
      paste0(substr(query@elementMetadata$seq, 1,250 + (loc)), nt.alt, substr(s, 252 + (loc),528))
    
    seqs.mut <- vector()
    for (i in 18:44) {
      seqs.mut <- append(seqs.mut, substr(s.mut, (i - 17),(i - 17) + 500))
    }
    
    queryAllPoints <- do.call("c",as.list(rep(query, 27)))
    queryAllPoints <- do.call("c",as.list(rep(queryAllPoints, 2)))
    
    queryAllPoints@elementMetadata$status <- c(rep("REF", length(seqs)),
                                               rep("ALT", length(seqs.mut)))
    
    queryAllPoints@elementMetadata$seq <- c(seqs, seqs.mut)
    queryAllPoints@elementMetadata$to_3prime_point <- rep(rep(44:18,each = length(query@elementMetadata$id)),2)
  }else{
    queryAllPoints <- do.call("c",as.list(rep(query, 27)))
    queryAllPoints@elementMetadata$status <- rep("REF", length(seqs))
    queryAllPoints@elementMetadata$seq <- seqs
    queryAllPoints@elementMetadata$to_3prime_point <- rep(44:18,each = length(query@elementMetadata$seq))
  }
  
  queryAllPoints@elementMetadata$to_5prime_point <- (queryAllPoints@elementMetadata$to_3prime + queryAllPoints@elementMetadata$to_5prime) -
    queryAllPoints@elementMetadata$to_3prime_point
  
  testSite <- queryAllPoints@ranges@start
  posStrand <- which(as.logical(queryAllPoints@strand == "+"))
  negStrand <- which(as.logical(queryAllPoints@strand == "-"))
  testSite[which(as.logical(queryAllPoints@strand == "+"))] <- (queryAllPoints@ranges@start + queryAllPoints@ranges@width -1)[which(as.logical(queryAllPoints@strand == "+"))]
  testSite[posStrand] <- (testSite + queryAllPoints@elementMetadata$to_3prime - queryAllPoints@elementMetadata$to_3prime_point)[posStrand]
  testSite[negStrand] <- (testSite - queryAllPoints@elementMetadata$to_3prime + queryAllPoints@elementMetadata$to_3prime_point)[negStrand]
  queryAllPoints@elementMetadata$test_site <- testSite

  #get sequence identity at position -5 to +5 relative to testing point
  queryAllPoints@elementMetadata$seq_pos0 <-
    factor(substr(queryAllPoints$seq,251,251), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_pos1 <-
    factor(substr(queryAllPoints$seq,252,252), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_pos2 <-
    factor(substr(queryAllPoints$seq,253,253), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_pos3 <-
    factor(substr(queryAllPoints$seq,254,254), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_pos4 <-
    factor(substr(queryAllPoints$seq,255,255), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$ seq_pos5 <-
    factor(substr(queryAllPoints$seq,256,256), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_neg1 <-
    factor(substr(queryAllPoints$seq,250,250), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_neg2 <-
    factor(substr(queryAllPoints$seq,249,249), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_neg3 <-
    factor(substr(queryAllPoints$seq,248,248), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_neg4 <-
    factor(substr(queryAllPoints$seq,247,247), levels = c("A","C","G","T"))
  queryAllPoints@elementMetadata$seq_neg5 <-
    factor(substr(queryAllPoints$seq,246,246), levels = c("A","C","G","T"))

  #find canonical AG splice dinucleotides
  f <- gregexpr("AG",substr(queryAllPoints@elementMetadata$seq, 252,501),perl = TRUE)

  if (useParallel) {
    cluster <- parallel::makeCluster(cores)
    canonHits <- parallel::parLapply(cluster,f, getCanonical3SS)
       queryAttributes <- cbind(queryAttributes, seq = df$seq)
    pyra <-
      parallel::parApply(cluster,queryAttributes, 1, getPPT)
    parallel::stopCluster(cluster)
  }else{
    canonHits <- lapply(f, getCanonical3SS)
    pyra <- lapply(queryAllPoints, getPPT)
  }

  canon <- matrix(unlist(canonHits), ncol = 5, byrow = TRUE)
  canon <- as.data.frame(canon, stringsAsFactors=FALSE)
  colnames(canon) <-
    c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
  
  queryAllPoints@elementMetadata <- cbind(queryAllPoints@elementMetadata, canon)
  queryAllPoints@elementMetadata$ppt_start <- unlist(lapply(pyra, "[[", 1))
  queryAllPoints@elementMetadata$ppt_run_length <- unlist(lapply(pyra, "[[", 2))
  
  queryAllPoints@elementMetadata$seq <- NULL

  
  return(queryAllPoints)
}
