#' Get locations of the first five AG 3' splice site motifs
#'
#' Takes a variable length vector of distances to the AG motif,
#' sorts and returns the first five.
#' If there are less than five elements in the vector, returns the sorted vector
#' and fills the remainder of the values with 300.
#' @param ag Vector of distances to the AG splice site motif.
#' @return Locations of the first five AG dinucleotides
#' @examples
#' sequence <- paste(c("A","T","C","G")[sample(1:4, 501, replace=TRUE)] , collapse="")
#' AGSeqPositions <- gregexpr("AG",sequence, 252,501)
#' canonHits <- getCanonical3SS(AGSeqPositions[[1]])
#' @author Beth Signal

getCanonical3SS <- function(ag){
  ag <- sort(ag)
  #get first five matches
  ag <- ag[1:5]
  #if less than 5 matches, replace remainder with 300
  ag[is.na(ag)] <- 300
  
  return(ag)
}

#' Get the best polypyramidine tract
#'
#' Takes a query genomic sequence, finds all potential polypyramidine tracts (PPTs)
#' between the test site and the annotated 3'exon.
#' Returns the distance to the start of the longest PPT, and its length.
#' @param attributesLine line from a query attributes GenomicRanges
#' @return distance to the start of the longest PPT, and its length
#' @import plyr
#' @examples
#' attributesLine <- c("id", "1000","27", paste(c("A","G",rep(c("C","T"),4))[sample(1:10, 501, replace=TRUE)] , collapse=""))
#' pyra <- getPPT(attributesLine)
#' @author Beth Signal

getPPT <- function(attributesLine){
  
  dist3prime <- as.numeric(attributesLine@elementMetadata$to_3prime_point)
  
  #get sequence between the tested site and 3'exon
  seq <- substr(attributesLine@elementMetadata$seq, 251, (250 + dist3prime))
  seq <- unlist(strsplit(seq, ""))
  
  pyramidines <- which(seq == "T" | seq == "C")
  
  if(length(pyramidines) >1){
    pyraDist <- vector()
    
    #get distance to the next pyramidine
    for(p in 2:length(pyramidines)){
      pyraDist[p-1] <- pyramidines[p] - pyramidines[p-1]
    }
    
    longestRun <- 0
    bestRun <- 0
    bestStart <- 0
    percentPyra <- 0
    start <- 1
    startVec <- vector()
    runVec <- vector()
    percentPyraVec <- vector()
    
    #find best PPT
    for(p in seq_along(pyraDist)){
      
      #For the first instance of at least 2 sequential pyramidines
      if(pyraDist[p] == 1 & longestRun == 0){
        longestRun <- 1
        start <- p
        seq.p <- seq[(pyramidines[start]):(pyramidines[start+longestRun])]
        percentPyra <- length(which(seq.p=="C"|seq.p=="T"))/length(seq.p)
        
        if(longestRun > bestRun){
          bestRun <- longestRun
          bestStart <- start
          bestPercent <- percentPyra
        }
        
        
        #allowing gaps with one purine, grow the tract
      }else if(pyraDist[p] == 1 | pyraDist[p] == 2){
        longestRun <- longestRun + 1
        seq.p <- seq[(pyramidines[start]):(pyramidines[start+longestRun])]
        percentPyra <- length(which(seq.p == "C"|seq.p == "T"))/length(seq.p)
        
        if(longestRun > bestRun){
          bestRun <- longestRun
          bestStart <- start
          bestPercent <- percentPyra
        }
        
        #once too many purines are encontered, save tract
      }else{
        startVec <- append(startVec, bestStart)
        runVec <- append(runVec, longestRun)
        percentPyraVec <- append(percentPyraVec, percentPyra)
        longestRun <- 0
        start <- p + 1
        bestStart <- start
      }
    }
    
    #add the last PPT
    startVec <- append(startVec, start)
    runVec <- append(runVec, longestRun)
    percentPyraVec <- append(percentPyraVec, percentPyra)
    
    #make into data.frame
    pyra.df <- data.frame(set=c(1:length(runVec)), runVec,startVec,percentPyraVec)
    pyra.df <- plyr::arrange(pyra.df, plyr::desc(runVec))
    
    bestStart <- pyra.df$startVec[which.max(pyra.df$runVec)]
    bestRun <- pyra.df$runVec[which.max(pyra.df$runVec)]
    
    #If PPTs are long enough, shorten to get pyramidine content above 80%
    pyra.df <- pyra.df[pyra.df$runVec >=10,]
    if(dim(pyra.df)[1] != 0){
      run.extend <- vector()
      start.extend <- vector()
      polyper.extend <- vector()
      
      for(py in seq_along(length(pyra.df[,1]))){
        seq.p <- seq[(pyramidines[pyra.df$startVec[py]]):(pyramidines[pyra.df$startVec[py]+pyra.df$runVec[py]])]
        startp <- 1
        
        run <- length(seq.p)
        
        polyper.extend[py] <- pyra.df$percentPyraVec[py]
        while(polyper.extend[py] < 0.8 & run >= 10){
          purine <- which(seq.p == "A"| seq.p == "G")
          purine.neg <- (length(seq.p) +1 -purine)*-1
          purine.both <- c(purine,purine.neg)
          minPurine <- which.min(abs(purine.both))
          
          point <- purine.both[minPurine]
          
          if(point < 0){
            rm <- length(seq.p)-point*-1
            newSeq <- seq.p[1:rm]
            run <- length(newSeq)
          }else{
            newSeq <- seq.p[(point+1):length(seq.p)]
            run <- length(newSeq)
            startp <- startp+point
          }
          polyper.extend[py] <- length(which(newSeq == "C" | newSeq == "T")) / length(newSeq)
          seq.p <- newSeq
        }
        
        start.extend[py] <- (pyramidines[pyra.df$startVec[py]])+startp-1
        run.extend[py] <- run
      }
      
      bestStart <- start.extend[which.max(run.extend)]
      bestRun <- run.extend[which.max(run.extend)]
    }
    
    
    line <- c(bestStart,bestRun)
    
    #if only 1 pyramidine
  }else if(length(pyramidines) == 1){
    line <- c(pyramidines, 1)
    #if no pyramidines
  }else{
    line <- c(0,0)
  }
  return(line)
}


#' Get branchpoint sequence features
#' Gets intronic sequence covering the branchpoint window and extracts predictive features
#' @param query branchpointer query GenomicRanges
#' @param rmChr remove "chr" before chromosome names before writing bed file.
#' Required if genome sequence names do not contain "chr"
#' @param queryType type of branchpointer query. "SNP" or "region".
#' @param genome .fa genome file location
#' @param bedtoolsLocation bedtools binary location (which bedtools)
#' @param uniqueId unique string identifier for intermediate .bed and .fa files.
#' @param workingDirectory directory where intermediate .bed and .fa are located
#' @param useParallel use parallelisation to speed up code?
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @param BSgenome BSgenome object
#' @return GenomicRanges with all features required to predict branchpoint probability scores
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom Biostrings complement
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(smallExons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,queryType = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' queryAttributes <- getBranchpointSequence(query,
#' queryType = "SNP",
#' BSgenome = genome)
#'
#' query <- makeRegions("ENSE00003541068.1", "exon_id", exons)
#' queryAttributes <- getBranchpointSequence(query,
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
  
  queryAllPoints@elementMetadata$seq <- stringr::str_sub(queryAllPoints@elementMetadata$seq,251 + 
                                                           (queryAllPoints@elementMetadata$to_3prime_point - 50),250 + queryAllPoints@elementMetadata$to_3prime_point)

  return(queryAllPoints)
}

#' Predict branchpoint probability scores
#'
#' predicts branchpoint probability scores for each query site.
#' @param query branchpointer query GenomicRanges
#' @param rmChr remove "chr" before chromosome names before writing bed file.
#' Required if genome sequence names do not contain "chr"
#' @param queryType type of branchpointer query. "SNP" or "region".
#' @param genome .fa genome file location
#' @param bedtoolsLocation bedtools binary location (which bedtools)
#' @param uniqueId unique string identifier for intermediate .bed and .fa files.
#' @param workingDirectory directory where intermediate .bed and .fa are located
#' @param useParallel use parallelisation to speed up code?
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @param BSgenome BSgenome object
#' @return GenomicRanges object with branchpoint probaility scores for each site in query
#' @export
#' @import caret
#' @importFrom kernlab predict
#' @importClassesFrom kernlab ksvm
#' @import gbm
#' @import parallel
#' @importFrom stringr str_sub
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.exons.small.txt",
#' package = "branchpointer")
#' exons <- readExonAnnotation(smallExons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,query_type = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' predictions <- predictBranchpoints(query,
#' queryType = "SNP",
#' useBSgenome = TRUE,
#' BSgenome = genome)
#' @author Beth Signal

predictBranchpoints <- function(query, uniqueId = "test",
                                     queryType,
                                     workingDirectory = ".",
                                     genome = NA,
                                     bedtoolsLocation,
                                     BSgenome = NULL,
                                     useParallel = FALSE,
                                     cores = 1,
                                     rmChr = FALSE){
  if(useParallel){

    maxCores <- parallel::detectCores()

    if(maxCores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", maxCores,")"))
      message(paste0("using all available cores"))
      cores <- maxCores

    }

  }
  
  queryAttributes <- getBranchpointSequence(query, uniqueId = uniqueId,
                                            queryType = queryType,
                                            workingDirectory = workingDirectory,
                                            genome = genome,
                                            bedtoolsLocation = bedtoolsLocation,
                                            BSgenome = BSgenome,
                                            useParallel = useParallel,
                                            cores = cores,
                                            rmChr = rmChr)

  #convert to data.frame for kernlab/caret
  
  queryAttributes.forModel <- as.data.frame(queryAttributes@elementMetadata)
  
  colNamesForDataFrame <- HCN_names[-1]
  colNamesForDataFrame <- gsub("new_ID", "id", colNamesForDataFrame)
  colNamesForDataFrame <- gsub("dist.1", "to_5prime_point", colNamesForDataFrame)
  colNamesForDataFrame <- gsub("dist.2", "to_3prime_point", colNamesForDataFrame)
  
  queryAttributes.forModel <- queryAttributes.forModel[,colNamesForDataFrame]
  queryAttributes.forModel <- cbind(set="UNKNOWN", queryAttributes.forModel)
  colnames(queryAttributes.forModel) <-  HCN_names

  #remove any rows with NA values
  is_na <- apply(queryAttributes.forModel[,-c(1,2)], 2, is.na)
  is_na_v <- apply(is_na,1, any)
  rm <- which(is_na_v == TRUE)
  if(length(rm) > 0){
    queryAttributes.forModel <- queryAttributes.forModel[-rm,]
    queryAttributes <- queryAttributes[-rm,]
  }

  #convert sequences to dummy vars
  queryAttributes.dummies <- predict(dummies, newdata=queryAttributes.forModel[,-2])
  queryAttributes.forModel <- cbind(set=queryAttributes.forModel$set, queryAttributes.dummies)
  queryAttributes.forModel <- apply(queryAttributes.forModel[,-1],2,as.numeric)

  #make sure all dummy vars are accounted for
  if(any(is.na(match(names(preProcValues$mean), colnames(queryAttributes.forModel))))){
    newColNames <- names(preProcValues$mean)[which(is.na(match(names(preProcValues$mean),
                                                                 colnames(queryAttributes.forModel))))]
    dfAdd <- as.data.frame(matrix(data=0,nrow=nrow(queryAttributes.forModel),
                                   ncol=length(newColNames)))
    colnames(dfAdd) <- newColNames
    queryAttributes.forModel <- cbind(queryAttributes.forModel, dfAdd)
    queryAttributes.forModel <- queryAttributes.forModel[,match(names(preProcValues$mean),
                                                                  colnames(queryAttributes.forModel))]
  }

  #pre-process values
  queryAttributes.forModel <- predict(preProcValues, queryAttributes.forModel)
  queryAttributes.forModel <- cbind(queryAttributes.forModel, Class="UNKNOWN")
  queryAttributes.forModel <- as.data.frame(queryAttributes.forModel, stringsAsFactors=FALSE)
  queryAttributes.forModel$Class <- as.factor(queryAttributes.forModel$Class)
  for(n in 1:(length(colnames(queryAttributes.forModel)) -1)){
    queryAttributes.forModel[,n] <- as.numeric(queryAttributes.forModel[,n])
  }

  #SVM prediction feature
  newFeat <- kernlab::predict(branchpointer.svm,queryAttributes.forModel, "probabilities")
  newFeat <- newFeat[,1]

  queryAttributes.forModel <- cbind(queryAttributes.forModel, newFeat)

  #gbm prediction
  p <- predict(object=branchpointer.gbm, queryAttributes.forModel[,predictorNames],"prob")

  #reconfigure
  
  queryAttributes@elementMetadata$branchpoint_prob <- p[,1]

  #### U2 binding energy###
  m <- match(colnames(U2_binding_df)[-c(1:3)],colnames(queryAttributes@elementMetadata))

  if(useParallel){
    cluster <- parallel::makeCluster(cores)
    U2Eightmers <- parallel::parApply(cluster,as.data.frame(queryAttributes@elementMetadata[,m]),1,paste, collapse="")
    parallel::stopCluster(cluster)
  }else{
    #needs as.data.frame to treat rows as vectors
    U2Eightmers <- apply(as.data.frame(queryAttributes@elementMetadata[,m]),1,paste, collapse="")
  }

  x <- match(U2Eightmers, U2_binding_df$eightmers)
  queryAttributes@elementMetadata$U2_binding_energy <- U2_binding_df$energy[x]


  #branchpointPred$id <- stringr::str_sub(branchpointPred$id,1,-8)
  #branchpointPred <- branchpointPred[order(branchpointPred[,1], branchpointPred[,5], branchpointPred[,4]),]
  #colnames(branchpointPred)[which(colnames(branchpointPred) == "seq_pos0")] <- "nucleotide"

  return(queryAttributes)
}
