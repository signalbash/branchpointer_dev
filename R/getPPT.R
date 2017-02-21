#' Get the best polypyramidine tract
#'
#' Takes a query genomic sequence, finds all potential polypyramidine tracts (PPTs)
#' between the test site and the annotated 3'exon.
#' Returns the distance to the start of the longest PPT, and its length.
#' @param attributesLine line from a query attributes data.frame (or a character vector)
#' must have 3'SS distance at 3,
#' and nucleotide sequence at position 4.
#' @return distance to the start of the longest PPT, and its length
#' @export
#' @import plyr
#' @examples
#' attributes_line <- c("id", "1000","27", paste(c("A","G",rep(c("C","T"),4))[sample(1:10, 501, replace=TRUE)] , collapse=""))
#' pyra_df <- getPPT(attributes_line)
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
