#' Get the best polypyramidine tract
#'
#' Takes a query genomic sequence, finds all potential polypyramidine tracts (PPTs)
#' between the test site and the annotated 3'exon.
#' Returns the distance to the start of the longest PPT, and its length.
#' @param attributes_line line from a query attributes data.frame (or a character vector)
#' must have 3'SS distance at 3,
#' and nucleotide sequence at position 4.
#' @return distance to the start of the longest PPT, and its length
#' @export
#' @import plyr
#' @examples
#' attributes_line <- c("id", "1000","27", paste(c("A","G",rep(c("C","T"),4))[sample(1:10, 501, replace=TRUE)] , collapse=""))
#' pyra_df <- getPPT(attributes_line)
#' @author Beth Signal

getPPT <- function(attributes_line){

  dist_3prime <- as.numeric(attributes_line[3])

  #get sequence between the tested site and 3'exon
  seq <- substr(attributes_line[4], 251, (250 + dist_3prime))
  seq <- unlist(strsplit(seq, ""))

  pyramidines <- which(seq == "T" | seq == "C")

  if(length(pyramidines) >1){
    pyra_dist <- vector()

    #get distance to the next pyramidine
    for(p in 2:length(pyramidines)){
      pyra_dist[p-1] <- pyramidines[p] - pyramidines[p-1]
    }

    longest_run <- 0
    best_run <- 0
    best_start <- 0
    percent_pyra <- 0
    start <- 1
    start_vec <- vector()
    run_vec <- vector()
    percent_pyra_vec <- vector()

    #find best PPT
    for(p in seq_along(pyra_dist)){

      #For the first instance of at least 2 sequential pyramidines
      if(pyra_dist[p] == 1 & longest_run == 0){
        longest_run <- 1
        start <- p
        seq_p <- seq[(pyramidines[start]):(pyramidines[start+longest_run])]
        percent_pyra <- length(which(seq_p=="C"|seq_p=="T"))/length(seq_p)

        if(longest_run > best_run){
          best_run <- longest_run
          best_start <- start
          best_percent <- percent_pyra
        }


      #allowing gaps with one purine, grow the tract
      }else if(pyra_dist[p] == 1 | pyra_dist[p] == 2){
        longest_run <- longest_run + 1
        seq_p <- seq[(pyramidines[start]):(pyramidines[start+longest_run])]
        percent_pyra <- length(which(seq_p == "C"|seq_p == "T"))/length(seq_p)

        if(longest_run > best_run){
          best_run <- longest_run
          best_start <- start
          best_percent <- percent_pyra
        }

      #once too many purines are encontered, save tract
      }else{
        start_vec <- append(start_vec, best_start)
        run_vec <- append(run_vec, longest_run)
        percent_pyra_vec <- append(percent_pyra_vec, percent_pyra)
        longest_run <- 0
        start <- p + 1
        best_start <- start
      }
    }

    #add the last PPT
    start_vec <- append(start_vec, start)
    run_vec <- append(run_vec, longest_run)
    percent_pyra_vec <- append(percent_pyra_vec, percent_pyra)

    #make into data.frame
    pyra_df <- data.frame(set=c(1:length(run_vec)), run_vec,start_vec,percent_pyra_vec)
    pyra_df <- plyr::arrange(pyra_df, plyr::desc(run_vec))

    best_start <- pyra_df$start_vec[which.max(pyra_df$run_vec)]
    best_run <- pyra_df$run_vec[which.max(pyra_df$run_vec)]

    #If PPTs are long enough, shorten to get pyramidine content above 80%
    pyra_df <- pyra_df[pyra_df$run_vec >=10,]
    if(dim(pyra_df)[1] != 0){
      run_extend <- vector()
      start_extend <- vector()
      polyper_extend <- vector()

      for(py in seq_along(length(pyra_df[,1]))){
        seq_p <- seq[(pyramidines[pyra_df$start_vec[py]]):(pyramidines[pyra_df$start_vec[py]+pyra_df$run_vec[py]])]
        startp <- 1

        run <- length(seq_p)

        polyper_extend[py] <- pyra_df$percent_pyra_vec[py]
        while(polyper_extend[py] < 0.8 & run >= 10){
          purine <- which(seq_p == "A"| seq_p == "G")
          purine_neg <- (length(seq_p) +1 -purine)*-1
          purine_both <- c(purine,purine_neg)
          min_purine <- which.min(abs(purine_both))

          point <- purine_both[min_purine]

          if(point < 0){
            rm <- length(seq_p)-point*-1
            new_seq <- seq_p[1:rm]
            run <- length(new_seq)
          }else{
            new_seq <- seq_p[(point+1):length(seq_p)]
            run <- length(new_seq)
            startp <- startp+point
          }
          polyper_extend[py] <- length(which(new_seq == "C" | new_seq == "T")) / length(new_seq)
          seq_p <- new_seq
        }

        start_extend[py] <- (pyramidines[pyra_df$start_vec[py]])+startp-1
        run_extend[py] <- run
      }

      best_start <- start_extend[which.max(run_extend)]
      best_run <- run_extend[which.max(run_extend)]
    }


    line <- c(best_start,best_run)

    #if only 1 pyramidine
  }else if(length(pyramidines) == 1){
    line <- c(pyramidines, 1)
    #if no pyramidines
  }else{
    line <- c(0,0)
  }
  return(line)
}
