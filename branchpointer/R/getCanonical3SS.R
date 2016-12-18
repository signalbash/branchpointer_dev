#' Get locations of the first five AG 3'Splice site motifs
#'
#' Takes a variable length vector of distances to the AG motif,
#' sorts and returns the first five.
#' If there are less than five elements in the vector, returns the sorted vector
#' and fills the remainder of the values with 300.
#' @param ag Vector of distances to the AG splice site motif.
#' @return Locations of the first five AG dinucleotides
#' @export
#' @examples
#' sequence <- paste(c("A","T","C","G")[sample(1:4, 501, replace=TRUE)] , collapse="")
#' AG_seq_positions <- gregexpr("AG",sequence, 252,501)
#' canon_hits <- getCanonical3SS(AG_seq_positions[[1]])
#' @author Beth Signal

getCanonical3SS <- function(ag){
  ag <- sort(ag)
  #get first five matches
  ag <- ag[1:5]
  #if less than 5 matches, replace remainder with 300
  ag[is.na(ag)] <- 300

  return(ag)
}
