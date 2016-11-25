#' Predict branchpoint probability scores
#'
#' predicts branchpoint probability scores for each query site.
#' @param query_attributes data.frame produced by getBranchpointSequence
#' containing features for branchpoint prediction
#' @param use_parallel use parallelisation to speed up code?
#' (reccomended for > 500 query entries)
#' @param cores number of cores to use in parallelisation (default = \code{1})
#' @return data.frame with branchpoint probaility scores for each site in query_attributes
#' @export
#' @import caret
#' @importFrom kernlab predict
#' @importClassesFrom kernlab ksvm
#' @import gbm
#' @import parallel
#' @importFrom stringr str_sub
#' @examples
#' branchpoint_predicitons <- predictBranchpoints(query_attributes)
#' @author Beth Signal

predictBranchpoints <- function(query_attributes,
                             use_parallel=FALSE,
                             cores=1){
  if(use_parallel){

    max_cores <- parallel::detectCores()

    if(max_cores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", max_cores,")"))
      message(paste0("using all available cores"))
      cores <- max_cores

    }

  }

  query_attributes_formodel <- cbind(set="UNKNOWN", query_attributes)
  query_attributes_formodel <- query_attributes_formodel[,c(1,2,9:28)]

  #match to model variable names
  colnames(query_attributes_formodel) <-  HCN_names

  #remove any rows with NA values
  is_na <- apply(query_attributes_formodel[,-c(1,2)], 2, is.na)
  is_na_v <- apply(is_na,1, any)
  rm <- which(is_na_v == TRUE)
  if(length(rm) > 0){
    query_attributes_formodel <- query_attributes_formodel[-rm,]
    query_attributes <- query_attributes[-rm,]
  }

  #convert sequences to dummy vars
  query_attributes_dummies <- predict(dummies, newdata=query_attributes_formodel[,-2])
  query_attributes_formodel <- cbind(set=query_attributes_formodel$set, query_attributes_dummies)
  query_attributes_formodel <- apply(query_attributes_formodel[,-1],2,as.numeric)

  #make sure all dummy vars are accounted for
  if(any(is.na(match(names(preProcValues$mean), colnames(query_attributes_formodel))))){
    new_col_names <- names(preProcValues$mean)[which(is.na(match(names(preProcValues$mean),
                                                                 colnames(query_attributes_formodel))))]
    df_add <- as.data.frame(matrix(data=0,nrow=nrow(query_attributes_formodel),
                                   ncol=length(new_col_names)))
    colnames(df_add) <- new_col_names
    query_attributes_formodel <- cbind(query_attributes_formodel, df_add)
    query_attributes_formodel <- query_attributes_formodel[,match(names(preProcValues$mean),
                                                                  colnames(query_attributes_formodel))]
  }

  #pre-process values
  query_attributes_formodel <- predict(preProcValues, query_attributes_formodel)
  query_attributes_formodel <- cbind(query_attributes_formodel, Class="UNKNOWN")
  query_attributes_formodel <- as.data.frame(query_attributes_formodel, stringsAsFactors=FALSE)
  query_attributes_formodel$Class <- as.factor(query_attributes_formodel$Class)
  for(n in 1:(length(colnames(query_attributes_formodel)) -1)){
    query_attributes_formodel[,n] <- as.numeric(query_attributes_formodel[,n])
  }

  #SVM prediction feature
  newFeat <- kernlab::predict(branchpointer.svm,query_attributes_formodel, "probabilities")
  newFeat <- newFeat[,1]

  query_attributes_formodel <- cbind(query_attributes_formodel, newFeat)

  #gbm prediction
  p <- predict(object=branchpointer.gbm, query_attributes_formodel[,predictorNames],"prob")

  #reconfigure
  branchpoint_pred <- data.frame(id=query_attributes$id,branchpoint_prob=p[,1],
                                 seq_pos0=query_attributes$seq_pos0, distance=query_attributes$to_3prime)
  branchpoint_pred$allele_status <- stringr::str_sub(branchpoint_pred$id, -3,-1)

  m <- match(branchpoint_pred$id, query_attributes$id)
  branchpoint_pred <- cbind(branchpoint_pred,
                         chromosome=query_attributes$chromosome[m],
                         strand=query_attributes$strand[m],
                         end=query_attributes$end[m],
                         exon_3prime=query_attributes$exon_3prime[m],
                         exon_5prime=query_attributes$exon_5prime[m])


  #### U2 binding energy###
  m <- match(colnames(U2_binding_df)[-c(1:3)],colnames(query_attributes))

  if(use_parallel){
    cluster <- parallel::makeCluster(cores)
    u2_eightmers <- parallel::parApply(cluster,query_attributes[,m],1,paste, collapse="")
    parallel::stopCluster(cluster)
  }else{
    u2_eightmers <- apply(query_attributes[,m],1,paste, collapse="")
  }

  x <- match(u2_eightmers, U2_binding_df$eightmers)
  branchpoint_pred$U2_binding_energy <- U2_binding_df$energy[x]


  branchpoint_pred$id <- stringr::str_sub(branchpoint_pred$id,1,-8)
  branchpoint_pred <- arrange(branchpoint_pred, id,allele_status, distance)
  colnames(branchpoint_pred)[which(colnames(branchpoint_pred) == "seq_pos0")] <- "nucleotide"

  return(branchpoint_pred)
}
