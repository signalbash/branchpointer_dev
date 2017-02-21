#' Predict branchpoint probability scores
#'
#' predicts branchpoint probability scores for each query site.
#' @param queryAttributes data.frame produced by getBranchpointSequence
#' containing features for branchpoint prediction
#' @param useParallel use parallelisation to speed up code?
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
#' branchpointPredictions <- predictBranchpoints(query_attributes)
#' @author Beth Signal

predictBranchpoints <- function(queryAttributes,
                             useParallel=FALSE,
                             cores=1){
  if(useParallel){

    maxCores <- parallel::detectCores()

    if(maxCores < cores){

      message(paste0("specified cores (", cores,") is greater than available cores(", maxCores,")"))
      message(paste0("using all available cores"))
      cores <- maxCores

    }

  }

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
    queryAttributes.forModel <- cbind(queryAttributes.forModel, df_add)
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
