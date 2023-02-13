#' Clustering of TIPC metrics
#'
#' Clustering of TIPC metrics using
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @param root_dir A directory path containing the TIPC_metrics, i.e. normalized
#'   TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; nrow = total. of tumors, ncol = 6
#'   (TIPC metrics) x 5 (directions).
#' @param min_k An integer indicating the minimum number of clusters.
#' @param max_k An integer indicating the maximum number of clusters.
#' @param output_bnm A character string appended to output folder name;
#'   sub-folders are created for different k from min_k to max_k.
#' @examples
#' root_dir <- "C:/Users/Mai Chan Lau/Desktop/TIPC_package/test_run/TIPC_hexLen100"
#' consensus_clustering(root_dir=root_dir)
#' @export
#' @importFrom grDevices pdf
#' @importFrom utils write.csv
consensus_clustering <- function(min_k = 2, max_k = 6,
                                 root_dir = NULL, output_bnm = 'test'){

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  ## ======================
  ## load TIPC_metrics
  ## ======================
  TIPC_metrics_holder <- load(file.path(root_dir,'TIPC_metrics.RData'))
  TIPC_metrics = get(TIPC_metrics_holder)

  ## ======================
  ## create output directory
  ## ======================
  output_dir_bnm <- paste0('ConsensusClusterPlus_',output_bnm)
  res_subdir <- file.path(root_dir, output_dir_bnm)
  dir.create(res_subdir)

  ## ======================
  ## data scaling
  ## ======================
  df <- as.data.frame(scale(TIPC_metrics))

  ## ======================
  ## caling ConsensusClusterPlus
  ## ======================
  clustering_res = ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(t(df)),maxK=max_k,reps=50,
                                                              pItem=0.8,pFeature=1,title=res_subdir,
                                                              distance="pearson",seed=123,plot="pdf")

  setwd(res_subdir)
  icl = ConsensusClusterPlus::calcICL(clustering_res,title='cluster_item_consensus_plots',plot="pdf")

  for (k in c(min_k:max_k)){
    res_k <- data.frame('tumor'=names(clustering_res[[k]][["consensusClass"]]),
                        'cluster_no'=clustering_res[[k]][["consensusClass"]])

    setwd(res_subdir)
    sub_folderNm <- paste0('k',k)
    dir.create(sub_folderNm)
    setwd(sub_folderNm)
    write.csv(res_k, file=paste0('cluster_no_k',k,'.csv'), row.names = FALSE)

  }


}
