#' Heat-map of TIPC metric color-coded by clustering results
#'
#' Plotting heat-maps of TIPC metric clustering results of all k's unless it is
#' specified; using R package \code{\link[ComplexHeatmap]{ComplexHeatmap}}
#'
#' @param root_dir A directory path containing (1) the TIPC_metrics, i.e.
#'   normalized TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; (2) the result folder of TIPC
#'   clustering results output from \code{\link[TIPC]{consensus_clustering}}.
#' @param clustering_subfolder_nm A character string of sub-folder name generated
#'   during \code{\link[TIPC]{consensus_clustering}}, parked under \code{root_dir}.
#' @param one_k An integer specifying the cluster number selected for heat-map
#'   plotting; otherwise all k's found under \code{TIPC_cluster_dir} will be
#'   processed.
#'
#' @examples
#' root_dir <- "C:/Users/Mai Chan Lau/Desktop/TIPC_package/test_run/TIPC_hexLen100/"
#' clustering_subfolder_nm <- 'ConsensusClusterPlus_test'
#' plot_TIPC_heatmap(root_dir = root_dir)
#' @export
#' @importFrom grDevices pdf dev.off
#' @importFrom grid gpar unit
#' @importFrom utils read.csv
#' @importFrom randomcoloR distinctColorPalette
#' @import ComplexHeatmap
plot_TIPC_heatmap <- function(root_dir =  NULL, one_k = NULL, clustering_subfolder_nm = 'test') {

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  TIPC_cluster_dir <- file.path(root_dir, 'ConsensusClusterPlus_test')
  ## ======================
  ## load TIPC_metrics
  ## ======================
  TIPC_metrics_holder <- load(file.path(root_dir,'TIPC_metrics.RData'))
  TIPC_metrics = get(TIPC_metrics_holder)

  ## ======================
  ## find all k's folders
  ## ======================
  sub_folder_nms <- list.dirs(path = TIPC_cluster_dir, full.names = TRUE)
  ## subsetting for folders with names containing 'k'
  sub_folder_nms <- grep(x=sub_folder_nms, pattern = 'k[0-9]', value = TRUE)
  if(!is.null(one_k)) sub_folder_nms <- grep(x=sub_folder_nms, pattern = paste0('k',one_k,'$'), value = TRUE)

  ## ======================
  ## palette selection: distinct color coding
  ## ======================
  k_range <- sapply(strsplit(x=sub_folder_nms, split = '/k'), "[[", 2)
  set.seed(123)
  col_vec <- distinctColorPalette(max(k_range))

  ## ======================
  ## loop over each clustering results with different k's
  ## ======================
  for (kk in sub_folder_nms){
    cat(kk, '\n')

    ## ---------------
    ## load cluster results
    ## ---------------
    filenm <- list.files(path = kk, pattern = 'cluster_no_k')
    if(length(filenm) > 1) stop('More than 1 clustering result file is found\n')

    TIPC_clusters <- read.csv(file.path(kk,filenm), row.names = NULL, as.is = TRUE)

    ## ---------------
    ## order by clustering number ascendingly
    ## ---------------
    TIPC_clusters <- TIPC_clusters[order(TIPC_clusters$cluster_no),]

    ## ---------------
    ## check if tumor_ids consistency
    ## ---------------
    if(sum(rownames(TIPC_metrics) %in% TIPC_clusters$tumor) != nrow(TIPC_clusters))
      stop('Some tumor ids in TIPC metrics are not found in TIPC clustering results!\n')

    ## ---------------
    ## matching tumor ids ordering: TIPC metrics and TIPC clusters
    ## ---------------
    TIPC_metrics <- TIPC_metrics[match(TIPC_clusters$tumor, rownames(TIPC_metrics)),]
    if(! identical(rownames(TIPC_metrics), TIPC_clusters$tumor) )
      stop('Tumor ids between TIPC metrics and TIPC clusters cannot be matched\n')

    ## ---------------
    ## TIPC metric scaling
    ## ---------------
    scaled_TIPC_metrics <- as.data.frame(scale(TIPC_metrics))

    ## ---------------
    ## naming color vect with TIPC cluster ids
    ## ---------------
    col_vec_sub <- col_vec[1:max(TIPC_clusters$cluster_no)]
    names(col_vec_sub) <- seq(min(TIPC_clusters$cluster_no):max(TIPC_clusters$cluster_no))

    ## ---------------
    ## column annotation of TIPC cluster ids
    ## ---------------
    ha_column = HeatmapAnnotation(df = data.frame(CLUSTER = TIPC_clusters$cluster_no),
                                  col=list(CLUSTER = col_vec_sub),
                                  gap = unit(1.5, "mm"), height = unit(2,"cm"),
                                  annotation_legend_param = list(CLUSTER = list(title = "CLUSTER", title_gp = gpar(fontsize = 18),
                                                                                labels_gp = gpar(fontsize = 14))))

    ## ---------------
    ## Calling heatmap plotting: scaled TIPC metrics
    ## ---------------
    h1 = Heatmap(t(scaled_TIPC_metrics), name = "scaled", row_title = "TIPIC 6-metrics",
                  column_title = "tumor ids",
                  column_dend_reorder = FALSE,
                  top_annotation = ha_column,cluster_columns=FALSE,
                  show_column_names = FALSE,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                              legend_direction = "horizontal",
                                              legend_height = unit(8, "cm"),legend_width = unit(3, "cm"),
                                              title_position = "lefttop"))

    fileNm <- file.path(kk,'heatmap_TIPCmetrics_scaled.pdf')
    pdf(fileNm)
    draw(h1,heatmap_legend_side = "bottom")
    dev.off()

    ## ---------------
    ## Calling heatmap plotting: raw metrics
    ## ---------------
    h2 = Heatmap(t(TIPC_metrics), name = "raw", row_title = "TIPIC 6-metrics",
                  column_title = "tumor ids",
                  column_dend_reorder = FALSE,top_annotation = ha_column,cluster_columns=FALSE,
                  show_column_names = FALSE,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                              legend_direction = "horizontal",
                                              legend_height = unit(8, "cm"),legend_width = unit(3, "cm"),
                                              title_position = "lefttop"))

    fileNm <- file.path(kk,'heatmap_TIPCmetrics_raw.pdf')
    pdf(fileNm)
    draw(h2,heatmap_legend_side = "bottom")
    dev.off()
  }# end all k's


}
