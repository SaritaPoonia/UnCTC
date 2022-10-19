#------------------------------------------------------------------------
#-------------------- Plots DDLK clusters   ----------------------------
#------------------------------------------------------------------------

#' Pathway based principal component analysis plots
#' @description unCTC UMAP and PCA plots
#'
#' @param Pathway_score Pathway score matrix. Row names are pathway/ genesets
#' names and column names are samples/ cells
#' @param Pathway_metadata PathwayDDLK_clust output from DDLK_Clust function
#' @param colorby color by any column name from Pathway_metadata, if pathway
#' metadata is not given then default is "Data_id"
#' @param Color_cluster Any column name from Pathway_metadata,
#' default is "Clusters"
#' @param pairsplotLegend Legend position for pairsplot.
#'
#' @importFrom cowplot plot_grid
#' @importFrom PCAtools pairsplot
#' @importFrom PCAtools pca
#' @import magrittr
#' @import ggplot2
#' @import umap
#'
#' @return plots list of color by class, color by clusters, pairsplot
#'
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' Genesets = unCTC::c2.all.v7.2.symbols
#' Pathway_score = PathwayEnrichmentScore(data_list=Data_list,
#'                                         data_id= Data_Id,
#'                                         Genesets=Genesets,
#'                                         min.size=70,
#'                                         max.size=100)
#'
#' DDLK_Clusters = DDLK_Clust(PathwayScore = Pathway_score$Pathway_score,
#'                            PathwayMetaData=Pathway_score$Pathway_metadata,
#'                             n=3,
#'                             out.dir = paste0(getwd(),"/unCTC"))
#'
#' PathwayScore = DDLK_Clusters$Pathway_score
#' PathwayMetadata = DDLK_Clusters$PathwayDDLK_clust
#' Plots_output = unCTC_pathway_plots(Pathway_score=PathwayScore,
#'                            Pathway_metadata = PathwayMetadata,
#'                            colorby = "Data_id",
#'                            Color_cluster = "Clusters")
#' @export
unCTC_pathway_plots = function(Pathway_score,
                       Pathway_metadata,
                       colorby = "Data_id",
                       Color_cluster = "Clusters",
                       pairsplotLegend = c("left","right","none"))
{
 #Colorkey to visualize clusters
  ColorKeyDataID = c("peru","steelblue","darkolivegreen4","palevioletred4",
                     "darkcyan","darkorchid4","darkslategray","firebrick1",
                     "salmon3","paleturquoise1","mediumaquamarine",
                     "greenyellow","black","deepskyblue3","mediumblue",
                     "darkred","gold","gray50","hotpink","khaki3",
                     "yellow4","lavender","cornsilk4","orchid4",
                     "yellow3", "darkgreen","skyblue1","khaki4",
                     "tan4","pink")

  # Match sample names in Pathway_score and Pathway_metadata
  ins = intersect(colnames(Pathway_score),rownames(Pathway_metadata))

  # Order cells of pathway matrix according to pathway meta data
  Pathway_score = Pathway_score[,ins]
  Pathway_metadata = Pathway_metadata[ins,]

  # pca output of pathway_score matrix. For details kindly see pcatools::pca()
  p1 <- pca(Pathway_score, metadata = Pathway_metadata, removeVar = 0.1)

  #umap output of pathway matrix
  umap_path = umap(t(Pathway_score),n_neighbors= 50,
                   n_components= 2,metric= "euclidean",n_epochs= 100,
                   min_dist = 0.1,init ="spectral",random_state=123,
                   alpha=2.5,gamma=2.5,
                   negative_sample_rate=5)

  # Make data frame with UMAP embeddings and metadata
  umap_df = data.frame(UMAP1 = umap_path$layout[,1],
             UMAP2 = umap_path$layout[,2],
             Pathway_metadata)


  #Draw multiple bi-plots of PC1 to PC5
  pairplot = pairsplot(p1,colby=colorby,colkey = ColorKeyDataID,
                       legendPosition = pairsplotLegend)


  #----------------------- plot ---------------------------------
  #Make dataframe from PCA embeddings and metadata
  gsva_pca_df = data.frame(PC1 = p1$rotated$PC1,
                           PC2 = p1$rotated$PC2,
                           Pathway_metadata)
  rownames(gsva_pca_df) = rownames(Pathway_metadata)

  #PC1 and PC2 class wise plot
  labels = gsva_pca_df[,colorby]
  Plot_class = gsva_pca_df %>%
    ggplot(aes(x = PC1, y = PC2, col = labels)) +
    labs(title = paste0("PCA plot group by",colorby))+
    geom_point(size = 2, stroke = 0.2, shape = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="right") +
    theme(legend.text = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_x_continuous(minor_breaks = seq(-10, 10, 5)) +
    scale_y_continuous(minor_breaks = seq(-10, 10, 5))+
    scale_colour_manual(values=ColorKeyDataID)

  #UMAP1 and UMAP2 class wise plot
  labels = gsva_pca_df[,colorby]
  Plot_umap_class = umap_df %>%
    ggplot(aes(x = UMAP1, y = UMAP2, col = labels)) +
    labs(title = paste0("UMAP plot group by",colorby))+
    geom_point(size = 2, stroke = 0.2, shape = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="right") +
    theme(legend.text = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_x_continuous(minor_breaks = seq(-10, 10, 5)) +
    scale_y_continuous(minor_breaks = seq(-10, 10, 5))+
    scale_colour_manual(values=ColorKeyDataID)

  #PC1 and PC2 cluster wise plot
  Cluster_labels = gsva_pca_df[,Color_cluster]
  Plot_cluster = gsva_pca_df %>%
    ggplot(aes(x = PC1, y = PC2, col = Cluster_labels)) +
    labs(title = paste0("PCA plot group by ",
                        Color_cluster))+
    geom_point(size = 2, stroke = 0.2, shape = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="right") +
    theme(legend.text = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_x_continuous(minor_breaks = seq(-10, 10, 5)) +
    scale_y_continuous(minor_breaks = seq(-10, 10, 5))+
    scale_colour_manual(values=ColorKeyDataID)


  #UMAP1 and UMAP2 cluster wise plot
  Cluster_labels = gsva_pca_df[,Color_cluster]
  Plot_umap_cluster = umap_df %>%
    ggplot(aes(x = UMAP1, y = UMAP2, col = Cluster_labels)) +
    labs(title = paste0("UMAP plot group by ",Color_cluster))+
    geom_point(size = 2, stroke = 0.2, shape = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="right") +
    theme(legend.text = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title=element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_x_continuous(minor_breaks = seq(-10, 10, 5)) +
    scale_y_continuous(minor_breaks = seq(-10, 10, 5))+
    scale_colour_manual(values=ColorKeyDataID)


  plots <- list(
                group_by_Class_PCA = Plot_class,
                group_by_Cluster_PCA =Plot_cluster,
                group_by_Class_umap=Plot_umap_class,
                group_by_Cluster_umap =Plot_umap_cluster,
                p1_p5_Pairsplot = pairplot
                )
  return(plots)
}
