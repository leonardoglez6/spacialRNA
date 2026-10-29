##%######################################################%##
#                                                          #
####              Descarga el set de datos              ####
#                                                          #
##%######################################################%##

set.seed(6)

setwd("/home/lgonzalez/bioinfo3/spacialRNA")

# Usaremos el set de datos de cerebro de ratón de la tecnología Visium,
# disponible en: https://www.10xgenomics.com/datasets/mouse-kidney-section-coronal-1-standard-1-1-0

# Para descargar, ejecuta estos comandos en la terminal de R:
#dir.create("visium_data/")

# Descomprime los archivos

#untar(tarfile = "visium_data/V1_Mouse_Kidney_raw_feature_bc_matrix.tar.gz", 
#      exdir = "visium_data/")

#untar(tarfile = "visium_data/V1_Mouse_Kidney_spatial.tar.gz", 
#      exdir = "visium_data/")


##%######################################################%##
#                                                          #
####       Configura las instrucciones de Giotto        ####
#                                                          #
##%######################################################%##

library(Giotto)

instructions <- createGiottoInstructions(save_dir = "visium_plots",
                                         save_plot = TRUE,
                                         show_plot = TRUE,
                                         return_plot = FALSE,
                                         python_path = NULL # alternativamente, "/local/python/path/python" 
                                         )

##%######################################################%##
#                                                          #
####               Crea el objeto Giotto                ####
#                                                          #
##%######################################################%##

# data_path <- "visium_data/" #Aauí es donde va a guardar los resultados

# visium_kidney <- createGiottoVisiumObject(visium_dir = data_path,
#                                         expr_data = "raw",
#                                         png_name = "tissue_lowres_image.png",
#                                         gene_column_index = 2,
#                                         instructions = instructions)


##%######################################################%##
#                                                          #
####           Alternativa del objeto Giotto            ####
#                                                          #
##%######################################################%##

library(data.table)

expression_matrix <- get10Xmatrix("visium_data/raw_feature_bc_matrix/",
                                  gene_column_index = 2)

######################################################################
######################################################################
# Con ayuda de Claude se crearon los archivos utilizados en la clase #
######################################################################
#pos <- fread("visium_data/spatial/tissue_positions_list.csv", header = FALSE,
#             col.names = c("barcode", "in_tissue", "array_row", "array_col",
#                           "pxl_row", "pxl_col"))

# Mismo orden que las columnas de la matriz (Giotto exige que coincidan)
#pos <- pos[match(colnames(expression_matrix), barcode)]
#stopifnot(!anyNA(pos$barcode))

# spatial_locs.csv
#spatial_locs <- data.frame(cell_ID = pos$barcode,
#                           sdimx   = pos$pxl_col,
#                           sdimy   = -pos$pxl_row)   # signo negativo para que no salga volteado
#write.csv(spatial_locs, "spatial_locs.csv", row.names = FALSE)

# cell_metadata.csv (la primera columna se usa como row.names en tu código)
#cell_metadata <- data.frame(cell_ID   = pos$barcode,
#                            in_tissue = pos$in_tissue,
#                            array_row = pos$array_row,
#                            array_col = pos$array_col)
#write.csv(cell_metadata, "cell_metadata.csv", row.names = FALSE)
######################################################################
######################################################################

spatial_locs <- read.csv("spatial_locs.csv",
                         header = TRUE)

visium_kidney <- createGiottoObject(expression = expression_matrix,
                                   spatial_locs = spatial_locs,
                                   instructions = instructions)

cell_metadata <- read.csv("cell_metadata.csv",
                          header = TRUE)

visium_kidney <- addCellMetadata(visium_kidney,
                                new_metadata = cell_metadata,
                                by_column = TRUE)

##%######################################################%##
#                                                          #
####                 Explora el objeto                  ####
#                                                          #
##%######################################################%##

# Visión general del objeto
visium_kidney

# Expresión
showGiottoExpression(visium_kidney)

# Metadatos de célula
showGiottoCellMetadata(visium_kidney)
pDataDT(visium_kidney)

getCellMetadata(visium_kidney)
pDataDT(visium_kidney)

# Metadatos de los genes/proteínas
showGiottoFeatMetadata(visium_kidney)

# Imágenes
showGiottoImageNames(visium_kidney)

# Reducción de dimensiones
showGiottoDimRed(visium_kidney)

# Instrucciones
instructions(visium_kidney)

##%######################################################%##
#                                                          #
####          Filtra los spots sobre el tejido          ####
#                                                          #
##%######################################################%##

spatPlot2D(gobject = visium_kidney, 
           cell_color = "in_tissue", 
           point_size = 1,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = FALSE, 
           image_name = "image")

cell_metadata <- pDataDT(visium_kidney)
in_tissue_barcodes <- cell_metadata[in_tissue == 1]$cell_ID

visium_kidney <- subsetGiotto(visium_kidney, 
                             cell_ids = in_tissue_barcodes)

spatPlot2D(gobject = visium_kidney, 
           cell_color = "in_tissue", 
           point_size = 1,
           cell_color_code = c("0" = "lightgrey", "1" = "blue"), 
           show_image = FALSE, 
           image_name = "image")

##%######################################################%##
#                                                          #
####                 Control de calidad                 ####
#                                                          #
##%######################################################%##

# Crea un objecto con el número de genes y transcritos por célula/spot
visium_kidney_statistics <- addStatistics(gobject = visium_kidney,
                                         expression_values = "raw")

spatPlot2D(gobject = visium_kidney_statistics, 
           cell_color = "nr_feats", 
           point_size = 1,
           color_as_factor = FALSE)

# Distribución de genes por spot/célula
filterDistributions(gobject = visium_kidney_statistics, 
                    detection = "cells")

# Distribución del número de células con cierto número de genes
filterDistributions(gobject = visium_kidney_statistics, 
                    detection = "feats")

# Evalúa diferentes puntos de corte al filtrar
filterCombinations(gobject = visium_kidney_statistics,
                   expression_thresholds = c(1, 2, 3),
                   feat_det_in_min_cells = c(50, 100, 200),
                   min_det_feats_per_cell = c(500, 1000, 1500))

##%######################################################%##
#                                                          #
####                      Filtrado                      ####
#                                                          #
##%######################################################%##

visium_kidney <- filterGiotto(gobject = visium_kidney,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                             expression_values = "raw",
                             verbose = TRUE)

##%######################################################%##
#                                                          #
####                   Normalización                    ####
#                                                          #
##%######################################################%##

visium_kidney <- normalizeGiotto(gobject = visium_kidney, 
                                scalefactor = 6000, 
                                verbose = TRUE)

visium_kidney

##%######################################################%##
#                                                          #
####                    Estadísticas                    ####
#                                                          #
##%######################################################%##

visium_kidney <- addStatistics(gobject = visium_kidney)

spatPlot2D(gobject = visium_kidney, 
           show_image = FALSE, 
           point_alpha = 0.7,
           cell_color = "nr_feats", 
           color_as_factor = FALSE)

##%######################################################%##
#                                                          #
####                        HVF                         ####
#                                                          #
##%######################################################%##

# Regresión de Loess es usada cuando la relación entre la expresión promedio
# y la varianza son no lineales o pueden ser descritas por un modelo no paramétrico

visium_kidney <- calculateHVF(gobject = visium_kidney, 
                             method = "cov_loess",
                             HVFname = paste0("hvf_cov_loess"),
                             save_plot = TRUE)

# Residuales Pearson son usados para la estabilización de la varianza 
# (para lidiar con el ruido técnico) y resaltar los genes con alta dispersión

visium_kidney <- calculateHVF(gobject = visium_kidney, 
                             method = "var_p_resid",
                             HVFname = paste0("hvf_var_p_resid"),
                             save_plot = TRUE)

# Grupos de covarianza (default), son usagos cuando la variabilidad de la expresión 
# difiere a través de los niveles de expresión o las regiones espaciales, sin
# asumir una relación específica entre la expresión promedio y la varianza.

visium_kidney <- calculateHVF(gobject = visium_kidney, 
                             method = "cov_groups",
                             HVFname = paste0("hvf_cov_groups"),
                             save_plot = TRUE)

metodos <- c("cov_groups", "cov_loess", "var_p_resid")

fm <- fDataDT(visium_kidney)
names(fm)   # ¿existen hvf_cov_groups, hvf_cov_loess, hvf_var_p_resid?

hvf <- lapply(metodos, function(m) {
  col <- paste0("hvf_", m)
  if (!col %in% names(fm)) return(NA_character_)   # avisa en vez de fallar
  fm$feat_ID[fm[[col]] == "yes"]
})
names(hvf) <- metodos
sapply(hvf, length)   # número de HVG por método

visium_kidney <- calculateHVF(gobject   = visium_kidney,
                              method    = "cov_loess",
                              save_plot = TRUE)

##%######################################################%##
#                                                          #
####              Reducción de dimensiones              ####
#                                                          #
##%######################################################%##

visium_kidney <- runPCA(gobject = visium_kidney)

screePlot(gobject = visium_kidney, 
          ncp = 30)

plotPCA(gobject = visium_kidney)

##%######################################################%##
#                                                          #
####                        UMAP                        ####
#                                                          #
##%######################################################%##

visium_kidney <- runUMAP(visium_kidney, 
                        dimensions_to_use = 1:10)

plotUMAP(gobject = visium_kidney)

##%######################################################%##
#                                                          #
####                        tSNE                        ####
#                                                          #
##%######################################################%##

visium_kidney <- runtSNE(visium_kidney, 
                        dimensions_to_use = 1:10)

plotTSNE(gobject = visium_kidney)

##%######################################################%##
#                                                          #
####                  Vecinos cercanos                  ####
#                                                          #
##%######################################################%##

# sNN network (default)
visium_kidney <- createNearestNetwork(gobject = visium_kidney, 
                                     dimensions_to_use = 1:10, 
                                     k = 15)

# K-Nearest Neighbors (kNN)
visium_kidney <- createNearestNetwork(gobject = visium_kidney, 
                                     dimensions_to_use = 1:10, 
                                     k = 15,
                                     type = "kNN")

##%######################################################%##
#                                                          #
####                 Leiden clustering                  ####
#                                                          #
##%######################################################%##

visium_kidney <- doLeidenCluster(gobject = visium_kidney, 
                                resolution = 0.4, 
                                n_iterations = 1000)

plotPCA(gobject = visium_kidney,
        cell_color = "leiden_clus")

plotUMAP(gobject = visium_kidney,
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2.5)

plotTSNE(gobject = visium_kidney,
         cell_color = "leiden_clus", 
         show_NN_network = TRUE, 
         point_size = 2.5)

spatPlot2D(visium_kidney, 
           cell_color = "leiden_clus",
           point_size = 1,
           coord_fix_ratio = 1)

# Podemos resaltar un grupo de interés
spatPlot2D(visium_kidney, 
           cell_color = "leiden_clus",
           select_cell_groups = "5", 
           point_size = 1,
           coord_fix_ratio = 1, 
           show_other_cells = TRUE,
           cell_color_code = c("5" = "steelblue3"), 
           other_cell_color = "grey")

spatPlot2D(visium_kidney, 
           cell_color = "leiden_clus",
           select_cell_groups = "3", 
           point_size = 1,
           coord_fix_ratio = 1, 
           show_other_cells = TRUE,
           cell_color_code = c("3" = "yellow"), 
           other_cell_color = "grey")

##%######################################################%##
#                                                          #
####              Integración con scRNAseq              ####
#                                                          #
##%######################################################%##

data_path <- "/home/lgonzalez/spacial_RNA"

## Descarga el dataset de referencia
GiottoData::getSpatialDataset(dataset = "scRNA_mouse_brain", 
                              directory = data_path)

data_path <- "/home/lgonzalez/spacial_RNA/scRNA_mouse_brain"

## Lee el scRNAseq como objeto Giotto
sc_expression <- file.path(data_path, "brain_sc_expression_matrix.txt.gz")
sc_metadata <- file.path(data_path, "brain_sc_metadata.csv")

g_sc <- createGiottoObject(expression = sc_expression,
                           instructions = instructions)

## Agrega los metadatos
cell_metadata <- data.table::fread(sc_metadata)

g_sc <- addCellMetadata(g_sc, 
                        new_metadata = cell_metadata[, c("Class", "Subclass")])

## Normalización
g_sc <- normalizeGiotto(g_sc)

##%######################################################%##
#                                                          #
####         Enriquecimiento de tipos celulares         ####
#                                                          #
##%######################################################%##

## Encuentra los genes que caracterizan el comportamiento de cada tipo celular 
## en el scRNAseq de referencia

library(scran)

markers_scran <- findMarkers_one_vs_all(gobject = g_sc, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "Class", 
                                        min_feats = 3)

topgenes_scran <- markers_scran[, head(.SD, 10), by = "cluster"]

celltypes <- levels(factor(markers_scran$cluster)) 

sign_list <- list()
for (i in 1:length(celltypes)){
    sign_list[[i]] = topgenes_scran[which(topgenes_scran$cluster == celltypes[i]),]$feats
}

sign_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)

visium_kidney <- runPAGEEnrich(gobject = visium_kidney, 
                              sign_matrix = sign_matrix)

spatCellPlot2D(gobject = visium_kidney,
               spat_enr_names = "PAGE",
               cell_annotation_values = celltypes,
               cow_n_col = 3,
               coord_fix_ratio = 1, 
               point_size = 1, 
               show_legend = TRUE)

##%######################################################%##
#                                                          #
####                   Deconvolución                    ####
#                                                          #
##%######################################################%##

library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)

# Crea el objeto referencia desde el scRNAseq 
sc_counts <- getExpression(g_sc,
                           values = "raw",
                           output = "matrix")

cell_metadata <- pDataDT(g_sc)
sc_cell_types <- factor(cell_metadata$Class)
names(sc_cell_types) <- cell_metadata$cell_ID

reference_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = sc_counts),
  colData = data.frame(cell_type = sc_cell_types)
)

# Crea el objeto espacial desde el spatial transcriptomics
spatial_locs <- getSpatialLocations(visium_kidney,
                                    output = "data.table")

st_coords <- data.frame(x = spatial_locs$sdimx,
                        y = spatial_locs$sdimy)
rownames(st_coords) <- spatial_locs$cell_ID

st_counts <- getExpression(visium_kidney,
                           values = "raw",
                           output = "matrix")

spatial_spe <- SpatialExperiment::SpatialExperiment(
  assay = st_counts,
  spatialCoords = as.matrix(st_coords)
)

# Crea y ejecuta RCTD
myRCTD <- createRctd(spatial_experiment = spatial_spe, 
                     reference_experiment = reference_se
)

results_spe <- runRctd(myRCTD, rctd_mode = "full", max_cores = 1)

# Da formato a los resultados
results_full <- results_spe@assays@data$weights

x <- as.data.frame(t(as.matrix(results_full)))
norm_weights <- sweep(x, 1, rowSums(x), '/') 
norm_weights <- cbind(data.frame(cell_ID = rownames(results_spe@colData)),
                      norm_weights)


# Visualiza el resultado
library(ggplot2)
library(scatterpie)

rctd_proportions <- merge(spatial_locs, norm_weights)

ggplot() +
  geom_scatterpie(aes(sdimx, sdimy),
                  data = rctd_proportions,
                  cols = colnames(x),
                  pie_scale = 0.5,
                  linewidth = 0) +
  labs(fill = "") +
  theme_minimal()
ggsave("visium_results/deconvolution_plot.png")

##%######################################################%##
#                                                          #
####           Genes con patrones espaciales            ####
#                                                          #
##%######################################################%##

# createSpatialNetwork() crea una red espacial que conecta a las células/spots
# basados en su distancia física. Existen dos métodos: Delaunay y kNN.
# En el método Delaunay, los vecinos son elegidos por la triangulación Delaunay
# y un criterio de máxima distancia.
# En el método kNN, el número de vecinos puede ser determinado por el parámetro k,
# o por la máxima distancia entre células.

visium_kidney <- createSpatialNetwork(gobject = visium_kidney,
                                     method = "kNN", 
                                     k = 6,
                                     maximum_distance_knn = 400,
                                     name = "spatial_network")

spatPlot2D(gobject = visium_kidney,  
           show_network= TRUE,
           network_color = "blue", 
           spatial_network_name = "spatial_network")

# Usa the el método de ranking binario para ordenar los genes que muestran 
# un patrón espacial. Este paso puede tomar algunos minutos.

ranktest <- binSpect(visium_kidney, 
                     bin_method = "rank",
                     calc_hub = TRUE, 
                     hub_min_int = 5,
                     spatial_network_name = "spatial_network")

spatFeatPlot2D(visium_kidney, 
               expression_values = "scaled",
               feats = ranktest$feats[1:4], 
               cow_n_col = 2, 
               point_size = 1)

##%######################################################%##
#                                                          #
####              Módulos de co-expresión               ####
#                                                          #
##%######################################################%##

# Encuentra los genes con patrones espaciales similares y crea módulos de co-expresión

# Toma los primeros 500 genes con patrón espacial
ext_spatial_genes <- ranktest[1:500,]$feats

# Calcula las distacias entre los genes
spat_cor_netw_DT <- detectSpatialCorFeats(
  visium_kidney,
  method = "network",
  spatial_network_name = "spatial_network",
  subset_feats = ext_spatial_genes)

# Identifica los genes con mayor correlación a un cierto gene
top10_genes <- showSpatialCorFeats(spat_cor_netw_DT, 
                                   feats = "Mbp", 
                                   show_top_feats = 10)

spatFeatPlot2D(visium_kidney, 
               expression_values = "scaled",
               feats = as.character(top10_genes$variable[1:4]), 
               point_size = 1.5)

# Agrupa los genes con mayor correlación, define el número de clusters usando 
# el parámetro k 
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 20)

# Visualiza la correlación de los genes
heatmSpatialCorFeats(visium_kidney,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL))

# Ordena los clusters
netw_ranks <- rankSpatialCorGroups(
  visium_kidney,
  spatCorObject = spat_cor_netw_DT, 
  use_clus_name = "spat_netw_clus")

# Calcula el score de enriquecimiento de los metagenes por cluster
cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT, 
                                        use_clus_name = "spat_netw_clus", 
                                        show_top_feats = 1)

cluster_genes <- cluster_genes_DT$clus 
names(cluster_genes) <- cluster_genes_DT$feat_ID

visium_kidney <- createMetafeats(visium_kidney, 
                                feat_clusters = cluster_genes, 
                                name = "cluster_metagene")

spatCellPlot(visium_kidney,
             spat_enr_names = "cluster_metagene",
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, 
             cow_n_col = 5)

##%######################################################%##
#                                                          #
####              Calcula Leiden clusters               ####
####           usando la información espacial           ####
#                                                          #
##%######################################################%##

# Toma los primeros genes por cada módulo de co-expresión
coexpr_dt <- data.table::data.table(
  genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
  cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)

data.table::setorder(coexpr_dt, cluster)

top30_coexpr_dt <- coexpr_dt[, head(.SD, 30) , by = cluster]

spatial_genes <- top30_coexpr_dt$genes

##%######################################################%##
#                                                          #
####      Re-calcula la dimensión de reducciones y      ####
####      los clusters usando los genes espaciales      ####
#                                                          #
##%######################################################%##

visium_kidney <- runPCA(gobject = visium_kidney,
                       feats_to_use = spatial_genes,
                       name = "custom_pca")

visium_kidney <- runUMAP(visium_kidney, 
                        dim_reduction_name = "custom_pca", 
                        dimensions_to_use = 1:20,
                        name = "custom_umap")

visium_kidney <- createNearestNetwork(gobject = visium_kidney,
                                     dim_reduction_name = "custom_pca",
                                     dimensions_to_use = 1:20, 
                                     k = 5,
                                     name = "custom_NN")

visium_kidney <- doLeidenCluster(gobject = visium_kidney, 
                                network_name = "custom_NN",
                                resolution = 0.15, 
                                n_iterations = 1000,
                                name = "custom_leiden")

spatPlot2D(visium_kidney, 
           cell_color = "custom_leiden",
           point_size = 3)
