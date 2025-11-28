############################ scRNA-Seq Data Analysis #######################

############## Data is from 10X genomics platform ###################
########### 1k Brain cells from E18 mouse ################

### https://www.10xgenomics.com/datasets/1-k-brain-cells-from-an-e-18-mouse-2-standard-2-1-0 ####

###### Cells are from hippocampus, cortex & sub-ventricular zone #############

############## Seurat package is used for analysis #################


################ LET'S START THE ANALYSIS ########################

########### Import required libraries & packages #############################

library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)



###################### STEP-1: LOAD THE DATA & CREATE SEURAT OBJECT ############

######## Download & extract the raw counts matrix files from 10x genomics


counts <- Read10X(data.dir = "Data") #### Read in the data 

############### Create seurat object ####################

############ Pre-filtering ################

######### Selecting genes that are expressed in at least 3 cells #########

########## Excluding cells expressing less than 200 genes #############

seurat <- CreateSeuratObject(counts = counts, min.cells = 3,
                             min.features = 200,
                             project = "Mouse_Brain") 

str(seurat)

seurat #### 14393 features across 1876 samples 

############################ STEP-2: QUALITY CONTROL #####################################

View(seurat@meta.data) #### Contains cells & features 



################# PURPOSE OF QC #########################

###### Filter out cells with low & high number of genes (features)

##### Low gene count means cells have not been properly sequenced

##### High gene count represent doublets or mutiplets (two or multiple cells seq togather)

###### Filter out cells with high mitochondrial genes (dying cells)

grep("^mt", row.names(seurat)) ### Check mitochondrial genes

######### Calculate %age of mitochondrial genes ##########################

seurat[["percent.mt"]]<- PercentageFeatureSet(seurat, pattern = "^mt")

View(seurat@meta.data) #### %mit genes added


############ Plot the data before filtering ###################

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#### Looking at plot, approx 7500 genes & 70000 transcripts are in our data. 


############## Filtering the data on default parameters ###################

###### Selecting cells that have gene count between 200 to 2500 and 
###### mitochondrial genes less than 5%


seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                   percent.mt < 5)

seurat #### 14393 features across 834 samples 

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

################################ STEP-3: NORMALIZATION #########################

######### PURPOSE #################

########## To make the cells comparable by accounting for technical differences

############ Defeault method is Log Normalization ####################

############ Calculate gene expression across all samples
######### Divide gene expression by total expression
######### Multiply it by scaling factor (default=10,000)

seurat <- NormalizeData(seurat)


########################### STEP-4: IDENTIFY HIGHLY VARIABLE FEATURES #############

########### PURPOSE ##############

###### Focus our analysis on genes that actually differ between cells 
###### Ignore the technical noise, batch effects and even biological variation 
###### biological sources of variation include cell cycle stage

###### Seurat calculates the average expression & dispersion for each gene
###### Places genes into bins
###### Then calculates z-scores for genes

seurat <- FindVariableFeatures(seurat, nfeatures = 2000) ### Selecting default 2000 variable features

######## Identify top10 higly variable features 

top10 <- head(VariableFeatures(seurat), 10)

######### Plot the variable features 

plot1 <- VariableFeaturePlot(seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


########################### STEP-5: SCALE THE DATA ###########################

####### PURPOSE
##### Transform the data so that each gene contributes equally to the analysis
##### Remove unwanted sources of variation by using var.to.regress function

seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))


########################## STEP-6: LINEAR DIMENSIONALITY REDUCTION #############

######### PURPOSE 
###### Makes the data compact for downstream analysis
##### Principal component analysis (PCA): which features define the dataset

seurat <- RunPCA(seurat, npcs = 50)

########### Plot the PCs using elbow plot 

####### Elbow plot show how much variance each PC explains
###### Elbow plot can be used to narrow down how many PCs to use in downstream analysis

ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

##### After 20 PCs, the PCs don't contribute much towards explaining the data

#### Hence, we'll select first 20 PCs


##################### STEP-7: NON-LINEAR DIMENSIONALITY REDUCTION ##############

######## PURPOSE 
###### Visualize data (PCs) in high dimensional space 
##### Methods: UMAP or t-SNE

seurat <- RunUMAP(seurat, dims = 1:20) #### UMAP on first 20 PCs

seurat <- RunTSNE(seurat, dims = 1:20) #### t-SNE on 20PCs

######### Plot UMAP & t-SNE 

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2


############## Check whether certain cell markers exist in dataset or not

plot1 <- FeaturePlot(seurat, c("Olig1", "Olig2", "Mbp", "Cd11b", "Gad1", "Gad2", "Gja1", "Cldn5", "Sox10", "Sox2", "Cd68"),
                     ncol=3, reduction = "tsne")

plot2 <- FeaturePlot(seurat, c("Olig1", "Olig2", "Mbp", "Cd11b", "Gad1", "Gad2", "Gja1", "Cldn5", "Sox10", "Sox2", "Cd68"),
                     ncol=3, reduction = "umap")

plot1 / plot2

########################## STEP-8: CLUSTER THE CELLS ###########################

###### PURPOSE
###### Cluster the cells together wit similar expression patterns
###### helps us understand cellular heterogeneity

seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = c(0.1,0.3, 0.5, 0.7, 1))

View(seurat@meta.data)

########## Visualize cell clusters at different resolution

DimPlot(seurat, group.by = "RNA_snn_res.0.1", label = TRUE) #### 4 clusters

DimPlot(seurat, group.by = "RNA_snn_res.0.3", label = TRUE) ##### 7 clusters

DimPlot(seurat, group.by = "RNA_snn_res.0.5", label = TRUE) #### 7 clusters

DimPlot(seurat, group.by = "RNA_snn_res.0.7", label = TRUE) ##### 8 clusters

DimPlot(seurat, group.by = "RNA_snn_res.1", label = TRUE) #### 9 clusters

###### Selecting 0.7 resolution, the clusters are distinct. 

Idents(seurat)  <- "RNA_snn_res.0.5"

########### Visualize clustering results in UMAP & t-SNE

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2


############################# STEP-9: ANNOTATE CELL CLUSTERS ###################

########## PURPOSE
######### Identify which cell types/states are present in each cluster
###### Multiple ways to annotate clusters

####### For this analysis, we'll identify clusters by alreday defined markers in cellmarker2.0

##################### Identify cluster markers for each cluster ###############

####### Identify regulated gene markers in clusters with log2FC as 0.25

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))

View(cl_markers)

cl_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
############## Visualize top cluster markers in heatmap ####################

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()

############ Save the cluster markers in CSV file ######################

write.csv(cl_markers, file = "cluster_markers.csv")


################ Let's plot individual markers in more detail #################


plot1 <- VlnPlot(seurat, features = c("Pcp4","Mcmbp", "Sox2", "Neurod6", 
                                      "Reln", "Gad1", "Gad2", "Slc32a1", 
                                      "Snap25", "Dcx", "Apoe", "Neurod1",
                                      "Tubb3", "Rbfox3", "Prox1", "Ascl1", 
                                      "Eomes", "Neurog2", "Pax6", "Vim", 
                                      "Cd9", "Foxg1"), pt.size = 0)

plot1

################## Prominent Markers from Each Cluster #########################
# Cluster0: APOE, Pax6, Vim, Cd9
# Cluster1: Neurod6, Tubb3, Vim
# Cluster2: Mcmbp, Neurod6, Snap25, Dcx, Tubb3, Foxg1
# Cluster3: Pcp4, Neurod6, Dcx, Neurod1, Tubb3, Rbfox3, Foxg1, Eomes, Neurog2, Pax6, Vim
# Cluster4: Sox2, gad1, gad2, Slc32a1, Snap25, Dcx, Tubb3, Foxg1
# Cluster5: Sox2, Reln, gad1, gad2, Dcx, Tubb3, Rbfox3, Prox1, Ascl1, Foxg1
# Cluster6: Tubb3, Vim 


####### Look up these markers in literature or databases to determine the cells they are most expressed in ################

##### For this project we used CellMarker2.0 & Pangloadb to determine the cells ##########3

############## Annotating the cells now #########################

####### Based on marker expression, following annotations can be given to clusters #####

# Cluster0: Neural Progenitor Cells
# Cluster1: Immature excitatory neurons
# Cluster2: Proliferating immature excitatory neurons
# Cluster3: Intermediate progenitors
# Cluster4: Immature inhibitory neurons
# Cluster5: Cajal-Retzius Cells
# Cluster6: Radial glial cells

new_ident <- setNames(c("Neural Progenitor Cells",
                        "Immature excitatory neurons",
                        "Proliferating immature excitatory neurons", 
                        "Intermediate progenitors",
                        "Immature inhibitory neurons", 
                        "Cajal-Retzius Cells",
                        "Radial glial cells"),
                      levels(seurat))

seurat <- RenameIdents(seurat, new_ident)

DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

############## Limitations: Used only some markers for cluster annotation #######

