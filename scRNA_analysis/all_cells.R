library(dplyr)
library(Seurat)

#load data
all.data <- Read10X(data.dir = "path/to/data")
all <- CreateSeuratObject(counts = all.data, project = "all_mouse")

#add gemgroup and names as metadata columns
classification.vec <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec) <- colnames(all@assays$RNA@data)
all <- AddMetaData(all, classification.vec, "gem.group")
head(all@meta.data)

#create a column with actual genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec1) <- colnames(all@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "2IVS_exp2"
classification.vec1[classification.vec1==2] <- "2LV_exp2"
classification.vec1[classification.vec1==3] <- "2LV_exp2"
classification.vec1[classification.vec1==4] <- "2RV_exp2"
classification.vec1[classification.vec1==5] <- "4RV_exp2"
classification.vec1[classification.vec1==6] <- "4LV_exp2"
classification.vec1[classification.vec1==7] <- "7IVS_exp2"
classification.vec1[classification.vec1==8] <- "7LV_exp2"
classification.vec1[classification.vec1==9] <- "7RV_exp2"
classification.vec1[classification.vec1==10] <- "8IVS_exp2"
classification.vec1[classification.vec1==11] <- "8LV_exp2"
classification.vec1[classification.vec1==12] <- "8RV_exp2"
classification.vec1[classification.vec1==13] <- "3IVS_exp1"
classification.vec1[classification.vec1==14] <- "3IVS_exp1"
classification.vec1[classification.vec1==15] <- "3RV_exp1"
classification.vec1[classification.vec1==16] <- "3LV_exp1"
classification.vec1[classification.vec1==17] <- "3LV_exp1"
classification.vec1[classification.vec1==18] <- "4IVS_exp1"
classification.vec1[classification.vec1==19] <- "4IVS_exp1"
classification.vec1[classification.vec1==20] <- "4IVS_exp1"
classification.vec1[classification.vec1==21] <- "4RV_exp1"
classification.vec1[classification.vec1==22] <- "4LV_exp1"
classification.vec1[classification.vec1==23] <- "IVS_BL6"
classification.vec1[classification.vec1==24] <- "IVS_BL6"
classification.vec1[classification.vec1==25] <- "LV_BL6"
classification.vec1[classification.vec1==26] <- "RV_BL6"
all <- AddMetaData(all, classification.vec1, "niceorder")
head(all@meta.data)

#create a column with actual genotype names 
classification.vec2 <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec2) <- colnames(all@assays$RNA@data)
classification.vec2[classification.vec2==1] <- "exp2"
classification.vec2[classification.vec2==2] <- "exp2"
classification.vec2[classification.vec2==3] <- "exp2"
classification.vec2[classification.vec2==4] <- "exp2"
classification.vec2[classification.vec2==5] <- "exp2"
classification.vec2[classification.vec2==6] <- "exp2"
classification.vec2[classification.vec2==7] <- "exp2"
classification.vec2[classification.vec2==8] <- "exp2"
classification.vec2[classification.vec2==9] <- "exp2"
classification.vec2[classification.vec2==10] <- "exp2"
classification.vec2[classification.vec2==11] <- "exp2"
classification.vec2[classification.vec2==12] <- "exp2"
classification.vec2[classification.vec2==13] <- "exp1"
classification.vec2[classification.vec2==14] <- "exp1"
classification.vec2[classification.vec2==15] <- "exp1"
classification.vec2[classification.vec2==16] <- "exp1"
classification.vec2[classification.vec2==17] <- "exp1"
classification.vec2[classification.vec2==18] <- "exp1"
classification.vec2[classification.vec2==19] <- "exp1"
classification.vec2[classification.vec2==20] <- "exp1"
classification.vec2[classification.vec2==21] <- "exp1"
classification.vec2[classification.vec2==22] <- "exp1"
classification.vec2[classification.vec2==23] <- "BL6"
classification.vec2[classification.vec2==24] <- "BL6"
classification.vec2[classification.vec2==25] <- "BL6"
classification.vec2[classification.vec2==26] <- "BL6"
all <- AddMetaData(all, classification.vec2, "batch")
head(all@meta.data)

#create a column with actual genotype names 
classification.vec3 <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec3) <- colnames(all@assays$RNA@data)
classification.vec3[classification.vec3==1] <- "hetIVS"
classification.vec3[classification.vec3==2] <- "hetLV"
classification.vec3[classification.vec3==3] <- "hetLV"
classification.vec3[classification.vec3==4] <- "hetRV"
classification.vec3[classification.vec3==5] <- "hetRV"
classification.vec3[classification.vec3==6] <- "hetLV"
classification.vec3[classification.vec3==7] <- "controlIVS"
classification.vec3[classification.vec3==8] <- "controlLV"
classification.vec3[classification.vec3==9] <- "controlRV"
classification.vec3[classification.vec3==10] <- "controlIVS"
classification.vec3[classification.vec3==11] <- "controlLV"
classification.vec3[classification.vec3==12] <- "controlRV"
classification.vec3[classification.vec3==13] <- "controlIVS"
classification.vec3[classification.vec3==14] <- "controlIVS"
classification.vec3[classification.vec3==15] <- "controlRV"
classification.vec3[classification.vec3==16] <- "controlLV"
classification.vec3[classification.vec3==17] <- "controlLV"
classification.vec3[classification.vec3==18] <- "hetIVS"
classification.vec3[classification.vec3==19] <- "hetIVS"
classification.vec3[classification.vec3==20] <- "hetIVS"
classification.vec3[classification.vec3==21] <- "hetRV"
classification.vec3[classification.vec3==22] <- "hetLV"
classification.vec3[classification.vec3==23] <- "IVS_BL6"
classification.vec3[classification.vec3==24] <- "IVS_BL6"
classification.vec3[classification.vec3==25] <- "LV_BL6"
classification.vec3[classification.vec3==26] <- "RV_BL6"
all <- AddMetaData(all, classification.vec3, "hetvscontrol")
head(all@meta.data)

#create a column with anatomical region names 
classification.vec4 <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec4) <- colnames(all@assays$RNA@data)
classification.vec4[classification.vec4==1] <- "IVS"
classification.vec4[classification.vec4==2] <- "LV"
classification.vec4[classification.vec4==3] <- "LV"
classification.vec4[classification.vec4==4] <- "RV"
classification.vec4[classification.vec4==5] <- "RV"
classification.vec4[classification.vec4==6] <- "LV"
classification.vec4[classification.vec4==7] <- "IVS"
classification.vec4[classification.vec4==8] <- "LV"
classification.vec4[classification.vec4==9] <- "RV"
classification.vec4[classification.vec4==10] <- "IVS"
classification.vec4[classification.vec4==11] <- "LV"
classification.vec4[classification.vec4==12] <- "RV"
classification.vec4[classification.vec4==13] <- "IVS"
classification.vec4[classification.vec4==14] <- "IVS"
classification.vec4[classification.vec4==15] <- "RV"
classification.vec4[classification.vec4==16] <- "LV"
classification.vec4[classification.vec4==17] <- "LV"
classification.vec4[classification.vec4==18] <- "IVS"
classification.vec4[classification.vec4==19] <- "IVS"
classification.vec4[classification.vec4==20] <- "IVS"
classification.vec4[classification.vec4==21] <- "RV"
classification.vec4[classification.vec4==22] <- "LV"
classification.vec4[classification.vec4==23] <- "IVS"
classification.vec4[classification.vec4==24] <- "IVS"
classification.vec4[classification.vec4==25] <- "LV"
classification.vec4[classification.vec4==26] <- "RV"
all <- AddMetaData(all, classification.vec4, "region")
head(all@meta.data)

#create a column with actual genotype names 
classification.vec5 <- as.numeric(gsub(".*-","", (colnames(all@assays$RNA@data))))
names(classification.vec5) <- colnames(all@assays$RNA@data)
classification.vec5[classification.vec5==1] <- "het"
classification.vec5[classification.vec5==2] <- "het"
classification.vec5[classification.vec5==3] <- "het"
classification.vec5[classification.vec5==4] <- "het"
classification.vec5[classification.vec5==5] <- "het"
classification.vec5[classification.vec5==6] <- "het"
classification.vec5[classification.vec5==7] <- "control"
classification.vec5[classification.vec5==8] <- "control"
classification.vec5[classification.vec5==9] <- "control"
classification.vec5[classification.vec5==10] <- "control"
classification.vec5[classification.vec5==11] <- "control"
classification.vec5[classification.vec5==12] <- "control"
classification.vec5[classification.vec5==13] <- "control"
classification.vec5[classification.vec5==14] <- "control"
classification.vec5[classification.vec5==15] <- "control"
classification.vec5[classification.vec5==16] <- "control"
classification.vec5[classification.vec5==17] <- "control"
classification.vec5[classification.vec5==18] <- "het"
classification.vec5[classification.vec5==19] <- "het"
classification.vec5[classification.vec5==20] <- "het"
classification.vec5[classification.vec5==21] <- "het"
classification.vec5[classification.vec5==22] <- "het"
classification.vec5[classification.vec5==23] <- "control"
classification.vec5[classification.vec5==24] <- "control"
classification.vec5[classification.vec5==25] <- "control"
classification.vec5[classification.vec5==26] <- "control"
all <- AddMetaData(all, classification.vec5, "genotype")
head(all@meta.data)

#QC
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#filter low quality cells
all <- subset(all, subset = nFeature_RNA > 1500 & nFeature_RNA < 7250 & nCount_RNA > 10000 & nCount_RNA < 50000)

#save the mega object
save(all, file = "allmouseregions_Seurat3.Robj")

#normalize data 
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)

#identify highly variable genes
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(all), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(all)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#cellcycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015 is read and segregated this list into markers of G2/M phase and markers of S phase
cc.genes <- readLines(con = "~/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Assign cell cycle scores
all <- CellCycleScoring(object = all, s.features = s.genes, g2m.features = g2m.genes,
                        set.ident = TRUE)
head(x = all@meta.data)
all@meta.data$CC.Difference <- all@meta.data$S.Score - all@meta.data$G2M.Score

#scale data
all.genes <- rownames(all)
#regress unwanted variables
all <- ScaleData(all, vars.to.regress = c("percent.mt", "CC.Difference"))

#dimensionality reduction
all <- RunPCA(all, features = VariableFeatures(object = all))

# Examine and visualize PCA results a few different ways
print(all[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all, dims = 1:2, reduction = "pca")
DimPlot(all, reduction = "pca")
ElbowPlot(all, ndims = 50)

#cluster cells
all <- FindNeighbors(all, dims = 1:40)
all <- FindClusters(all, resolution = c(0.5,0.6,0.7,0.8,0.9,1.0))

#umap
all <- RunUMAP(all, dims = 1:40)
Idents(all) <- "RNA_snn_res.0.5"
DimPlot(all, reduction = "umap", label = T, label.size = 6) + NoLegend()
DimPlot(object = all, reduction = "umap", group.by = "region", pt.size = 0.4, cols = c("magenta3","green3","navy"), shuffle = T) + NoLegend()
DimPlot(object = all, reduction = "umap", group.by = "genotype",cols = c("steelblue3","red3"), shuffle=T, pt.size = 0.4) + NoLegend()
DimPlot(object = all, reduction = "umap", group.by = "batch",pt.size = 0.4,cols = c("orange","blue","green4"), shuffle = T) + NoLegend()


#feature plot #endocardial- Nfatc1, Tie2
#Conduction system - Gja1, Hcn4
#cell-cell connection - Gja3
#epicardial cells - Wt1, Tbx18
#Nppa, Pln, Ryr2, Tecrl, Trh, Rspo3, Nkx2-5 - Tbx5 dependent genes
FeaturePlot(all, features = c("Tnnt2","Postn","Wt1","Plvap","Hbb-b2","C1qb"))


##rename clusters by celltype and add metadata column
Idents(all) <- "RNA_snn_res.0.5"
newclusterids <- c("Endothelial","CM","CM","Fibroblast","CM","Epicardial","Fibroblast","Endothelial","Endothelial","CM","Endothelial",
                   "CM","CM","Complement","Blood","Endothelial","CM")
names(newclusterids) <- levels(all)
all <- RenameIdents(all, newclusterids)
all$Celltype <- Idents(all)
DimPlot(object = all, reduction = "umap", group.by = "Celltype",pt.size = 0.4, label =T, label.size = 9) #+ NoLegend()
save(all, file = "allmouseregions_Seurat3_processed.Robj")


#save processed object
save(all, file = "allmouseregions_Seurat3_processed.Robj")
