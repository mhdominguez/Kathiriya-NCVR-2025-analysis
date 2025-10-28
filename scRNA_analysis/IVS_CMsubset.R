#subset cm clusters only from IVSall
#resolution set to 0.5
IVScms <- subset(IVSall, idents = c(0,2,5,6,10,14))

#create a column with anatomical region names 
classification.vec4 <- as.numeric(gsub(".*-","", (colnames(IVScms@assays$RNA@data))))
names(classification.vec4) <- colnames(IVScms@assays$RNA@data)
classification.vec4[classification.vec4==1] <- "IVS"
classification.vec4[classification.vec4==7] <- "IVS"
classification.vec4[classification.vec4==10] <- "IVS"
classification.vec4[classification.vec4==13] <- "IVS"
classification.vec4[classification.vec4==14] <- "IVS"
classification.vec4[classification.vec4==18] <- "IVS"
classification.vec4[classification.vec4==19] <- "IVS"
classification.vec4[classification.vec4==20] <- "IVS"
classification.vec4[classification.vec4==23] <- "IVS"
classification.vec4[classification.vec4==24] <- "IVS"
IVScms <- AddMetaData(IVScms, classification.vec4, "region")
head(IVScms@meta.data)

#create a column with actual genotype names 
classification.vec5 <- as.numeric(gsub(".*-","", (colnames(IVScms@assays$RNA@data))))
names(classification.vec5) <- colnames(IVScms@assays$RNA@data)
classification.vec5[classification.vec5==1] <- "mutant"
classification.vec5[classification.vec5==7] <- "control"
classification.vec5[classification.vec5==10] <- "control"
classification.vec5[classification.vec5==13] <- "control"
classification.vec5[classification.vec5==14] <- "control"
classification.vec5[classification.vec5==18] <- "mutant"
classification.vec5[classification.vec5==19] <- "mutant"
classification.vec5[classification.vec5==20] <- "mutant"
classification.vec5[classification.vec5==23] <- "control"
classification.vec5[classification.vec5==24] <- "control"
IVScms <- AddMetaData(IVScms, classification.vec5, "genotype")
head(IVScms@meta.data)

#create a column with actual genotype names 
classification.vec2 <- as.numeric(gsub(".*-","", (colnames(IVScms@assays$RNA@data))))
names(classification.vec2) <- colnames(IVScms@assays$RNA@data)
classification.vec2[classification.vec2==1] <- "exp2"
classification.vec2[classification.vec2==7] <- "exp2"
classification.vec2[classification.vec2==10] <- "exp2"
classification.vec2[classification.vec2==13] <- "exp1"
classification.vec2[classification.vec2==14] <- "exp1"
classification.vec2[classification.vec2==18] <- "exp1"
classification.vec2[classification.vec2==19] <- "exp1"
classification.vec2[classification.vec2==20] <- "exp1"
classification.vec2[classification.vec2==23] <- "exp3"
classification.vec2[classification.vec2==24] <- "exp3"
IVScms <- AddMetaData(IVScms, classification.vec2, "batch")
head(IVScms@meta.data)

#create a column with actual genotype names 
classification.vec3 <- as.numeric(gsub(".*-","", (colnames(IVScms@assays$RNA@data))))
names(classification.vec3) <- colnames(IVScms@assays$RNA@data)
classification.vec3[classification.vec3==1] <- "mutant_IVS"
classification.vec3[classification.vec3==7] <- "control_IVS"
classification.vec3[classification.vec3==10] <- "control_IVS"
classification.vec3[classification.vec3==13] <- "control_IVS"
classification.vec3[classification.vec3==14] <- "control_IVS"
classification.vec3[classification.vec3==18] <- "mutant_IVS"
classification.vec3[classification.vec3==19] <- "mutant_IVS"
classification.vec3[classification.vec3==20] <- "mutant_IVS"
classification.vec3[classification.vec3==23] <- "control_IVS"
classification.vec3[classification.vec3==24] <- "control_IVS"
IVScms <- AddMetaData(IVScms, classification.vec3, "genotype_region")
head(IVScms@meta.data)

#create a column with actual genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(IVScms@assays$RNA@data))))
names(classification.vec1) <- colnames(IVScms@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "mutant_IVS_replicate2"
classification.vec1[classification.vec1==7] <- "control_IVS_replicate2"
classification.vec1[classification.vec1==10] <- "control_IVS_replicate3"
classification.vec1[classification.vec1==13] <- "control_IVS_replicate1_techrep1"
classification.vec1[classification.vec1==14] <- "control_IVS_replicate1_techrep2"
classification.vec1[classification.vec1==18] <- "mutant_IVS_replicate1,techrep1"
classification.vec1[classification.vec1==19] <- "mutant_IVS_replicate1,techrep2"
classification.vec1[classification.vec1==20] <- "mutant_IVS_replicate1,techrep3"
classification.vec1[classification.vec1==23] <- "control_IVS_replicate4_techrep1"
classification.vec1[classification.vec1==24] <- "control_IVS_replicate4_techrep2"
IVScms <- AddMetaData(IVScms, classification.vec1, "genotype_region_replicate_techrep")
head(IVScms@meta.data)

# #normalize data 
IVScms <- NormalizeData(IVScms, normalization.method = "LogNormalize", scale.factor = 10000)
# #identify highly variable genes
IVScms <- FindVariableFeatures(IVScms, selection.method = "vst", nfeatures = 1500)
# # Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(IVScms), 20)
# # plot variable features with and without labels
plot1 <- VariableFeaturePlot(IVScms)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
# #scale data
all.genes <- rownames(IVScms)
# #regress unwanted variables
IVScms <- ScaleData(IVScms, vars.to.regress = c("percent.mt","CC.Difference"))

#dimensionality reduction
IVScms <- RunPCA(IVScms)

# Examine and visualize PCA results a few different ways
print(IVScms[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(IVScms, dims = 1:2, reduction = "pca")
DimPlot(IVScms, reduction = "pca")
DimHeatmap(IVScms, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(IVScms, ndims = 40)

#cluster cells using future to parallilize FindClusters and running multiple resolutions
library(future)
plan()
plan("multiprocess", workers = 4)
IVScms <- FindNeighbors(IVScms, dims = 1:30)
IVScms <- FindClusters(IVScms, resolution = c(0.1,0.2,0.3,0.4))
IVScms <- FindClusters(IVScms, resolution = c(0.5,0.6,0.7,0.8,0.9,1.0))
IVScms <- FindClusters(IVScms, resolution = c(1.1,1.2,1.3,1.4,1.5))
IVScms <- FindClusters(IVScms, resolution = c(1.6,1.7,1.8,1.9,2.0))

#umap
IVScms <- RunUMAP(IVScms, dims = 1:30)
Idents(IVScms) <- "SCT_snn_res.0.8"
DimPlot(IVScms, reduction = "umap", label =T, label.size=6, pt.size = 0.8) + NoLegend()
DimPlot(object = IVScms, group.by = "genotype", cols = c("blue","red"), pt.size = 0.6) + NoLegend()
DimPlot(object = IVScms, reduction = "umap", group.by = "batch")
FeaturePlot(IVScms, features = c("Rspo3","Bmp2","Robo1","Unc5b","Cited1","Nppa","Slit2","Ntn1"))


#Cluster Tree
IVScms <- BuildClusterTree(IVScms)
PlotClusterTree(IVScms)


##pseudobulk in ivs cms
Idents(IVScms) <- "genotype"
IVScms_contvsmutant_pseudobulk <- FindMarkers(IVScms, ident.1= "control", ident.2= "mutant")
write.csv(IVScms_contvsmutant_pseudobulk, file = "IVScms_contvsmutant_pseudobulk.csv")
##make a heatmap for pseudobulk
genes <- IVScms_contvsmutant_pseudobulk[['x']]
p <- DoHeatmap(IVScms, features = genes, group.by = "ident") + scale_x_discrete(limits=c('control','het')) + theme(axis.text.y = element_text(hjust = 1,size =12, face = 'bold')) + scale_fill_gradientn(colors = c("blue", "white", "green"))


##cluster-cluster comparison
cluster10vs5_res0.8 <- FindMarkers(IVScms, ident.1 = 10, ident.2 = 5)
write.csv(cluster10vs5_res0.8, file = "cluster10vs5_res0.8.csv")


##dotplot -for manuscript Fig 5H
library(ggplot2)
Idents(IVScms) <- "SCT_snn_res.0.8"
cluster5and10 <- subset(IVScms, idents = c(5,10))
Idents(cluster5and10, cells=WhichCells(cluster5and10, idents = 5)) <- 'cluster5'
Idents(cluster5and10, cells=WhichCells(cluster5and10, idents = 10)) <- 'cluster10'
genes <- c("Ccnd2","Tbx5","Nkx2-5","Gata6","Unc5b","Robo1","Tbx20",'Tbx3','Myh6','Sfrp1','Id3','Pitx2','Fhl2','Acta1','Id2','Rspo3',
          "Irx3","Vcan","Hand1","Smyd2","Hopx","Angpt1","Cav3","Gja5","Ckm",'Pln','Smpx',"Myl2",'Myh7','Myl3','Gja1',"Ntn1","Slit2",'Nppb','Nppa','Cited1')
genes2 <- c("Igfbp5","Rspo3","Bmp2","Bambi","Sln","Id1","Id2","Id3","Tgfb2","Unc5b","Robo1","Tnnt1","Myl7","Myl1","Myh6","Tbx3","Tbx20","Gata6","Nkx2-5",
            "Tbx5","Cdh2","Bmp4","Cited2","Hand2")
genes3 <- c("Cited1","Nppa","Nppb","Gja1","Gja5","Slit2","Ntn1","Myl3","Myh7","Myl2","Smpx","Pln","Gyg","Hopx","Smyd2","Hand1","Irx3","Cxcl12","Ttn","Tnnt2","Mef2c",
            "Thbs4","Actc1","Hey2")
DotPlot(cluster5and10, features = genes3, dot.min = 0) + scale_y_discrete(expand=c(6,0), limits=c('cluster5','cluster10')) + scale_x_discrete(expand=c(0.2,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,size =12, face="italic")) + NoLegend()


##subset tdTomato+ by expression, control vs mutant markers
IVSred <- subset(IVScms, subset = `tdTomato-full` > 0.001)
Idents(IVSred) <- "genotype"
# IVScms_red_mutantvscont <- FindMarkers(IVSred, ident.1 = "mutant", ident.2 = "control")
# write.csv(IVScms_red_mutantvscont, file = "IVScms_red_mutantvscontrol.csv")
##dotplot
genes <- c("Cited1","Nppa","Nppb","Gja1","Slit2","Ntn1","Gja5","Pln","Unc5b","Id3","Id2","Myl7","Tagln","Igfbp5","Bmp2","Bambi","Rspo3")
DotPlot(IVSred, features = genes, dot.min = 0) + scale_y_discrete(expand=c(6,0), limits=c('control','mutant'))+ scale_x_discrete(expand = c(0.5,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,size =12, face="italic"))


##subset Zsgreen+ by expression, control vs mutant markers
IVSgreen <- subset(IVScms, subset = `ZsGreen-full` > 0.001)
Idents(IVSgreen) <- "genotype"
# IVScms_green_mutantvscont <- FindMarkers(IVSgreen, ident.1 = "mutant", ident.2 = "control")
# write.csv(IVScms_green_mutantvscont, file = "IVScms_green_mutantvscontrol.csv")
##dotplot
genes <- c("Cited1","Nppa","Nppb","Gja1","Slit2","Mest","Kcnk3","Casq1","Cxcl12","Id3","Id2","Myl7","Tagln","Igfbp5","Bmp2","Bambi","Rspo3")
DotPlot(IVSgreen, features = genes, dot.min = 0) + scale_y_discrete(expand=c(6,0), limits=c('control','mutant'))+ scale_x_discrete(expand = c(0.5,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,size =12, face="italic"))


#save object
save(IVScms, file = "IVScmsreclustered_stdpipeline_Seurat3.Robj")
