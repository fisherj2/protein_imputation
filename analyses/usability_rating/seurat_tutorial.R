library(Seurat)
library(ggplot2)
library(patchwork)
options(SeuratData.repo.use = "http://seurat.nygenome.org")

reference <- load('/scratch/jfisher2/protein_prediction/source_data/Hao_et_al_Seurat_obj.Rdata') %>% get()
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

library(SeuratData)
InstallData('pbmc3k')

pbmc3k <- LoadData('pbmc3k')
pbmc3k <- UpdateSeuratObject(pbmc3k)

pbmc3k <- SCTransform(pbmc3k, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc3k,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

pbmc3k <- MapQuery(
  anchorset = anchors,
  query = pbmc3k,
  reference = reference,
  refdata = list(
    celltype = "celltype",
    detailed_celltype = "detailed_celltype",
    predicted_ADT = "protein"
  ),
  reference.reduction = "pca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(pbmc3k, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc3k, reduction = "ref.umap", group.by = "predicted.detailed_celltype", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

FeaturePlot(pbmc3k, features = c("pDC", "CD16 Mono", "Treg"),  reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))


library(ggplot2)
plot <- FeaturePlot(pbmc3k, features = "CD16 Mono",  reduction = "ref.umap", cols = c("lightgrey", "darkred"))  + ggtitle("CD16 Mono") + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(color = "Prediction Score") +  xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 25)) 
#ggsave(filename = "../output/images/multimodal_reference_mapping.jpg", height = 7, width = 12, plot = plot, quality = 50)

treg_markers <- FindMarkers(pbmc3k, ident.1 = "Treg", only.pos = TRUE, logfc.threshold = 0.1)
print(head(treg_markers))

DefaultAssay(pbmc3k) <- 'predicted_ADT'
# see a list of proteins: rownames(pbmc3k)
FeaturePlot(pbmc3k, features = c("CD3-1", "CD45RA", "IgD"), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 3)

reference <- DietSeurat(reference, counts = FALSE, dimreducs = "spca")
pbmc3k <- DietSeurat(pbmc3k, counts = FALSE, dimreducs = "ref.spca")

#merge reference and query
reference$id <- 'reference'
pbmc3k$id <- 'query'
refquery <- merge(reference, pbmc3k)
refquery[["spca"]] <- merge(reference[["spca"]], pbmc3k[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)


# Both datasets are available through SeuratData
library(SeuratData)
#load reference data
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
#load query data
InstallData('hcabm40k')
hcabm40k <- LoadData(ds = "hcabm40k")

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
DimPlot(bm, group.by = "detailed_celltype", reduction = "wnn.umap") 


bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

library(dplyr)
library(SeuratData)
InstallData('hcabm40k')
hcabm40k.batches <- SplitObject(hcabm40k, split.by = "orig.ident")

hcabm40k.batches <- lapply(X = hcabm40k.batches, FUN = NormalizeData, verbose = FALSE)

anchors <- list()
for (i in 1:length(hcabm40k.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = bm,
    query = hcabm40k.batches[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(hcabm40k.batches)) {
  hcabm40k.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = hcabm40k.batches[[i]],
    reference = bm, 
    refdata = list(
      celltype = "detailed_celltype", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

p1 <- DimPlot(hcabm40k.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3)
p2 <- DimPlot(hcabm40k.batches[[2]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3)
p1 + p2 + plot_layout(guides = "collect")

# Merge the batches 
hcabm40k <- merge(hcabm40k.batches[[1]], hcabm40k.batches[2:length(hcabm40k.batches)], merge.dr = "ref.umap")
DimPlot(hcabm40k, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

p3 <- FeaturePlot(hcabm40k, features = c("rna_TRDC", "rna_MPO", "rna_AVP"), reduction = 'ref.umap', 
                  max.cutoff = 3, ncol = 3)

# cell type prediction scores
DefaultAssay(hcabm40k) <- 'prediction.score.celltype'
p4 <- FeaturePlot(hcabm40k, features = c("CD16 Mono", "HSC", "Prog-RBC"), ncol = 3, 
                  cols = c("lightgrey", "darkred"))

# imputed protein levels
DefaultAssay(hcabm40k) <- 'predicted_ADT'
p5 <- FeaturePlot(hcabm40k, features = c("CD45RA", "CD16", "CD161"), reduction = 'ref.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 3)
p3 / p4 / p5