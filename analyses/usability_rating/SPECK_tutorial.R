library(SPECK)
library(Seurat)
library(ggplot2)
library(gridExtra)

data("pbmc.rna.mat")
dim(pbmc.rna.mat)


speck.full <- speck(counts.matrix = pbmc.rna.mat, rank.range.end = 100,
                    min.consec.diff = 0.01, rep.consec.diff = 2,
                    manual.rank = NULL, max.num.clusters = 4,
                    seed.rsvd = 1, seed.ckmeans = 2)
speck.rank <- speck.full$rrr.rank
paste("Rank: ", speck.rank, sep = "")
#> [1] "Rank: 23"
plot(speck.full$component.stdev, ylab = "Stdev. of non-centered sample PCs",
     xlab = "Rank range", main = paste("Selected rank (k=", speck.rank, ")", sep=""))
abline(v = speck.rank, lty = 2, col = "red")

head(speck.full$clust.num); table(speck.full$clust.num)
#> number.clusters
#> MIR1302-2HG 4
#> FAM138A 1
#> OR4F5 2
#> AL627309.1 1
#> AL627309.3 1
#> AL627309.2 3
#>
#> 1 2 3 4
#> 25658 5883 1470 527
head(speck.full$clust.max.prop)
#> proportion.max.clust
#> MIR1302-2HG 6.6
#> FAM138A 0.0
#> OR4F5 0.0
#> AL627309.1 0.0
#> AL627309.3 0.0
#> AL627309.2 0.0
speck.output <- speck.full$thresholded.mat
paste("# of samples in RRR object:", dim(speck.output)[1])
#> [1] "# of samples in RRR object: 1000"
paste("# of genes in RRR object:", dim(speck.output)[2])
#> [1] "# of genes in RRR object: 33538"
SPECK_assay <- CreateAssayObject(counts = t(speck.output))
pbmc.rna.seurat <- CreateSeuratObject(counts = t(as.matrix(pbmc.rna.mat)))
pbmc.rna.seurat[["SPECK"]] <- SPECK_assay


