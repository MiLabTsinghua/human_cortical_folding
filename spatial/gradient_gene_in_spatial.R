library(Seurat)
library(AUCell)

mat <- mat[mat$flag == 'Pass',]
genelist <- split(mat$gene, mat$cluster)

sob_sp <- get_spatial_sob(tp = tp)
s <- subset(sob_sp, unsup_layer_cl2r %in% layerlist)
if ((tp == 'GW13') & (tag == 'EX')) {
    s <- subset(s, ang <=5)
}
mat <- as.matrix(s@assays$SCT@data)

cells_rankings <- AUCell_buildRankings(mat, nCores = 1, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(genelist, cells_rankings)
t <- getAUC(cells_AUC)
for (i in 1:nrow(t)) {
    s@meta.data[,rownames(t)[i]] <- t[i,colnames(s)]
}
