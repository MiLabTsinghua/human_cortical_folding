suppressMessages(library(Seurat))
suppressMessages(library(monocle))

s2 <- subset(sob, anno2tp_fin %in% ctlist)
genelist <- rownames(s2)
genelist <- setdiff(genelist, c(genelist[grepl(pattern = '^A[CLP][0-9]', genelist)], genelist[grepl(pattern = '^LINC', genelist)]))
s2 <- subset(s2, features = genelist)
hvglist <- lapply(unique(s2$batch), function(b) {
    s2 <- subset(s2, batch == b)
    DefaultAssay(s2) <- 'RNA'
    s2 <- FindVariableFeatures(s2, nfeatures = 5000, assay = 'RNA')
    return(VariableFeatures(s2))
})
s2 <- subset(s2, features =  unique(unlist(hvglist)))
s2@meta.data$pseudoage <- pseudoage_folding[colnames(s2),'PC1']
s2 <- subset(s2, pseudoage != 'NA')

x_monocle <- Seurat::as.CellDataSet(s2, assay = 'SCT', reduction = 'umap')
x_monocle$Pseudotime <- s2@meta.data[colnames(x_monocle), 'pseudoage']
x_monocle <- estimateSizeFactors(x_monocle)
x_monocle <- estimateDispersions(x_monocle)
x_monocle$CellType <- as.character(s2@meta.data[colnames(x_monocle),'anno2tp_fin'])
x_monocle$State <- as.character(s2@meta.data[colnames(x_monocle),'location'])
DE <- differentialGeneTest(x_monocle, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1, verbose = FALSE)
cds_subset <- x_monocle[intersect(rownames(x_monocle), rownames(DE)[DE$qval < 0.05 & DE$status == 'OK']), ]
