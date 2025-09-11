library(Seurat)

### normalize
tag = 'GW11'
s <- Load10X_Spatial(paste0('/home/baiyq/projects/fold_dev/data/10X_Visium/', tag, '/outs'))
# normalize1
s <- PercentageFeatureSet(s, pattern = "^MT-", col.name = "percent.mt")
s <- subset(s, subset = (nFeature_Spatial > 800) & (percent.mt < 5))
s <- SCTransform(s, assay = 'Spatial', new.assay.name = 'SCT', vars.to.regress = c('percent.mt'), verbose = FALSE)
# normalize2
s <- CellCycleScoring(s, assay = 'SCT', s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
s <- SCTransform(s, assay = 'Spatial', new.assay.name = 'SCT', vars.to.regress = c('percent.mt','S.Score','G2M.Score'), verbose = FALSE)

### format
get_spatial_sob <- function(s, tp) {
    ### prepare sob
    sob_sp <- s
    sob_sp$unsup_layer_cl2r <- unlist(lapply(sob_sp$unsup_layer_cl2, function(item) strsplit(item, '[_]')[[1]][1]))
    # filter sob
    s <- sob_sp
    if (tp == 'GW11') {
        s <- subset(s, cells = colnames(s)[s@images$slice1@coordinates$col < 100])
        anglist <- lapply(1: nrow(s@meta.data), function(i) {
            x = max(s@images$slice1@coordinates$col)
            y = 37
            A = c(x,y)
            B = c(x,min(s@images$slice1@coordinates$row))
            C = c(s@images$slice1@coordinates$col[i], s@images$slice1@coordinates$row[i])
            D = dist(rbind(rbind(A,B),C))
            D = as.matrix(D)
            a = acos((D['A','B']^2 + D['C','A']^2 - D['C','B']^2) / (2 * D['A','B'] * D['C','A']))
            return(a)
        })
        s@meta.data$ang <- unlist(anglist)
    }
    if (tp == 'GW13') {
        s <- subset(s, cells = colnames(s)[s@images$slice1@coordinates$col < 90])
        anglist <- lapply(1: nrow(s@meta.data), function(i) {
            x = max(s@images$slice1@coordinates$col)
            y = 50
            A = c(x,y)
            B = c(x,min(s@images$slice1@coordinates$row))
            C = c(s@images$slice1@coordinates$col[i], s@images$slice1@coordinates$row[i])
            D = dist(rbind(rbind(A,B),C))
            D = as.matrix(D)
            a = acos((D['A','B']^2 + D['C','A']^2 - D['C','B']^2) / (2 * D['A','B'] * D['C','A']))
            return(a)
        })
        s@meta.data$ang <- unlist(anglist)
    }
    s <- subset(s, subset = unsup_layer_cl2r != 'NA')

    s$ang_ori <- s$ang
    t <- data.frame(y = c(1,6), ang = c(min(s$ang_ori), max(s$ang_ori)))
    model <- glm(y~ang, data = t)
    s$ang <- predict(model, data.frame(ang = s$ang_ori))

    ### renorm
    Folding_dat <- s
    DefaultAssay(Folding_dat) = "Spatial"
    Folding_dat = NormalizeData(Folding_dat, assay = 'Spatial')
    Folding_dat[["percent.mt"]] <- PercentageFeatureSet(Folding_dat, pattern = "^MT-", assay = 'Spatial')
    Folding_dat = CellCycleScoring(Folding_dat,
                                    s.features = cc.genes.updated.2019$s.genes,
                                    g2m.features = cc.genes.updated.2019$g2m.genes,
                                    set.ident = FALSE, assay = 'Spatial')
    Folding_dat = SCTransform(Folding_dat,assay = 'Spatial',
                                vars.to.regress = c('percent.mt','S.Score','G2M.Score'),
                                verbose = FALSE)
    return(Folding_dat)
}
