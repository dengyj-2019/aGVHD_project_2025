library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Hs.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)


load('/data1/02.private/dengyj/analysis/database/human_cc_genes.rds')
load('/data1/02.private/dengyj/analysis/aGVHD/doublets/aGVHD_step_1.rds')
#load('/data1/02.private/dengyj/analysis/aGVHD/doublets/aGVHD_relieved_step_1.rds')
library(harmony)
library(openxlsx)

library(patchwork)

source('/data1/02.private/dengyj/analysis/mycode/Cal_enrichment.R')

library(AUCell)
library(SCopeLoomR)
library(SCENIC)

library(clusterProfiler)

#library(rliger)
library(SeuratWrappers)


mycolors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
'#968175')

cc_genes <- readRDS("/data1/02.private/dengyj/analysis/thymus_multiomics/support/cc_genes/cc_genes.rds")

SM <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/SM/humanSM.rds')
TF <- readRDS('/data1/02.private/dengyj/analysis/database/TF_SM/human/TF/humanTF.rds')

dim(aGVHD_step_1)





aGVHD_step_1[["percent.mt"]] <- PercentageFeatureSet(aGVHD_step_1, pattern = "MT-")

options(repr.plot.width=7, repr.plot.height=7)
VlnPlot(aGVHD_step_1, features = 'percent.mt', pt.size = 0)

aGVHD_step_1 <- subset(aGVHD_step_1, percent.mt < 20)





ctrl_val <- min(sapply(cc_genes, function(y)  length(intersect(y, rownames(aGVHD_step_1)))))

aGVHD_step_1 <- NormalizeData(aGVHD_step_1, verbose = T)
aGVHD_step_1 <- FindVariableFeatures(aGVHD_step_1,selection.method = "vst", nfeatures = 2000)
aGVHD_step_1 <- CellCycleScoring(aGVHD_step_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
aGVHD_step_1$CC.Difference <- aGVHD_step_1$S.Score - aGVHD_step_1$G2M.Score 
aGVHD_step_1@meta.data[, names(cc_genes)] <- NULL
for(i in names(cc_genes)){
    aGVHD_step_1 <- AddModuleScore(aGVHD_step_1, features = cc_genes[i], name = i, 
                                  ctrl = ctrl_val
                                  )
    colnames(aGVHD_step_1@meta.data) <- gsub(paste0(i, '1'), i, 
                                          colnames(aGVHD_step_1@meta.data))
}

aGVHD_step_1 <- ScaleData(aGVHD_step_1, vars.to.regress = names(cc_genes))
aGVHD_step_1 <- RunPCA(aGVHD_step_1, npcs = 50)

aGVHD_step_1$time <- 'first'
aGVHD_step_1$time[grep('6|7|8', aGVHD_step_1$batch)] <- 'second'
aGVHD_step_1$time[grep('DHY1008', aGVHD_step_1$batch)] <- 'third'
aGVHD_step_1$time[aGVHD_step_1$batch %in% c('WXL0702','DHY1022','LWM0818', 'YX0708')] <- 'fourth'
aGVHD_step_1$time[aGVHD_step_1$batch %in% c('P8_resistance')] <- 'fifth'



aGVHD_step_1$type <- 'null'
aGVHD_step_1$type[grep('C|LWM0818|YX0708', aGVHD_step_1$batch)] <- 'C'
aGVHD_step_1$type[grep('P|DHY1008', aGVHD_step_1$batch)] <- 'P'
aGVHD_step_1$type[aGVHD_step_1$batch %in% c('P8_resistance')] <- 'null'

table(aGVHD_step_1$type, aGVHD_step_1$batch)

aGVHD_step_1$hormone_type <- 'null'
#sensitive
aGVHD_step_1$hormone_type[aGVHD_step_1$type %in% 'P' & aGVHD_step_1$batch %in% c('P6', 'DHY1008')] <- 'S'
#resistent
aGVHD_step_1$hormone_type[aGVHD_step_1$type %in% 'P' & !aGVHD_step_1$batch %in% c('P6', 'DHY1008')] <- 'R'

aGVHD_step_1$organ_type <- 'null'
#intestine
aGVHD_step_1$organ_type[aGVHD_step_1$type %in% 'P' & aGVHD_step_1$batch %in%  c('P1', 'P7', 'P8')] <- 'Intestine'
#intestine+other target
aGVHD_step_1$organ_type[aGVHD_step_1$type %in% 'P' & !aGVHD_step_1$batch %in% c('P1', 'P7', 'P8')] <- 'Other'

aGVHD_step_1$time_type <- 'null'
#former
aGVHD_step_1$time_type[aGVHD_step_1$batch %in% c('P8','DHY1008')] <- 'former'
#latter
aGVHD_step_1$time_type[aGVHD_step_1$batch %in% c('WXL0702','DHY1022')] <- 'latter'

aGVHD_step_1 <- RunHarmony(aGVHD_step_1, c('batch', 'time', 'type'))

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(aGVHD_step_1, ndims = 50, reduction = 'harmony')

print(aGVHD_step_1[['harmony']], dims = 1:30, nfeatures = 30)


options(repr.plot.width=6, repr.plot.height=6)
for(i in 1:(30/2)){
    iter_num <- c((2*i-1),2*i)
    p <- DimPlot(aGVHD_step_1, reduction = "harmony", group.by = "Phase", dims = iter_num)
    print(p)
}


aGVHD_step_1 <- RunUMAP(aGVHD_step_1, reduction = "harmony", dims = 1:20)

mycolors <- c('#9E6BAB','#CFE3C6','#E3CEE4','#F3746C','#86382A','#ABA300','#D6D4EB','#B9A96B','#FF68A1',
'#EAA944','#7AAF93','#7A8DBB','#7673AE','#396F68','#ECAFCF','#EACF68','#F7DDD4','#EBDBE4',
'#66CEF6','#F8F4A8','#C35338','#EF5276','#A0D7C9','#63B472','#F9DBE5','#0CB702','#F48930',
'#6B6A6B','#27BDCF','#F8BFAF','#F5C785','#DEEAB1','#217CB8','#31FEB3','#74517B','#588198',
'#CAA57D','#9C99C4','#2D563D','#FF77AB','#9F8CFF','#D5E7F7','#22A3FF','#00E8F7','#BB4A94',
'#69B4CE','#C9BDB2','#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87',
'#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA',
'#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398','#AA9A59','#E63863','#E39A35','#C1E6F3',
'#6778AE','#91D0BE','#B53E2B','#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180',
'#3A6963','#968175')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(aGVHD_step_1, reduction = 'umap', group.by = 'Phase')
DimPlot(aGVHD_step_1, reduction = 'umap', group.by = 'batch')+
scale_color_manual(values = mycolors)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(aGVHD_step_1, reduction = 'umap', group.by = 'type')

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(aGVHD_step_1, features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(aGVHD_step_1, features = c('MPO', 'ELANE','AZU1'),
            cols = c("lightgrey", 'red'))

orig_meta <- readRDS('./old/220401/old_metadata.rds')
aGVHD_step_1 <- AddMetaData(aGVHD_step_1, orig_meta[, 'Identity', drop =F], col.name = 'tmp_idents')
aGVHD_step_1$tmp_idents <-as.character(aGVHD_step_1$tmp_idents)
aGVHD_step_1$tmp_idents[is.na(aGVHD_step_1$tmp_idents)] <- 'none'

options(repr.plot.width=12, repr.plot.height=9)
DimPlot(aGVHD_step_1[, !aGVHD_step_1$batch %in% c('DHY1008','WXL0702','DHY1022','LWM0818', 'YX0708','P8_resistance')], 
        reduction = 'umap', group.by = 'tmp_idents', label = T)+
scale_color_manual(values = mycolors)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(aGVHD_step_1, reduction = 'umap', group.by = 'batch', label = T)+
scale_color_manual(values = mycolors)

options(repr.plot.width=9, repr.plot.height=9)
DimPlot(aGVHD_step_1, reduction = 'umap', group.by = 'Phase', label = T)

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(aGVHD_step_1, features = c('CD3E', 'CD19', 'SPI1', 'TYROBP'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(aGVHD_step_1, features = c('HBA1', 'PPBP', 'CD34', 'ELANE'),
            cols = c("lightgrey", 'red'))





aGVHD_step_1 <- FindNeighbors(aGVHD_step_1, reduction = "harmony", dims = 1:20)

aGVHD_step_1 <- FindClusters(aGVHD_step_1, resolution = seq(0.1, 1, 0.05))

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(aGVHD_step_1, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 

options(repr.plot.width=12, repr.plot.height=8)
VlnPlot(aGVHD_step_1, features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'), group.by = 'RNA_snn_res.0.2', 
       pt.size = 0)



Idents(aGVHD_step_1) <- aGVHD_step_1$RNA_snn_res.0.2

options(repr.plot.width=7.5, repr.plot.height=7.5)
DimPlot(aGVHD_step_1, reduction = 'umap',label = T)



#######extract T and NK for further clustering
T_NK <- aGVHD_step_1[, Idents(aGVHD_step_1) %in% c('0','3', '1')]

orig_T_NK_idents <- readRDS('old/220808/orig_T_NK_idents.rds')

T_NK <- AddMetaData(T_NK, metadata = orig_T_NK_idents, col.name = 'orig_T_NK_idents')


T_NK <- NormalizeData(T_NK)
T_NK <- FindVariableFeatures(T_NK, selection.method = "vst", nfeatures = 2000)
T_NK <- ScaleData(T_NK, vars.to.regress = c(names(cc_genes)))
T_NK <- RunPCA(T_NK, npcs = 50)

T_NK <- RunHarmony(T_NK, c('batch', 'time', 'type'), max.iter.harmony = 10)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(T_NK, ndims = 50, reduction = 'harmony')

T_NK <- RunUMAP(T_NK, reduction = "harmony", dims = 1:20)
T_NK <- FindNeighbors(T_NK, reduction = "harmony", dims = 1:20)

T_NK <- FindClusters(T_NK, resolution = seq(0.1, 1, 0.05))

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK, reduction = 'umap', group.by = 'Phase')
DimPlot(T_NK, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'
                        ))     

options(repr.plot.width=8, repr.plot.height=8)
p <- DimPlot(T_NK, reduction = 'umap', label = T, group.by = 'orig_T_NK_idents')
p

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'CD3E', 'CD3D', 'CD3G', 'CD247', 'CD8A', 'CD4', 'CD8B'
                        ))

options(repr.plot.width=16, repr.plot.height=16)####
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'PF4', 'PPBP',  'TUBB1', 'ITGA2B'
                        ))


options(repr.plot.width=16, repr.plot.height=16)####
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'ISG15', 'MX1', 'IFIT1', 'OAS1'
                        ))

options(repr.plot.width=16, repr.plot.height=16)####细胞比较多的类群，容易造成contamination
FeaturePlot(T_NK, reduction = 'umap', cols = c('lightgrey', 'red'),
            features = c( 'PF4', 'PPBP',  'HBA1', 'MPO', 'GATA1', 'VCAN', 'CSF3R', 'ELANE', 
                         'CD14', 'FCGR3A','CD1C','SPI1', 
                         'CD247', 'CD3E', 'CD34', 'SPINK2'
                        ))






Idents(T_NK) <- T_NK$RNA_snn_res.0.1
T_NK_2 <- T_NK[, !Idents(T_NK) %in% c('3')]###removing doublets
T_NK_2 <- NormalizeData(T_NK_2)
T_NK_2 <- FindVariableFeatures(T_NK_2, selection.method = "vst", nfeatures = 2000)
T_NK_2 <- ScaleData(T_NK_2, vars.to.regress = c(names(cc_genes)))
T_NK_2 <- RunPCA(T_NK_2, npcs = 50)

T_NK_2 <- RunHarmony(T_NK_2, c('batch', 'time', 'type'), max.iter.harmony = 10)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(T_NK_2, ndims = 50, reduction = 'harmony')

T_NK_2 <- RunUMAP(T_NK_2, reduction = "harmony", dims = 1:25)
T_NK_2 <- FindNeighbors(T_NK_2, reduction = "harmony", dims = 1:25)

T_NK_2 <- FindClusters(T_NK_2, resolution = seq(0.1, 1, 0.05))

options(repr.plot.width=7.5, repr.plot.height=7.5)
for(i in seq(0.1, 1, 0.05)){
    group_iter <- paste0('RNA_snn_res.', i)
    p <- DimPlot(T_NK_2, reduction = 'umap',label = T, group.by = group_iter) + 
    ggtitle(paste0('Resolution: ', i))
    print(p)
} 

options(repr.plot.width=10, repr.plot.height=8)
p <- DimPlot(T_NK_2, reduction = 'umap', label = T, group.by = 'orig_T_NK_idents')
p+scale_color_manual(values=mycolors)
p


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CD4', 'CD40LG', 'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G', 'CD247'))

options(repr.plot.width=8, repr.plot.height=8)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CD34'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CDK1', 'MKI67','CCNA2'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('IFIT1', 'MX1','ISG15', 'OAS3'))

options(repr.plot.width=16, repr.plot.height=16)###
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CCR7', 'LEF1', 'SELL', 'SOX4', 'ANXA2', 'ANXA1', 'ACTN1', 
                        'CD55', 'FAS'))


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c( 'GZMK','GZMB', 'CD27'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('HSPD1', 'HSPE1','HSP90AB1'))

options(repr.plot.width=16, repr.plot.height=16)###NK/NKT
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('TYROBP', 'NCR1', 'NCAM1', 'TRDC'))

options(repr.plot.width=16, repr.plot.height=16)###gdT
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('TRDC', 'TRDV1', 'TRDV2', 'TRGV9'))

options(repr.plot.width=16, repr.plot.height=16)###
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('GZMK', 'GZMB', 'GNLY', 'CCL5', 'PRF1', 'GZMA', 'GZMH'))

options(repr.plot.width=16, repr.plot.height=16)###MAIT
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('SLC4A10', 'RORC', 'RORA', 'NCR3', 'CD8A', 'CD8B', 
                        'ZBTB16', 'KLRB1', 'CD4'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CD3D'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CD3G'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('CD3E'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c('KLRB1'))

Idents(T_NK_2) <- T_NK_2$RNA_snn_res.0.9


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(T_NK_2, reduction = 'umap', label = T)


####extract naive T-cell to distinguish CCD4 and CD8
TN_like <- T_NK_2[, Idents(T_NK_2) %in% c('4','11','18')]
TN_like <- NormalizeData(TN_like)
TN_like <- FindVariableFeatures(TN_like, selection.method = "vst", nfeatures = 2000)
TN_like <- ScaleData(TN_like, vars.to.regress = c(names(cc_genes), 'batch'))
TN_like <- RunPCA(TN_like, npcs = 50)



options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(TN_like, ndims = 50, reduction = 'pca')

# options(repr.plot.width=16, repr.plot.height=16)
# for(i in 1:(20/2)){
#     iter_num <- c((2*i-1),2*i)
#     p <- FeaturePlot(TN_like, reduction = "pca", features = c('CD4', 'CD8A', 'CD8B', 'CD40LG'),
#                      dims = iter_num)
#     print(p)
# }

TN_like <- RunUMAP(TN_like, reduction = "pca", dims = c(1:10))
TN_like <- FindNeighbors(TN_like, reduction = "pca", dims = c(1:10))

TN_like <- FindClusters(TN_like, resolution = 0.7)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(TN_like, reduction = 'umap', group.by = 'Phase')
DimPlot(TN_like, reduction = 'umap', label=T)
DimPlot(TN_like, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(TN_like, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(TN_like, cols = c('lightgrey', 'red'),
            features = c( 'CD4', 'CD40LG', 'CD8A','CD8B','CD27', 'CCR7', 'LEF1', 'TCF7', 'GZMK'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(TN_like, cols = c('lightgrey', 'red'),
            features = c( 'ACTN1', 'FAS', 'LRRN3'))


options(repr.plot.width=7, repr.plot.height=7)
TN_like$tmp_idents <- as.character(TN_like$seurat_clusters)
TN_like$tmp_idents[TN_like$seurat_clusters %in% c('2','4','5','6')] <- 'CD8+TN'
TN_like$tmp_idents[!TN_like$seurat_clusters %in% c('2','4','5','6')] <- 'CD4+TN'
Idents(TN_like) <- TN_like$tmp_idents
DimPlot(TN_like, label =T)

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c( 'IL2', 'IL3', 'CSF2'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c( 'IFNG', 'TBX21', 'GZMK', 'TNF'))



options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c( 'RORC','CCR6','IL23R', 'TBX21','IFNG',
             'MIR155HG','TNF'))


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(T_NK_2, cols = c('lightgrey', 'red'),
            features = c( 'ITGB1', 'PASK',  'AHNAK',  'KLRB1'))


Th_like <- T_NK_2[, Idents(T_NK_2) %in% c('9')]

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Th_like, ndims = 50, reduction = 'harmony')

Th_like <- RunUMAP(Th_like, reduction = "harmony", dims = 1:10)
Th_like <- FindNeighbors(Th_like, reduction = "harmony", dims = 1:10)


Th_like <- FindClusters(Th_like, resolution = 0.8)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(Th_like, reduction = 'umap', label=T)
DimPlot(Th_like, reduction = 'umap', group.by = 'Phase')
DimPlot(Th_like, reduction = 'umap', group.by = 'batch')


options(repr.plot.width=7, repr.plot.height=7)
p <- DimPlot(Th_like, reduction = 'umap', group.by = 'orig_T_NK_idents')
p
p1 <- p
p1$data <- p1$data[p1$data$orig_T_NK_idents != 'none', ]
p1+
scale_color_manual(values = mycolors)

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'TBX21', 'IFNG', 'GZMK', 'GZMB'))  

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'CD8A', 'CD8B', 'CD4', 'CD40LG'))  

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'LIMS1', 'PASK', 'ITGB1', 'AHNAK'))  

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'TCF7', 'CCR7', 'LEF1', 'PASK', 'SELL', 'IL7R', 'CCL5', 'NKG7'))  

options(repr.plot.width=8, repr.plot.height=8)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'FAS','PRF1','KLRC1'))  



options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(Th_like, cols = c('lightgrey', 'red'),
            features = c( 'CCL5', 'GZMK', 'NKG7', 'ANXA1'))  

options(repr.plot.width=7, repr.plot.height=7)
Th_like$tmp_idents <- as.character(Th_like$seurat_clusters)
Th_like$tmp_idents[!Th_like$seurat_clusters %in% c('1', '3')] <- 'CD4+CXCR3+T'
Th_like$tmp_idents[Th_like$seurat_clusters %in% c('1','3')] <- 'CD4+TCM'
Idents(Th_like) <- Th_like$tmp_idents
DimPlot(Th_like, label =T)







CD8_TEM_like <- T_NK_2[, Idents(T_NK_2) %in% c('1','6','10','13','14')]
CD8_TEM_like <- NormalizeData(CD8_TEM_like)
CD8_TEM_like <- FindVariableFeatures(CD8_TEM_like, selection.method = "vst", nfeatures = 2000)
CD8_TEM_like <- ScaleData(CD8_TEM_like, vars.to.regress = c(names(cc_genes)))
CD8_TEM_like <- RunPCA(CD8_TEM_like, npcs = 50)

CD8_TEM_like <- RunHarmony(CD8_TEM_like, c('batch', 'time', 'type'), max.iter.harmony = 10)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD8_TEM_like, ndims = 50, reduction = 'harmony')

CD8_TEM_like <- RunUMAP(CD8_TEM_like, reduction = "harmony", dims = 1:20)
CD8_TEM_like <- FindNeighbors(CD8_TEM_like, reduction = "harmony", dims = 1:20)

CD8_TEM_like <- FindClusters(CD8_TEM_like, resolution = 0.4)#0.8

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_TEM_like, reduction = 'umap', label=T)
DimPlot(CD8_TEM_like, reduction = 'umap', group.by = 'Phase')
DimPlot(CD8_TEM_like, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'CD8A', 'CD8B', 'CD4', 'CD40LG'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'CCR7', 'LEF1', 'GZMK', 'TCF7', 'CD27', 'SELL', 'IL7R', 'CCL5'))  


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'C1orf162','ITGB1','MAL','CAPG','TNFSF10','GATA3','LMNA','RCBTB2','FXYD7'))  

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'CX3CR1', 'GZMB', 'FGFBP2'))  


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'IFNG', 'TBX21', 'CXCR3'))  


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEM_like, cols = c('lightgrey', 'red'),
            features = c( 'PASK', 'AHNAK', 'ITGB1'))  


options(repr.plot.width=7, repr.plot.height=7)
p <- DimPlot(CD8_TEM_like, reduction = 'umap', group.by = 'tmp_idents')
p
p1 <- p
p1$data <- p1$data[p1$data$tmp_idents != 'none', ]
p1+
scale_color_manual(values = mycolors)





options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_TEM_like, reduction = 'umap', label = T)





options(repr.plot.width=7, repr.plot.height=7)
CD8_TEM_like$tmp_idents <- as.character(CD8_TEM_like$seurat_clusters)
CD8_TEM_like$tmp_idents[CD8_TEM_like$seurat_clusters %in% c('6')] <- 'doublets'
#CD8_TEM_like$tmp_idents[CD8_TEM_like$seurat_clusters %in% c('0')] <- 'CD8+TEFF'
CD8_TEM_like$tmp_idents[CD8_TEM_like$seurat_clusters %in% c('1')] <- 'CD8+TCM'
CD8_TEM_like$tmp_idents[CD8_TEM_like$seurat_clusters %in% c('7', '5')] <- 'low_quality'
CD8_TEM_like$tmp_idents[CD8_TEM_like$seurat_clusters %in% c('0','2','4','7','3')] <- 'CD8+TEM'
Idents(CD8_TEM_like) <- CD8_TEM_like$tmp_idents
DimPlot(CD8_TEM_like, label =T)

options(repr.plot.width=7, repr.plot.height=7)

DotPlot(CD8_TEM_like, 
        features =c('GZMK', 'CXCR3', 'CXCR6', 'NKG7', 'HLA-DRB1','HLA-DRA',
                    'CCR5', 'GZMM', 'IFNG', 'TBX21', 'CCL5', 'PDCD1', 'LAG3', 
                   'CCR7', 'TCF7', 'LEF1'))+
theme(axis.text.x =element_text(angle=45, hjust = 1))






gdT <- T_NK_2[, Idents(T_NK_2) %in% c('12')]

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(gdT, ndims = 50, reduction = 'harmony')

gdT <- RunUMAP(gdT, reduction = "pca", dims = 1:30)
gdT <- FindNeighbors(gdT, reduction = "pca", dims = 1:30)


options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(gdT, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(gdT, cols = c('lightgrey', 'red'),
            features = c( 'TRDV1', 'NCR1','TRDV2', 'TRGV9'))

gdT <- FindClusters(gdT, resolution = 0.1)#0.8

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(gdT, reduction = 'umap', label=T)
DimPlot(gdT, reduction = 'umap', group.by = 'Phase')
DimPlot(gdT, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=7, repr.plot.height=7)
gdT$tmp_idents <- as.character(gdT$seurat_clusters)
gdT$tmp_idents[gdT$seurat_clusters %in% c('0')] <- 'Vδ2_T'
gdT$tmp_idents[!gdT$seurat_clusters %in% c('0')] <- 'Vδ1_T'

Idents(gdT) <- gdT$tmp_idents
DimPlot(gdT, label =T)



NK_like <- T_NK_2[, Idents(T_NK_2) %in% c('15')]
# NK_like <- NormalizeData(NK_like)
# NK_like <- FindVariableFeatures(NK_like, selection.method = "vst", nfeatures = 2000)
# NK_like <- ScaleData(NK_like, vars.to.regress = c(names(cc_genes), 'batch'))
# NK_like <- RunPCA(NK_like, npcs = 50)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(NK_like, ndims = 50, reduction = 'harmony')

NK_like <- RunUMAP(NK_like, reduction = "harmony", dims = 1:20)
NK_like <- FindNeighbors(NK_like, reduction = "harmony", dims = 1:20)

NK_like <- FindClusters(NK_like, resolution = 0.4)#0.8


options(repr.plot.width=7, repr.plot.height=7)
DimPlot(NK_like, reduction = 'umap', label=T)
DimPlot(NK_like, reduction = 'umap', group.by = 'Phase')
DimPlot(NK_like, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(NK_like, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(NK_like, cols = c('lightgrey', 'red'),
            features = c( 'CD3D', 'CD3G', 'CD3E', 'CD247'))

options(repr.plot.width=6, repr.plot.height=6)##
DotPlot(NK_like, cols = c('lightgrey', 'red'),
            features = c( 'CD3D', 'CD3G', 'CD3E', 'CD247'))

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(NK_like, cols = c('lightgrey', 'red'),
            features = c( 'NCR1', 'TYROBP', 'ITGAX', 'KLRB1'))



CD8_TEFF_like <- RunHarmony(CD8_TEFF_like, c('batch', 'time', 'type'), max.iter.harmony = 10)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(CD8_TEFF_like, ndims = 50, reduction = 'harmony')

CD8_TEFF_like <- RunUMAP(CD8_TEFF_like, reduction = "harmony", dims = 1:5)
CD8_TEFF_like <- FindNeighbors(CD8_TEFF_like, reduction = "harmony", dims = 1:5)

CD8_TEFF_like <- FindClusters(CD8_TEFF_like, resolution = 0.1)

options(repr.plot.width=7, repr.plot.height=7)
DimPlot(CD8_TEFF_like, reduction = 'umap', label=T)
DimPlot(CD8_TEFF_like, reduction = 'umap', group.by = 'Phase')
DimPlot(CD8_TEFF_like, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEFF_like, cols = c('lightgrey', 'red'),
            features = c( 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))  

options(repr.plot.width=16, repr.plot.height=16)##
FeaturePlot(CD8_TEFF_like, cols = c('lightgrey', 'red'),
            features = c( 'GZMB', 'GZMK', 'FGFBP2', 'CX3CR1'))  

options(repr.plot.width=7, repr.plot.height=7)
Idents(CD8_TEFF_like) <- CD8_TEFF_like$seurat_clusters
CD8_TEFF_like <- RenameIdents(CD8_TEFF_like, '0' = 'CD8+TEFF', '1'='CD8+TEM')
DimPlot(CD8_TEFF_like, label =T)   





options(repr.plot.width=7, repr.plot.height=7)
DimPlot(T_NK_2, label=T)

options(repr.plot.width=7, repr.plot.height=7)
Idents(T_NK_2) <- T_NK_2$`RNA_snn_res.0.9`
T_NK_2 <- RenameIdents(T_NK_2, '0' = 'CD8+TEFF',#'2' = 'CD8+TEFF',
                       '3' = 'CD8+TEFF', '5' = 'CD4+TNF+T','7' = 'CD4+TNF+T',
                      '8' = 'Treg', '17'='MAIT','16'='IFN_T')
DimPlot(T_NK_2, label =T)   

Tcell <- merge(T_NK_2[, Idents(T_NK_2) %in% c('CD8+TEFF', 'CD4+TNF+T','IFN_T', 'Treg', 'MAIT')], 
               list(Th_like,#[, !Idents(Th_like) %in% c('doublets')],                 
                    TN_like,
                    gdT, 
                    NK_like[, Idents(NK_like) %in% c('NKT')],
                    CD8_TEFF_like,
                    CD8_TEM_like[, !Idents(CD8_TEM_like) %in% c('doublets', 'low_quality')]))

Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = c(names(cc_genes)))
Tcell <- RunPCA(Tcell, npcs = 50)


Tcell <- RunHarmony(Tcell, c('batch', 'time', 'type'), max.iter.harmony = 10)

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(Tcell, ndims = 50, reduction = 'harmony')

Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:30)

options(repr.plot.width=7, repr.plot.height=7)
#DimPlot(Tcell, reduction = 'umap', label=T)
DimPlot(Tcell, reduction = 'umap', group.by = 'Phase')
DimPlot(Tcell, reduction = 'umap', group.by = 'batch')
DimPlot(Tcell, reduction = 'umap', group.by = 'type')

cluster_colors <- c("CD4+TCM"="#E5D2DD",#"CD8+TEM2"="#53A85F",
                    #"CD8+TEM3"="#F1BB72",#"CD8+TEM4"="#F3B1A0",
                    "Vδ2_T"="#E59CC4",
                    "CD8+TEFF"="#BD956A",
"CD8+TEM"="#8C549C","Th17"="#B53E2B","Treg"="#9FA3A8",'Vδ1_T'="#E0D4CA",
                    "CD4+TN"="#956da1","Th1"="#C5DEBA",
"CD8+TN"="#58A4C3","NKT"="#E4C755","CD8+TCM"="#F7F398","MAIT"="#AA9A59","CD4+TEM"="#E63863",
                    "IFN_T"="#E39A35",
                   'CD8+TEM_1' = '#de6f54', 'CD8+TEM_2' = '#4e77a4', 'CD8+TEM_3'='#04bec3', 
                   'CD4+CXCR3+T' = '#C5DEBA', 'CD4+TNF+T'= '#B53E2B')

levels(Tcell) <- c('CD8+TN','CD8+TCM','CD8+TEM','CD8+TEFF','MAIT',
'CD4+TN','CD4+TCM','CD4+TNF+T','CD4+CXCR3+T','Treg','IFN_T','Vδ1_T','Vδ2_T',
'NKT')


options(repr.plot.width=9, repr.plot.height=9)
DimPlot(Tcell, reduction = 'umap', label=F, repel=F)+
scale_color_manual(values = cluster_colors)

Tcell$Identity <- as.character(Idents(Tcell))



saveRDS(Tcell, file = 'Tcell.rds')

save(T_NK, T_NK_2, Th_like,#Th_like_2,
     TN_like,
     gdT,NK_like, CD8_TEM_like,CD8_TEFF_like,
     file = 'T_NK.rds')

saveRDS(aGVHD_step_1, file = 'aGVHD_step_1.rds')
