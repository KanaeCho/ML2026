library(tidyverse)
library(Seurat)
library(harmony)
library(SoupX)
library(patchwork)
library(scDblFinder)

################## load data ########################
LoadData <- function(data_path=NULL, assay='ATAC'){
    if (assay=='ATAC'){
        atac_assay=t(readMM(paste0(data_path,'/matrix.mtx')))
        barcodes=readLines(paste0(data_path,'/barcodes.tsv'))
        peaks=readLines(paste0(data_path,'/peak.bed'))
        row.names(atac_assay)=peaks
        colnames(atac_assay)=barcodes
        #
        return (atac_assay)
    }
    if (assay=='RNA'){
        rna_assay=readMM(paste0(data_path,'/matrix.mtx'))
        barcodes=readLines(paste0(data_path,'/barcodes.tsv'))
        features=readLines(paste0(data_path,'/genes.tsv'))
        row.names(rna_assay)=features
        colnames(rna_assay)=barcodes
        #
        return (rna_assay)
    }
    if (assay=='GeneScore'){
        rna_assay=readMM(paste0(data_path,'/matrix.mtx'))
        barcodes=readLines(paste0(data_path,'/barcodes.tsv'))
        features=readLines(paste0(data_path,'/genes.tsv'))
        row.names(rna_assay)=features
        colnames(rna_assay)=barcodes
        #
        return (rna_assay)
    }
}

####################################### pre-processing #######################################

### remove Ambient RNAs by contamination
runSoupX_Droplets <- function(matrix_path=NULL, expression_matrix=NULL, droplet_matrix=NULL){
  # loading data
  if (!is.null(matrix_path)){
    toc = Seurat::Read10X(file.path(matrix_path, "04.Matrix", "FilterMatrix"), gene.column = 1)
    tod = Seurat::Read10X(file.path(matrix_path, "02.cDNAAnno", "RawMatrix"), gene.column = 1)
  }
  else{
    toc=expression_matrix
    tod=droplet_matrix
  }
  #
  genes=intersect(rownames(toc),rownames(tod))
  toc=toc[genes,]
  tod=tod[genes,]
  sc = SoupChannel(tod, toc)
  # get clusters info
  obj <- CreateSeuratObject(counts = toc)
  obj=Process_RNA(obj)
  soupx_groups = Idents(obj)
  #
  sc = setClusters(sc, soupx_groups)
  sc = autoEstCont(sc, doPlot=FALSE)
  out = adjustCounts(sc,roundToInt = TRUE)
  return (out)
}

runSoupX_NoDroplets <- function(matrix_path=NULL, expression_matrix=NULL){
  # loading data
  if (!is.null(matrix_path)){
    toc = Seurat::Read10X(matrix_path, gene.column = 1)
  }
  else{toc=expression_matrix}
  #
  scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
  # Calculate soup profile
  soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
  scNoDrops = setSoupProfile(scNoDrops, soupProf)
  # get clusters info
  obj <- CreateSeuratObject(counts = toc)
  obj <- Process_RNA(obj)
  soupx_groups = Idents(obj)
  # Set cluster information in SoupChannel
  scNoDrops = setClusters(scNoDrops, soupx_groups)
  # Estimate contamination fraction
  scNoDrops  = autoEstCont(scNoDrops, doPlot=FALSE)
  # Infer corrected table of counts and rount to integer
  out = adjustCounts(scNoDrops, roundToInt = TRUE)
  return (out)
}

############# 创建Seurat对象及QC ###########
# import scRNA data and create seurat object
Create_scRNA_object <- function(data.dir=NULL, rna.assay=NULL){
  #
  if (is.null(data.dir)){
    if (is.null(rna.assay)){stop('Please provide the RNA expression matrix')}
  }
  else{
    rna.assay=Read10X(data.dir = data.dir, gene.column = 1)
    #rna.assay=Read10X(data.dir = data.dir, gene.column = 2)
  }
  #
  proj <- CreateSeuratObject(counts = rna.assay)
  #
  return(proj)
}

# get overlap genes from multiple seurat objects
Overlap_Seurat_Genes <- function(list_seurat=NULL){
  if (!inherits(x = list_seurat, what = "list")) {
    cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }
  gene_list <- lapply(X = list_seurat, FUN = function(x) {
    x <- rownames(x)
  })
  overlap_genes <- purrr::reduce(gene_list, function(x, y) {
    intersect(x, y)
  })
}

Create_scATAC_object <- function(counts=NULL,fragment_file=NULL){
    #
    atac.assay <- CreateChromatinAssay(counts, fragments = fragment_file)
    atac.obj <- CreateSeuratObject(atac.assay, assay = "ATAC")
    #
    return (atac.obj)
}
# Do annotation for seurat object
annotate_ATAC_obj <- function(proj=proj, annotation='hg38'){
    #
    if (annotation=='hg38'){
        hg38 = "/media/AnalysisTempDisk/Caipengfei/workspace/zzj/genome/hg38/genes/genes.gtf"
        gtf <- rtracklayer::import(hg38)
        gene.coords <- gtf[gtf$type == 'gene']
        seqlevelsStyle(gene.coords) <- 'UCSC'
        gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
        gene.coords$gene_biotype <- gene.coords$gene_type
    }
    #
    if (annotation=='mm10'){
        mm10 = "/media/AnalysisTempDisk/Caipengfei/workspace/zzj/genome/mm10/genes.gtf"
        gtf <- rtracklayer::import(mm10)
        gene.coords <- gtf[gtf$type == 'gene']
        seqlevelsStyle(gene.coords) <- 'UCSC'
        gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
        gene.coords$gene_biotype <- gene.coords$gene_type
    }
    # add the gene information to the object
    Annotation(proj) <- gene.coords
    #
    return (proj)
}



## 合并数据框
Merge_Seurat_List <- function(list_seurat=NULL){
  #seurat.obj <- Merge_Seurat_List(list_seurat = seurat.list)
  seurat.obj <- merge(x = seurat.list[[1]], y = seurat.list[c(1:length(seurat.list))[-1]], add.cell.ids = sample_names)
  return(seurat.obj)
}

# Compute QC index of seurat object of RNA
QC_RNA <- function(proj=proj, genome='hg38', assay='RNA'){
  DefaultAssay(proj) <- assay
  #
  if (genome=='hg38'){
    proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
    proj[["percent.ribo"]] <- PercentageFeatureSet(proj, pattern = "^RP[SL]")
  }
  if (genome=='mm10'){
    features=grep(pattern = "^MT-", x = rownames(proj), value = TRUE,ignore.case = TRUE)
    proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^mt-",features=features)
    proj[["percent.ribo"]] <- PercentageFeatureSet(proj, pattern = "^Rp[sl]")
  }
  #
  return(proj)
}

# detect outliers by MAD method
is_outlier <- function(seurat_obj, metric=NULL, nmads=NULL) {
  # Check if metric is a valid column name
  if (!metric %in% colnames(seurat_obj@meta.data)) {
    stop("metric must be a valid column name")
  }
  # Get the values of the metric column
  M <- seurat_obj@meta.data[,metric]
  
  # Calculate the median of the M values
  median <- median(M)
  
  # Calculate the MAD of the M values
  mad <- mad(M)
  
  # Check if any of the M values are more than nmads MADs away from the median
  outlier <- (M < median - nmads * mad) | (median + nmads * mad < M)
  
  # Return the outlier vector
  return(outlier)
}

# remove unwanted genes for clustering purposes
GetBlackListGenes <- function(proj=proj, MT=TRUE, Ribo=TRUE, Cellcycle=TRUE, Sex=TRUE){
  # mitochondrial:
  mt.genes <- grep(pattern = "^MT-", x = rownames(proj), value = TRUE)
  # ribosome genes:
  ribosome.genes <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA", rownames(proj), value=TRUE) 
  #  Cell cycle:
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  # X/Y chromosome genes:
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  geneGR <- GenomicFeatures::genes(txdb)
  sexGenesGR <- geneGR[seqnames(geneGR) %in% c("chrY", "chrX")]
  matchedGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db, keys = sexGenesGR$gene_id, 
                                              columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  sexChr.genes <- matchedGeneSymbols$SYMBOL
  # Genes to ignore (just for clustering purposes)
  vector_genes=c(MT,Ribo,Sex,Cellcycle,Cellcycle)
  order_num=c()
  for (i in 1:length(vector_genes)){
    if (isTRUE(vector_genes[i])){order_num<-c(order_num, i)}
  }
  #print (order_num)
  blacklist.genes <- list(mt.genes, ribosome.genes, sexChr.genes, s.genes, g2m.genes)[order_num]
  blacklist.genes <- do.call(c, blacklist.genes)
  #
  return (blacklist.genes)
}

# Plot QC index of seurat object
Plot_QC_RNA <- function(proj=proj, group=NULL, width=12, height=4, add_stats=TRUE){
  options(repr.plot.width=width, repr.plot.height=height)
  #
  if(add_stats){
    p=VlnPlot(object = proj,features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.ribo'),
              ncol = 2, pt.size = 0,group.by = group, combine = FALSE) 
    stats=stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95)
    #
    p1=p[[1]]+stats + theme(legend.position = 'none')
    p2=p[[2]]+stats + theme(legend.position = 'none')
    p3=p[[3]]+stats + theme(legend.position = 'none')
    p4=p[[4]]+stats + theme(legend.position = 'none')
    return (p1+p2+p3+p4)
  }
  else{
    VlnPlot(object = proj,features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.ribo'),
            ncol = 4, pt.size = 0,group.by = group, combine = TRUE)
  }
}

# Plot QC index of seurat object
Plot_QC_RNA2 <- function(proj=proj, group=NULL, width=12, height=4, add_stats=TRUE){
  options(repr.plot.width=width, repr.plot.height=height)
  #
  if(add_stats){
    p=VlnPlot(object = proj,features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.ribo'),
              ncol = 1, pt.size = 0,group.by = group, combine = FALSE) 
    stats=geom_boxplot(width=0.1, fill="white", outlier.size = 0, outlier.stroke = 0) 
    #
    p1=p[[1]]+ stats + theme(legend.position = 'none') #+ scale_x_discrete(limits = sample_order)
    p2=p[[2]]+ stats + theme(legend.position = 'none') #+ scale_x_discrete(limits = sample_order)
    p3=p[[3]]+ stats + theme(legend.position = 'none') #+ scale_x_discrete(limits = sample_order)
    p4=p[[4]]+ stats + theme(legend.position = 'none') #+ scale_x_discrete(limits = sample_order)
    return (p1+p2+p3+p4)
  }
  else{
    VlnPlot(object = proj,features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.ribo'),
            ncol = 4, pt.size = 0,group.by = group, combine = TRUE)
  }
}

runDoubletFinder <- function(obj, dims, estDubRate=0.075, ncores=1){
  # Run DoubletFinder on a provided (preprocessed) Seurat object
  # Return the seurat object with the selected pANN parameter and the 
  # DoubletFinder doublet classifications
  
  ### pK Identification (parameter-sweep) ###
  # "pK ~ This defines the PC neighborhood size used to compute pANN (proportion of artificial nearest neighbors), 
  # expressed as a proportion of the merged real-artificial data. 
  # No default is set, as pK should be adjusted for each scRNA-seq dataset"
  
  sweep.res.list <- paramSweep_v3(obj, PCs=dims, sct=FALSE, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  
  # Get expected doublets (DF.classify will identify exactly this percent as doublets!)
  nExp_poi <- round(estDubRate * length(Cells(obj)))
  
  # DoubletFinder:
  obj <- doubletFinder_v3(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL
  
  return(obj)
}

#寻找双细胞

# one sample/library
runscDblFinder <- function(obj, assay='ATAC'){
    #
    counts <- GetAssayData(obj[[assay]],slot = "counts")
    sce <- SingleCellExperiment(list(counts=counts))
    #
    sce <- scDblFinder(sce, clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
    dbl_score <- sce@colData[,c('scDblFinder.class','scDblFinder.score')]
    obj=AddMetaData(obj, metadata = as.data.frame(dbl_score[colnames(obj),c('scDblFinder.class', 'scDblFinder.score')]))
    #
    return (obj)
}

# multiple samples/libraries
runscDblFinder_MultiSamples <- function(obj, assay='ATAC', group='batch'){
    #
    counts <- GetAssayData(obj[[assay]],slot = "counts")
    sce <- SingleCellExperiment(list(counts=counts), colData=DataFrame(obj@meta.data[group]))
    #
    sce <- scDblFinder(sce, samples=group, clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
    dbl_score <- sce@colData[,c('scDblFinder.class','scDblFinder.score')]
    obj=AddMetaData(obj, metadata = as.data.frame(dbl_score[colnames(obj),c('scDblFinder.class', 'scDblFinder.score')]))    
    return (obj)
}

Find_doublet1 <- function(data,dims,Doubletrate=0.05){
  # 寻找最优pk值
  sweep.res.list <- paramSweep(data, PCs = dims, sct = FALSE) # 若使用SCT方法 标准化则'sct=T'
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk <-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  message(sprintf("Using pK = %s...", pk))
  
  # 期望doublet数量
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) #可使用注释好的细胞类型
  nExp_poi <- round(Doubletrate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # 鉴定doublets
  data <- doubletFinder(data, PCs = dims, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  pann <- grep(pattern="^pANN", x=names(data@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(data@meta.data), value=TRUE)
  data$pANN <- data[[pann]]
  data[[pann]] <- NULL
  data
}

Find_doublet2 <- function(data,dims,Doubletrate=0.05){
  # 寻找最优pk值
  sweep.res.list <- paramSweep_v3(data, PCs = dims, sct = FALSE) # 若使用SCT方法 标准化则'sct=T'
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk <-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  message(sprintf("Using pK = %s...", pk))
  
  # 期望doublet数量
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) #可使用注释好的细胞类型
  nExp_poi <- round(Doubletrate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # 鉴定doublets
  data <- doubletFinder_v3(data, PCs = dims, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  pann <- grep(pattern="^pANN", x=names(data@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(data@meta.data), value=TRUE)
  data$pANN <- data[[pann]]
  data[[pann]] <- NULL
  data
}

Doubletfind <- function(proj=proj, nfeatures=2000, pc.nums=50, dims.use=1:50, pca.name = 'pca', umap.name= 'umap', res = 0.8){
  DefaultAssay(proj) <- "RNA"
  #
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- ScaleData(proj)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = pca.name)
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- Find_doublet(proj)
  #
  return(proj)
}

####################################### scRNA processing #######################################

### run clusters using counts table
Process_RNA <- function(proj=proj, nfeatures=2000, pc.nums=50, dims.use=1:50, pca.name = 'pca', umap.name= 'umap', res = 0.8){
  DefaultAssay(proj) <- "RNA"
  #
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- ScaleData(proj)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = pca.name)
  proj <- FindNeighbors(proj, reduction = pca.name, dims = dims.use, verbose = FALSE)
  proj <- FindClusters(proj, resolution = res, algorithm = 3, verbose = FALSE)
  proj <- RunUMAP(proj, reduction = pca.name, dims = dims.use, reduction.name = umap.name)
  #
  return(proj)
}

# integration by Seurat method (normalization and integration methods)
Integration_Seurat <- function(proj=proj, group_by="batch", normalization_method="vst", integration_method='cca', nfeatures=2000, blacklist.genes=NULL){
  #
  DefaultAssay(proj) <- "RNA"
  seurat.list <- SplitObject(proj, split.by = group_by)
  
  ###
  if(normalization_method=="vst"){
    # normalize and identify variable features for each dataset independently
    seurat.list <- lapply(X = seurat.list, FUN = function(x){
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = nfeatures)
    features <- features[!features %in% blacklist.genes]
    # select integration method
    if(integration_method=="rpca"){
      seurat.list <- lapply(X = seurat.list, FUN = function(x){
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")
      integrated <- IntegrateData(anchorset = integration.anchors)
    }
    if(integration_method=="cca"){
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
      integrated <- IntegrateData(anchorset = integration.anchors)
    }
  }
  
  ###
  if(normalization_method=='SCT'){
    #
    #seurat.list <- lapply(X = seurat.list, FUN = SCTransform)
    seurat.list <- lapply(X = seurat.list, FUN = SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = nfeatures)
    features <- features[!features %in% blacklist.genes]
    seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
    #
    if(integration_method=="rpca"){
      seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca")
      integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")
    }
    if(integration_method=="cca"){
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features)
      integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")
    }
  }
  
  ###
  if(normalization_method=='SCT_V2'){
    seurat.list <- lapply(X = seurat.list, FUN = function(x) {x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)})
    features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = nfeatures)
    features <- features[!features %in% blacklist.genes]
    seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
    #
    if(integration_method=="rpca"){
      seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca")
      integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")
    }
    if(integration_method=="cca"){
      integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features)
      integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")
    }
  }
  #
  return(integrated)
}

####################################### basic functions #######################################

# count cell numbers of scATAC/scRNA data
count_cells <- function(data_path=NULL, ATAC=FALSE, RNA=FALSE, standard_ATAC=FALSE, standard_RNA=FALSE){
  if (isTRUE(ATAC)){
    barcode_path1=paste0(data_path,'/out/Peak/barcodes.tsv')
    print (system(paste0("wc ", barcode_path1), intern=TRUE))
  }
  if (isTRUE(RNA)){
    barcode_path1=paste0(data_path,'/04.Matrix/FilterMatrix/barcodes.tsv.gz')
    print (system(paste0("zcat ", barcode_path1, ' | wc '), intern=TRUE))
  }
  if (isTRUE(standard_ATAC)){
    barcode_path1=paste0(data_path,'/barcodes.tsv')
    print (system(paste0("wc ", barcode_path1), intern=TRUE))
  }
  if (isTRUE(standard_RNA)){
    barcode_path1=paste0(data_path,'/barcodes.tsv.gz')
    print (system(paste0("zcat ", barcode_path1, ' | wc '), intern=TRUE))
  }
}

# import sparse matrix for scATAC-seq and scRNA-seq data
LoadData <- function(data_path=NULL, assay='ATAC'){
  if (assay=='ATAC'){
    atac_assay=readMM(paste0(data_path,'/matrix.mtx'))
    barcodes=readLines(paste0(data_path,'/barcodes.tsv'))
    peaks=readLines(paste0(data_path,'/peak.bed'))
    row.names(atac_assay)=peaks
    colnames(atac_assay)=barcodes
    #
    return (atac_assay)
  }
  if (assay=='RNA'){
    rna_assay=readMM(paste0(data_path,'/matrix.mtx'))
    barcodes=readLines(paste0(data_path,'/barcodes.tsv'))
    features=readLines(paste0(data_path,'/features.tsv'))
    row.names(rna_assay)=features
    colnames(rna_assay)=barcodes
    #
    return (rna_assay)
  }
  if (assay=='RNAVelocity'){
    spliced_assay=readMM(paste0(data_path,'/spliced.mtx.gz'))
    unspliced_assay=readMM(paste0(data_path,'/unspliced.mtx.gz'))
    #
    barcodes=readLines(paste0(data_path,'/barcodes.tsv.gz'))
    features=readLines(paste0(data_path,'/features.tsv.gz'))
    
    row.names(spliced_assay)=features
    colnames(spliced_assay)=barcodes
    
    row.names(unspliced_assay)=features
    colnames(unspliced_assay)=barcodes
    #
    return (c(spliced_assay,unspliced_assay))
  }
}

# export 10X format matrix of single-cell data
ExportData_10X_format <- function(data=NULL, out_dir=NULL, assay='RNA', if_compress=TRUE){
  #
  if (is.null(data)){stop('please supply the data')}
  #
  if(assay=='ATAC'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/peaks.bed'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
      R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
      R.utils::gzip(paste0(out_dir,'/peaks.bed'))
      R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  if(assay=='RNA'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/features.tsv'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
      R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
      R.utils::gzip(paste0(out_dir,'/features.tsv'))
      R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  print('Success to output data')
}


# import scRNA data and create seurat object
Create_RNAVelocity_object <- function(spliced_assay=NULL, unspliced_assay=NULL){
  #
  if (is.null(spliced_assay)){
    stop('Please provide the spliced matrix')
  }
  if (is.null(unspliced_assay)){
    stop('Please provide the unspliced matrix')
  }
  #
  proj <- CreateSeuratObject(counts = spliced_assay, assay = "spliced")
  proj[["unspliced"]] <- CreateAssayObject(counts = unspliced_assay)
  #
  return(proj)
}

####################################### Visualization functions #######################################

# plot two variables' overlaps (i.e. celltypes by samples) in a Seurat object
Plot_consistency <- function(proj=proj, group1=NULL, group2=NULL,plot.width=9.5,plot.height=8){
  #
  if (is.null(group1)){stop('please specify the column name of the group1')}
  if (is.null(group2)){stop('please specify the column name of the group2')}
  #
  cluster1=proj@meta.data[,group1]
  names(cluster1)=row.names(proj[[group1]])
  cluster2=proj@meta.data[,group2]
  names(cluster2)=row.names(proj[[group2]])
  predictions <- table(cluster1,cluster2)
  predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
  predictions <- as.data.frame(predictions)
  #
  options(repr.plot.width=plot.width, repr.plot.height=plot.height)
  p1=ggplot(predictions, aes(cluster1, cluster2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", 
                                                                                                   low = "#ffffc8", high = "#7d0025") + xlab("group1") + ylab("group2") + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(p1)
}

confusionMatrix <- function(i = NULL, j = NULL){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}
#
jaccardIndex <- function(mat, i, j){
  # Calculate Jaccard Index between row i and column j in matrix mat
  # (Matrix is an intersection matrix of categories in rows i and columns j)
  AiB <- mat[i,j]
  AuB <- sum(mat[i,]) + sum(mat[,j]) - AiB
  AiB/AuB
}

Plot_confusionMatrix <- function(obj, group1 = NULL, group2 = NULL, plot.width = 8, plot.height = 8){
  #
  if (is.null(group1)){stop('please specify the column name of the group1')}
  if (is.null(group2)){stop('please specify the column name of the group2')}
  #
  cM <- confusionMatrix(paste0(obj@meta.data[,group1]), paste0(obj@meta.data[,group2]))
  cM <- cM / Matrix::rowSums(cM)
  #
  options(repr.plot.width=plot.width, repr.plot.height=plot.height)
  p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = colorRampPalette(c("white", "blue"))(100),
    border_color = "black",fontsize = 14)
  return (p)
}