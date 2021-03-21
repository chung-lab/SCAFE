#!/usr/bin/env Rscript
suppressMessages({
  library(docopt)
})

doc <- "Usage: get_cicero_cds.R [options]

options:
--count MATRIX       path to the matrix folder
--peaks PEAK         path to the peak bed
--min MIN            Minimum cells CRE must be present in
--binarize BIN       True or False: binarize matrix
--out OUT            Output folder
--help               Show this help text"

opt <- docopt(doc) # docopt parsing
#print(opt)

timestamp()
options(scipen=999)
cat('Loading packages ...\n')

suppressMessages({
  library(monocle3)
  library(cicero)
  library(Matrix)
  library(data.table)
  library(scales)
})


output.dir <- opt$out
if(!dir.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
count <- opt$count
peaks <- opt$peaks
min_cells <- as.numeric(opt$min)
binarize <- as.logical(opt$binarize)
min_cre = 1
sample_n = 100

cat('Loading Data ...\n')

# read matrix
indata <- readMM(file=paste0(count,'/matrix.mtx'))
colnames(indata) <- scan(paste0(count,'/barcodes.tsv'), what='character')
rownames(indata) <- read.table(paste0(count,'/genes.tsv'), header=F)$V1
barcodes <- colnames(indata)

# binarize the matrix
if (binarize) {
  indata@x[indata@x > 0] <- 1
}

# format peak info
peakinfo <- fread(peaks) 
setkey(peakinfo, V4)
peakinfo <- as.data.frame(peakinfo[rownames(indata), 1:4, with = F])
names(peakinfo) <- c("chr","bp1","bp2","name")
row.names(peakinfo) <- peakinfo$name

#Cell metadata
cellinfo <- data.frame(cells = barcodes)
rownames(cellinfo) <- barcodes

#Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata, cell_metadata = cellinfo, gene_metadata = peakinfo))

#Remove peaks in less than n cells
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds) > 0) > min_cells,]

#remove cells less than n cre
input_cds <- input_cds[, Matrix::colSums(exprs(input_cds) > 0) >= min_cre] 

timestamp()
cat('Running UMAP ...\n')
# monocle3 umap construction
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")

# Create metacells
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
meta_cell_expr <- as.matrix(assay(cicero_cds, 1))
meta_cell_expr <- round(meta_cell_expr, digits = 5)
write.table(meta_cell_expr, file=gzfile(paste0(output.dir, '/meta_cell_expr.tsv.gz')), quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

saveRDS(cicero_cds, paste0(output.dir,"/cicero_cds.RDS"))
