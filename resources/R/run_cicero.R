#!/usr/bin/env Rscript
suppressMessages({
  library(docopt)
})

doc <- "Usage: scafe_cicero.R [options]

options:
--base BASE          SCAFE run folder
--out OUT            Output folder
--sample SAMPLE      Sample name
--cell CELL          Cell type
--min MIN            Minimum cells CRE must be present in
--binarize BIN       True or False: binarize matrix
--set SET            Which CRE set (e.g. default, lenient, robust)
--chr CHR            File with chromsome sizes
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
  library(Seurat)
  library(data.table)
  library(scales)
})


set <- opt$set
scafe_base <- opt$base
cell <- opt$cell
sample <- opt$sample
chr_sizes <- opt$chr
out <- opt$out
mtx <- paste0(scafe_base,'/count/',set,'/',sample ,'/matrix/')
peaks <- paste0(scafe_base,'/annotate/',set,'/',sample,'/bed/',cell,'.CRE.coord.bed.gz')
anno <- paste0(scafe_base,'/annotate/',set,'/',sample,'/log/',cell,'.CRE.info.tsv.gz')
min_cells <- as.numeric(opt$min)
binarize <- as.logical(opt$binarize)

cat('Loading Data ...\n')
# read matrix
indata <- Read10X(mtx)
barcodes <- colnames(indata)

# binarize the matrix
if (binarize) {
  indata@x[indata@x > 0] <- 1
}

# format peak info
peakinfo <- read.table(peaks) 
peakinfo <- peakinfo[, 1:4]
names(peakinfo) <- c("chr","bp1","bp2","name")
row.names(peakinfo) <- peakinfo$name

# cell metadata
cellinfo <- data.frame(cells = barcodes)
row.names(cellinfo) <- barcodes

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

#Remove peaks in less than n cells
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds) > 0) > min_cells,] 

# monocle3 umap construction
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")

# Create metacells
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

cat('Running cicero ...\n')
# Run cicero
conns <- run_cicero(cicero_cds, chr_sizes) 

# Output gzip file
fwrite(conns, paste0(out, '/cicero_',set,'.tsv.gz'), sep='\t')

#################
cat('Parsing results ...\n')

cc.rna <- as.data.table(conns)
cc.rna <- cc.rna[coaccess >= 0.2]

tss_anno <- fread(anno)
setkey(tss_anno, CREID)
cc.rna$location1 <- tss_anno[cc.rna$Peak1, proximity]
cc.rna$location2 <- tss_anno[cc.rna$Peak2, proximity]
cc.rna[, type := paste(location1, location2, sep='-')]
cc.rna[type == 'distal-proximal', type := 'proximal-distal']
cc.rna[, table(RNA_connections=type) / 2] 
cc.rna[, table(RNA_connections_percent=type)/.N]

fwrite(cc.rna[location1 == 'proximal' & location2 =='distal', 
              .(Peak1, Peak2, coaccess)], 
       paste0(out, '/cicero_',set,'.prox-dist.csv'), 
       row.names = F, col.names = T, sep=',')

#######################

f <- function(file, outfile, desc){
  
  a <- fread(cmd=paste0("tr ',_' '\t' < ", file), skip = 1)
  setnames(a, c('chr.1','st.1','end.1','chr.2','st.2','end.2','score'))
  # maybe filter st.1 < st.2 if not done before to avoid double
  x <- a[chr.1 != 'chrM' & score >= 0.2, .(chr.1, 
                                           floor((st.1+end.1)/2), 
                                           floor((st.2+end.2)/2),
                                           paste(chr.1, st.2, end.1, st.2, end.2, sep='_'),
                                           rescale(score, to = c(0,100)),
                                           '+',
                                           floor((st.1+end.1)/2), 
                                           floor((st.2+end.2)/2),
                                           '.',
                                           2,
                                           '1,1',
                                           paste(0,floor((st.2+end.2)/2) - floor((st.1+end.1)/2), sep=',')
  )]
  x[V2 < V3, V6 := '-']
  
  cat(paste0('track ','name=',desc,' graphType=junctions\n'), file = outfile, append = FALSE)
  fwrite(x, file=outfile, row.names = F, col.names = F, append = T, quote = F, sep='\t')
}

f(paste0(out, '/cicero_',set,'.prox-dist.csv'), 
  paste0(out, '/cicero_',set,'.prox-dist.junc.bed'), 
  'cicero_links')

timestamp()
