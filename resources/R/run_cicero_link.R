#!/usr/bin/env Rscript
suppressMessages({
  library(docopt)
})

doc <- "Usage: run_cicero_link.R [options]

options:
--chr_size_TSV CHR   File with chromsome sizes
--cutoff CUTOFF 		cutoff for CCAN
--cicero_cds_RDS CDS R object of cicero cds
--chr_list CHR       comma delimited list of chr to run, use all to run all [default: all]
--out OUT            Output folder
--help               Show this help text"

opt <- docopt(doc) # docopt parsing
#print(opt)

timestamp()
options(scipen=999)
cat('Loading packages ...\n')

suppressMessages({
  library(cicero)
  library(data.table)
})


#cicero_cds_RDS <- '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/test_run/workflow.sc.solo/pbmc_stim/link/pbmc_stim/cicero/cicero_cds.RDS'
#chr_size_TSV <- '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/resources/genome/hg19.gencode_v32lift37/tsv/chrom.sizes.tsv'
#chr_list <- 'chr18'
#chr_list <- 'all'
#out <- '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/test_run/workflow.sc.solo/pbmc_stim/link/pbmc_stim/cicero/test/'

cutoff <- as.numeric(opt$cutoff)
cicero_cds_RDS <- opt$cicero_cds_RDS
chr_size_TSV <- opt$chr_size_TSV
chr_list <- opt$chr_list
out <- opt$out

if(!dir.exists(out)) dir.create(out, recursive=TRUE)
sample_n = 100

cat('Loading Data ...\n')
chr_size <- read.table(chr_size_TSV)
cicero_cds <- readRDS(cicero_cds_RDS)

if (chr_list != 'all') {
	chr_list <- unlist(strsplit(chr_list, ","))
	chr_size <- chr_size[chr_size$V1 %in% chr_list ,]
}

cat('Running cicero ...\n')
# Run cicero
conns <- run_cicero(cicero_cds, chr_size, sample_num = sample_n) 

CCAN_assigns <- generate_ccans(conns, coaccess_cutoff_override=cutoff)
# Output gzip file
fwrite(conns, paste0(out, '/cicero_link.tsv.gz'), sep='\t')
fwrite(CCAN_assigns, paste0(out, '/cicero_cis_coactivity_network.tsv.gz'), sep='\t')
