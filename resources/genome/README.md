## Reference Genome file
SCAFE relies on parsed reference genome files for data processing. Please use ./scripts/download.resources.genome to download the genome corresponding to your datasets. The demo data is on hg19.gencode_v32lift37. Alternatively, you can use ./scripts/tool.cm.prep_genome to build your custom genome.

### download.resources.genome [[top]](#0)<a name="21"></a>
   This script download reference genome data and save in ./resources/genome.

```
 Usage:
   download.resources.genome --genome
   
   --genome <required> [string] name of genome reference, currently available genomes:
                                hg19.gencode_v32lift37
                                hg38.gencode_v32
                                mm10.gencode_vM25
                                TAIR10.AtRTDv2

 Dependencies:
   wget
   tar

 To demo run, cd to SCAFE dir and run:
   ./scripts/download.resources.genome \
   --genome=hg19.gencode_v32lift37
```

### tool.cm.prep\_genome [[top]](#0)<a name="12"></a>
   This tool prepares a reference genome assembly and its gene models for others tools in scafe.

```
 Usage:
   tool.cm.prep_genome [options] --gtf_path --fasta_path --chrom_list_path --mask_bed_path --outputPrefix --outDir
   
   --gtf_path         <required> [string] gtf of the gene models
   --fasta_path       <required> [string] fasta of the genome assembly
   --chrom_list_path  <required> [string] list of <chromosome name><\t><alternative chromosome name> 
                                          e.g. <chr1><\t><1>
                                          chromosome name and alternative chromosome name could be the same
                                          alternative chromosome name is necessary if the cellranger bam
                                          file uses alternative chromosome name that is different from those
                                          in $fasta_path
   --mask_bed_path    <required> [string] a bed file specific the CRE regions. For human or mouse, consider 
                                          using ENCODE CREs. for other species, consider using merged ATAC-seq
                                          from multiple tissues. If ATAC is not available, use the +/- 500nt of 
                                          gene model 5'end.
   --outputPrefix     <required> [string] prefix for the output files (should be name of the genome reference)
   --outDir           <required> [string] directory for the output files (should be resource dir in scafe dir)
   --overwrite        (optional) [yes/no] erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   samtools

 To demo run, cd to SCAFE dir and run:
   ./scripts/tool.cm.prep_genome \
   --overwrite=yes \
   --gtf_path=./demo/input/genome/TAIR10.AtRTDv2.gtf.gz \
   --fasta_path=./demo/input/genome/TAIR10.genome.fa.gz \
   --chrom_list_path=./demo/input/genome/TAIR10.chrom_list.txt \
   --mask_bed_path=./demo/input/genome/TAIR10.ATAC.bed.gz \
   --outputPrefix=TAIR10.AtRTDv2 \
   --outDir=./demo/output/genome/
```
