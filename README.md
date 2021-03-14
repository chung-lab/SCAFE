<h1 align="center"> <i>SCAFE</i> (Single Cell Analysis of Five'Ends)</h1>

 ```shell
           5'-O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AAA-3'
                        O~~~AA      O~~         O~       O~~~~~~~AO~~~~~~~~A
                      O~~    O~~ O~~   O~~     O~O~~     O~~      O~~       
                       O~~      O~~           O~  O~~    O~~      O~~       
                         O~~    O~~          O~~   O~~   O~~~~~AA O~~~~~~A  
                            O~~ O~~         O~~~~~A O~~  O~~      O~~       
                      O~~    O~~ O~~   O~~ O~~       O~~ O~~      O~~       
                        O~~~~A     O~~~   O~~         O~~O~~      O~~~~~~~AA
       ┌─ᐅ 5'-O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-3'
 ...===┴========================================================================================...
 ```

*SCAFE* (Single Cell Analysis of Five'Ends) provides an end-to-end solution for processing of single cell 5’end RNA-seq data. It takes a read alignment file (**.bam*) from single-cell RNA-5’end-sequencing (e.g. 10xGenomics Chromimum®), precisely maps the cDNA 5'ends (i.e. transcription start sites, TSS), filters for the artefacts and identifies genuine TSS clusters using logistic regression. Based on the TSS clusters, it defines transcribed cis-regulatory elements (tCRE) and annotated them to gene models. It then counts the UMI in tCRE in single cells and returns a tCRE UMI/cellbarcode matrix ready for downstream analyses, e.g. cell-type clustering, linking promoters to enhancers
*etc* .

## Citing *SCAFE*

Profiling of transcribed cis-regulatory elements in single cells. *bioRxiv*, 2021, [XXXXXX](https://XXXXXXXXXX/)

## What does *SCAFE* do?
<div style="text-align:center"><img src=".github/images/tCRE_definition.png?" width="860"></div>

### *SCAFE* extracts transcribed cis-regulatory elements from single-cell RNA-5’end-sequencing data
Profiling of cis-regulatory elements (CREs, mostly promoters and enhancers) in single cells allows us to interrogate the cell-type specific contexts of gene regulation and genetic predisposition to diseases. Single-cell RNA-5’end-sequencing methods (sc-end5-seq, available from [10xGenomics Chromimum®](https://kb.10xgenomics.com/hc/en-us/articles/360000939852-What-is-the-difference-between-Single-Cell-3-and-5-Gene-Expression-libraries-)) theorectically captures the 5'end of cDNA, which represents transcription start sites (TSS). Measuring the RNA output at TSS allows us to precisely locate transcribed CREs (tCREs) on the genome, enabling the quantification of promoter and enhancer activities in single cells. **Figure (a)** shows the sc-end5-seq signal at the two promoters of gene *DHX30*. It highlights the consistency between sc-end5-seq and sc-ATAC-seq data, as well as the dynamic alternaitve TSS usage between cell states (i.e. resting and stimulated immune cells). *SCAFE* identify genuine TSS information from sc-end5-seq data. **Figure (b)** illustrates the stretagies of SCAFE to defined tCRE from TSS information. **Figure (c)** shows the proximal and distal tCREs defined by *SCAFE* at *PTGER4* locus. **Figure (d)** shows the activities of the proximal and distal tCREs at *PTGER4* locus in resting and simulated immune cells, demonstarting the capability of *SCAFE* to study CRE activities in single cells.

## How does *SCAFE* do it?
<div style="text-align:center"><img src=".github/images/flowchart.png?" width="860"></div>

### *SCAFE* Core Tools and Workflows
*SCAFE* consists of [a set of perl programs](scripts/) for processing of sc-end5-seq data. Major tools are listed in **Figure (a)**. *SCAFE* accepts read alignment in *\*.bam* format from 10xGenomics Chromimum® software *cellranger*. Tool ***bam\_to\_ctss*** extracts the 5’ position of reads, taking the 5’ unencoded-Gs into account. Tool ***remove\_strand\_invader*** removes read 5’ends that are strand invasion artifacts by aligning the TS oligo sequence to the immediate upstream sequence of the read 5’end. Tool *cluster* performs clustering of read 5’ends using 3rd-party tool *Paraclu*. Tool ***filter*** extracts the properties of TSS clusters and performs multiple logistic regression to distinguish genuine TSS clusters from artifacts. Tool ***annotate*** define tCREs by merging closely located TSS clusters and annotate tCREs based on their proximity to known genes. Tool ***count*** counts the number of UMI within each tCRE in single cells and generates a tCRE-Cell UMI count matrix. SCAFE tools were also implemented workflows for processing of individual samples or pooling of multiple samples.

*P.S.* *SCAFE* also accepts bulk 5'end RNA-Seq data (e.g. bulk CAGE). See [scripts](scripts/) for details.

### *SCAFE* discovers *de novo* genunie TSS clusters and tCREs
A fraction of TSS identified based on read 5′ends from template switching (TS) reactions (used in 10xGenomics Chromimum®) may not be genuine, attributed to various artefacts including strand invasion and other sources. This results in excessive artifactual TSS misidentified along the gene body, collectively known as “exon painting”. While strand invasion artefacts can be specifically minimized by considering the complementarity of TSS upstream sequence to TS oligo sequence, a non-negligible fraction of artefactual TSS remains after filtering for strand invasion. To minimize the artifactual TSS, *SCAFE* examines the properties of TSS clusters, as shown in **Figure (b)**, and devised a classifier to identify genuine TSS based on multiple logistic regression. This classifier, i.e. logistic probability, achieved excellent performance with AUC>0.98 across sequencing depths and outperformed all individual metrics. This is implemented in the tool ***filter***.

## Dependencies
### perl
*SCAFE* is mainly written in perl (v5.24.1 or later). All scripts are standalone applications and **DO NOT require installations** of extra perl modules. Check whether perl is properly installed on your system.

```shell
#--- Check your perl version
perl --version
```

### R
*SCAFE* relies on R for logistic regression, ROC analysis and graph plotting. Rscript **(v3.5.1 or later)** and the following R packages have to be properly installed:

* [ROCR](https://cran.r-project.org/web/packages/ROCR/readme/README.html), [PRROC](https://cran.r-project.org/web/packages/PRROC/index.html), [caret](https://cran.r-project.org/web/packages/caret/index.html), [e1071](https://cran.r-project.org/web/packages/e1071/index.html), [ggplot2](https://ggplot2.tidyverse.org/), [scales](https://cran.r-project.org/web/packages/scales/index.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

```shell
#--- Check your Rscript version, must be >3.5.1
Rscript --version

#--- Check your R packages, install if missing
Rscript -e 'if (!require("ROCR")) install.packages("ROCR", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("PRROC")) install.packages("PRROC", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("caret")) install.packages("caret", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("e1071")) install.packages("e1071", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("ggplot2")) install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("scales")) install.packages("scales", repos = "http://cran.us.r-project.org")'
Rscript -e 'if (!require("reshape2")) install.packages("reshape2", repos = "http://cran.us.r-project.org")'
```
### Other 3rd party applications
*SCAFE* also relies on a number of 3rd party applications. The binaries and executables (Linux) of these applications are distributed with this reprository in the directory ./resources/bin and **DO NOT require installations**.

* [bigWigAverageOverBed](https://github.com/ENCODE-DCC/kentUtils), [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils),  [bedtools](https://bedtools.readthedocs.io/en/latest/), [samtools](http://www.htslib.org/), [paraclu](http://cbrc3.cbrc.jp/~martin/paraclu/), [paraclu-cut.sh](http://cbrc3.cbrc.jp/~martin/paraclu/)

### OS
SCAFE was developed and tested on Debian GNU/Linux 9. Running SACFE on other OS are not guranteed.

## Installing SCAFE
Once you ensured the above dependencies are met, you are ready to download SCAFE to your system.

### Clone this respository

```shell
#--- make a directory to install SCAFE
mkdir -pm 755 /my/path/to/install/
cd /my/path/to/install/

#--- Obtain SCAFE from github
git clone https://github.com/chung-lab/SCAFE
cd SCAFE

#--- making sure the scripts and binaries are executable
chmod 755 -R ./scripts/
chmod 755 -R ./resources/bin/
```
### Check the dependencies
SCAFE depends on perl, R and a number of 3rd party tools (as listed above). To ensure the dependencies, please run ./scripts/check.dependencies

```shell
#--- run check.dependencies to check 
./scripts/check.dependencies
```
If everything runs smoothly, you should see the following report on screen.

```shell

===============================
start checking
===============================
Check Type           Check Item                          Check Status
===============      ===============                     ===============
perl executables     workflow.sc.subsample               successful
perl executables     workflow.sc.solo                    successful
perl executables     workflow.sc.pool                    successful
perl executables     workflow.bk.subsample               successful
perl executables     workflow.bk.solo                    successful
perl executables     workflow.bk.pool                    successful
perl executables     tool.sc.subsample_ctss              successful
perl executables     tool.sc.pool                        successful
perl executables     tool.sc.count                       successful
perl executables     tool.sc.bam_to_ctss                 successful
perl executables     tool.cm.remove_strand_invader       successful
perl executables     tool.cm.prep_genome                 successful
perl executables     tool.cm.filter                      successful
perl executables     tool.cm.ctss_to_bigwig              successful
perl executables     tool.cm.cluster                     successful
perl executables     tool.cm.annotate                    successful
perl executables     tool.bk.subsample_ctss              successful
perl executables     tool.bk.pool                        successful
perl executables     tool.bk.count                       successful
perl executables     tool.bk.bam_to_ctss                 successful
perl executables     download.resources.genome           successful
perl executables     download.demo.input                 successful
perl executables     demo.test.run                       successful
perl executables     check.install.dependencies          successful
dropbox access       wget tar data                       successful
3rd-party apps       bedtools                            successful
3rd-party apps       samtools                            successful
3rd-party apps       paraclu                             successful
3rd-party apps       paraclu-cut.sh                      successful
3rd-party apps       bedGraphToBigWig                    successful
3rd-party apps       bigWigAverageOverBed                successful
R version            v3.5 or later                       successful
R packages           ROCR                                successful
R packages           PRROC                               successful
R packages           caret                               successful
R packages           e1071                               successful
R packages           ggplot2                             successful
R packages           scales                              successful
R packages           reshape2                            successful
R scripts            benchmark_roc.R                     successful
R scripts            build_glm.R                         successful
R scripts            predict_prob.R                      successful
===============================
Finished checking
===============================

Successful for all checks. SCAFE should run well.

```

## Getting started with demo data
Now you have enssured all dependencies and downloaded SCAFE, time to get the demo data and test a few runs on the demo data.

### Download demo data and reference genome
Demo data and reference genome must be downloaded for testing *SCAFE* on your system.

```shell
#--- download the demo data using script download.demo.input
./scripts/download.demo.input
	
#--- download the reference genome hg19.gencode_v32lift37 for testing demo data
./scripts/download.resources.genome --genome=hg19.gencode_v32lift37
```

### Test run a single cell solo workflow for demo data
Now, let's test *SCAFE* with a workflow (*workflow.sc.solo*) that processes one library of single cell 5'end RNA-seq data. First we check out the help message of *workflow.sc.solo*. Remember you can always check the help message for all scripts of *SCAFE* using the *--help* flag. 

```shell
#--- check out the help message of workflow.sc.solo
./scripts/workflow.sc.solo --help
```
It should print the help message as the followings:

```shell
          5'-O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AAA-3'
                       O~~~AA      O~~         O~       O~~~~~~~AO~~~~~~~~A
                     O~~    O~~ O~~   O~~     O~O~~     O~~      O~~       
                      O~~      O~~           O~  O~~    O~~      O~~       
                        O~~    O~~          O~~   O~~   O~~~~~AA O~~~~~~A  
                           O~~ O~~         O~~~~~A O~~  O~~      O~~       
                     O~~    O~~ O~~   O~~ O~~       O~~ O~~      O~~       
                       O~~~~A     O~~~   O~~         O~~O~~      O~~~~~~~AA
      ┌─ᐅ 5'-O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-3'
...===┴========================================================================================...

                     Single Cell Analysis of Five'End (SCAFE) Tool Suite 
                               ---> workflow.sc.solo <---
                 <--- workflow, single-cell mode, process a single sample --->

Description:
  This workflow process a single sample, from a cellranger bam file to tCRE UMI/cellbarcode count matrix

Usage:
  workflow.sc.solo [options] --run_bam_path --run_cellbarcode_path --genome --run_tag --run_outDir

  --run_bam_path         <required> [string]  bam file from cellranger, can be read 1 only or pair-end
  --run_cellbarcode_path <required> [string]  tsv file contains a list of cell barcodes,
                                              barcodes.tsv.gz from cellranger
  --genome               <required> [string]  name of genome reference, e.g. hg19.gencode_v32lift37
  --run_tag              <required> [string]  prefix for the output files
  --run_outDir           <required> [string]  directory for the output files
  --training_signal_path (optional) [string]  quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                              regions (e.g. annotated CRE, in bed format) used for training of logical 
                                              regression model If null, $usr_glm_model_path must be supplied for 
                                              pre-built logical regression model. It overrides usr_glm_model_path 
                                              (default=null)
  --testing_signal_path  (optional) [string]  quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                              regions (e.g. annotated CRE, in bed format) used for testing the performance 
                                              of the logical regression model. If null, annotated TSS from $genome will be 
                                              used as binary genomic regions. (default=null)
  --max_thread           (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                              avoid memory overflow (default=5)
  --overwrite            (optional) [yes/no]  erase run_outDir before running (default=no)

Dependencies:
  R packages: 'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales', 'reshape2'
  bigWigAverageOverBed
  bedGraphToBigWig
  bedtools
  samtools
  paraclu
  paraclu-cut.sh

For demo, cd to SCAFE dir and run,
  ./scripts/workflow.sc.solo \
  --overwrite=yes \
  --run_bam_path=./demo/input/sc.solo/demo.cellranger.bam \
  --run_cellbarcode_path=./demo/input/sc.solo/demo.barcodes.tsv.gz \
  --genome=hg19.gencode_v32lift37 \
  --run_tag=demo \
  --run_outDir=./demo/output/sc.solo/

```

The help message details the input options, noted some are ***\<required\>*** and ***(optional)***. *workflow.sc.solo* takes a bam file (*--run\_bam\_path=*), a cellbarcode list file (*--run\_cellbarcode\_path=*), the corresponding reference genome (*--genome=*) and the output prefix (*--run_tag=*) and directory (*--run_outDir=*). At the end of the message, it also prints the commands for running the script on demo data. Now, copy the command and run as the following:

```shell
#--- run the workflow on the demo.cellranger.bam, it'll take a couple of minutes 
./scripts/workflow.sc.solo \
--overwrite=yes \
--run_bam_path=./demo/input/sc.solo/demo.cellranger.bam \
--run_cellbarcode_path=./demo/input/sc.solo/demo.barcodes.tsv.gz \
--genome=hg19.gencode_v32lift37 \
--run_tag=demo \
--run_outDir=./demo/output/sc.solo/
```

If it finishes smoothly, you should see the following message on screen:

```shell
 #=====================================================
 Results of all tasks found. All tasks run successfully
 #=====================================================
```

You can also check out some of the outputs: 

```shell
#--- check the output of tCRE annotation tables
ls -alh ./demo/output/sc.solo/annotate/demo/log

#--- check the output of tCRE UMI count matrix
ls -alh ./demo/output/sc.solo/count/demo/matrix
```

### Test run all workflows on bulk and single cell data
```shell
#--- check out the help message of demo.test.run
./scripts/demo.test.run --help

#--- run the all six available workflows on the demo bulk and single
#--- it'll take a around 20 minutes
./scripts/demo.test.run \
--overwrite=yes \
--run_outDir=./demo/output/
```
If everything runs smoothly, you should see the following report on screen

```shell
#=======================#
Results of Demo Test Run.
#=======================#

workflow                       tool                           status    
==============                 ==============                 ==============
workflow.bk.pool               manager                        successful
workflow.bk.pool               tool.bk.pool                   successful
workflow.bk.pool               tool.cm.cluster                successful
workflow.bk.pool               tool.cm.filter                 successful
workflow.bk.pool               tool.cm.ctss_to_bigwig         successful
workflow.bk.pool               tool.cm.annotate               successful
workflow.bk.pool               tool.bk.count                  successful
workflow.bk.pool               tool.bk.count                  successful
workflow.bk.solo               manager                        successful
workflow.bk.solo               tool.bk.bam_to_ctss            successful
workflow.bk.solo               tool.cm.cluster                successful
workflow.bk.solo               tool.cm.filter                 successful
workflow.bk.solo               tool.cm.ctss_to_bigwig         successful
workflow.bk.solo               tool.cm.annotate               successful
workflow.bk.solo               tool.bk.count                  successful
workflow.bk.subsample          manager                        successful
workflow.bk.subsample          tool.bk.subsample_ctss         successful
workflow.bk.subsample          tool.cm.cluster                successful
workflow.bk.subsample          tool.cm.filter                 successful
workflow.bk.subsample          tool.cm.ctss_to_bigwig         successful
workflow.bk.subsample          tool.cm.annotate               successful
workflow.bk.subsample          tool.bk.count                  successful
workflow.sc.pool               manager                        successful
workflow.sc.pool               tool.sc.pool                   successful
workflow.sc.pool               tool.cm.remove_strand_invader  successful
workflow.sc.pool               tool.cm.cluster                successful
workflow.sc.pool               tool.cm.filter                 successful
workflow.sc.pool               tool.cm.ctss_to_bigwig         successful
workflow.sc.pool               tool.cm.annotate               successful
workflow.sc.pool               tool.sc.count                  successful
workflow.sc.pool               tool.sc.count                  successful
workflow.sc.solo               manager                        successful
workflow.sc.solo               tool.sc.bam_to_ctss            successful
workflow.sc.solo               tool.cm.remove_strand_invader  successful
workflow.sc.solo               tool.cm.cluster                successful
workflow.sc.solo               tool.cm.filter                 successful
workflow.sc.solo               tool.cm.ctss_to_bigwig         successful
workflow.sc.solo               tool.cm.annotate               successful
workflow.sc.solo               tool.sc.count                  successful
workflow.sc.subsample          manager                        successful
workflow.sc.subsample          tool.sc.subsample_ctss         successful
workflow.sc.subsample          tool.cm.remove_strand_invader  successful
workflow.sc.subsample          tool.cm.cluster                successful
workflow.sc.subsample          tool.cm.filter                 successful
workflow.sc.subsample          tool.cm.ctss_to_bigwig         successful
workflow.sc.subsample          tool.cm.annotate               successful
workflow.sc.subsample          tool.sc.count                  successful
```

## Run *SCAFE* on your own data 

### Running *SCAFE* workflows with default options

We recommend most users to run *SCAFE* on their own data using workflows with default options. There are 3 types of workflow: **(1)** "solo" for processing of a single library, **(2)** "pool" for pooling of multiple libraries and **(3)** "subsample" for down-sampling a single library.


The basic input for *SCAFE* is a *\*.bam* file and a list of cellbarcode. Its 

```shell
run demo

```
### Running *SCAFE* individual tools with custom options 
Some users might want to run *SCAFE* using individual tools with default options. 

### Making a custom reference genome
Currently, four reference genomes ara available. See *./script/download.resources.genome* for downloading. Alternatively, some users might work on genomes of other organisms, or prefer to use custom gene models for annotating tCREs.  *tool.cm.prep_genome* converts user-supplied genome *\*.fasta* and gene model *\*.gtf* into necessary files for *SCAFE*. You can check out the help message for inputs of *tool.cm.prep_genome* and then test run a demo using TAIR10 genome with AtRTDv2 gene model.

```shell
#--- check out the help message of tool.cm.prep_genome
./scripts/tool.cm.prep_genome --help

#--- run the tool on the TAIR10 assembly with gene model AtRTDv2 
./scripts/tool.cm.prep_genome \
--overwrite=yes \
--gtf_path=./demo/input/genome/TAIR10.AtRTDv2.gtf.gz \
--fasta_path=./demo/input/genome/TAIR10.genome.fa.gz \
--chrom_list_path=./demo/input/genome/TAIR10.chrom_list.txt \
--mask_bed_path=./demo/input/genome/TAIR10.ATAC.bed.gz \
--outputPrefix=TAIR10.AtRTDv2 \
--outDir=./demo/output/genome/
```
### Running *SCAFE* with bulk CAGE data 









