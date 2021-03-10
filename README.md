
 ```shell
                        O~~~AA      O~~         O~       O~~~~~~~AO~~~~~~~~A
                      O~~    O~~ O~~   O~~     O~O~~     O~~      O~~       
                       O~~      O~~           O~  O~~    O~~      O~~       
                         O~~    O~~          O~~   O~~   O~~~~~AA O~~~~~~A  
                            O~~ O~~         O~~~~~A O~~  O~~      O~~       
                      O~~    O~~ O~~   O~~ O~~       O~~ O~~      O~~       
                        O~~~~A     O~~~   O~~         O~~O~~      O~~~~~~~AA
 ```

**SCAFE** (Single Cell Analysis of Five' End) is a tool suite for processing of single cell 5’end RNA-seq data. It takes a read alignment file (*.bam) from single-cell 5'end RNA sequencing, precisely identifies the read 5'ends and removes strand invasion artefacts, performs TSS clustering and filters for genuine TSS clusters using logistic regression, defines transcribed cis-regulatory elements (tCRE) and annotated them to gene models. It counts the UMI in tCRE in single cells and returns a tCRE UMI/cellbarcode matrix ready for downstream analyses. 

## Citing SCAFE

Profiling of transcribed cis-regulatory elements in single cells. _bioRxiv_, 2021, [XXXXXX](https://XXXXXXXXXX/)

## What does SCAFE do? It defines tCRE using RNA 5'ends
<div style="text-align:center"><img src="img/tCRE_definition.png?" width="640"></div>
<div style="text-align:center"><img src="img/zenbu.png?" width="400"></div>

Profiling of cis-regulatory elements (CREs, mostly promoters and enhancers) in single cells allows us to interrogate the cell-type specific contexts of gene regulation and genetic predisposition to diseases. Single-cell RNA-5’end-sequencing (sc-end5-seq) methods can detect transcribed CREs (tCREs), enabling the quantification of promoter and enhancer activities in single cells. To identify genuine tCREs, we implemented a workflow which effectively eliminates false positives from sc-end5-seq data.


## Core Tools and Workflows
<div style="text-align:center"><img src="img/flowchart.png?" width="860"></div>

SCAFE consists of a set of perl programs for processing of single cell 5’end RNA-seq data. Major tools are listed here. SCAFE accepts read alignment in .bam format from standard 10X GenomicsTM tool cellranger. Tool bam_to_ctss extracts the 5’ position of reads, taking the 5’ unencoded-Gs into account. Tool remove_strand_invader removes read 5’ends that are strand invasion artifacts by aligning the TS oligo sequence to the immediate upstream sequence of the read 5’end. Tool cluster performs clustering of read 5’ends using 3rd-party tool Paraclu. Tool filter extracts the properties of TSS clusters and performs multiple logistic regression to distinguish genuine TSS clusters from artifacts. Tool annotate define tCREs by merging closely located TSS clusters and annotate tCREs based on their proximity to known genes. Tool count counts the number of UMI within each tCRE in single cells and generates a tCRE-Cell UMI count matrix. SCAFE tools were also implemented workflows for processing of individual samples or pooling of multiple samples.

## Installation and Dependencies

Blablabla

## Getting Ref Genome data
Blablabla

 ```shell
run demo
 ```
## Getting Demo data
Blablabla

 ```shell
run demo
 ```

## Usage of SCAFE
Blablabla

 ```shell
run demo
 ```

## Run SCAFE

```shell
```
