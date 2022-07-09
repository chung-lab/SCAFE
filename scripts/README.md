## SCAFE Tools and Workflows<a name="0"></a>
This folder contains the following tools and workflows. A tool perform a single task and a workflow runs multiple tools. Some scripts are seperately implemented as bulk (.bk.) and single-cell (.sc.) mode, while others are common (.cm.) for both.

* [**scafe.workflow.sc.subsample**](#1) ---> workflow, single-cell mode, subsample ctss
* [**scafe.workflow.sc.solo**](#2) ---> workflow, single-cell mode, process a single sample
* [**scafe.workflow.cm.aggregate**](#3) ---> workflow, commond mode, aggregate ctss of multiple samples to define CRE
* [**scafe.workflow.bk.subsample**](#4) ---> workflow, bulk mode, subsample ctss
* [**scafe.workflow.bk.solo**](#5) ---> workflow, bulk mode, process a single sample
* [**scafe.tool.sc.subsample\_ctss**](#6) ---> tool, single-cell mode, subsample ctss
* [**scafe.tool.sc.count**](#7) ---> tool, single-cell mode, count of UMI within tCRE
* [**scafe.tool.sc.bam\_to\_ctss**](#8) ---> tool, single-cell mode, convert bam to ctss
* [**scafe.tool.cm.remove\_strand\_invader**](#9) ---> tool, common mode, remove strand invader artefact
* [**scafe.tool.cm.prep\_genome**](#10) ---> tool, common mode, prepare custom reference genome
* [**scafe.tool.cm.filter**](#11) ---> tool, common mode, filter for genuine TSS clusters
* [**scafe.tool.cm.directionality**](#12) ---> tool, common mode, calculate directionality of tCREs
* [**scafe.tool.cm.ctss\_to\_bigwig**](#13) ---> tool, common mode, convert ctss to bigwig
* [**scafe.tool.cm.cluster**](#14) ---> tool, common mode, cluster ctss
* [**scafe.tool.cm.annotate**](#15) ---> tool, common mode, define and annotate tCRE
* [**scafe.tool.cm.aggregate**](#16) ---> tool, common mode, aggregate ctss of multiple samples
* [**scafe.tool.bk.subsample\_ctss**](#17) ---> tool, bulk mode, subsample ctss
* [**scafe.tool.bk.count**](#18) ---> tool, bulk mode, count ctss within tCREs
* [**scafe.tool.bk.bam\_to\_ctss**](#19) ---> tool, bulk mode, convert bam to ctss bed
* [**scafe.download.resources.genome**](#20) ---> download, reference genome to resources dir
* [**scafe.download.demo.input**](#21) ---> download, demo input data for testing
* [**scafe.demo.test.run**](#22) ---> demo, run demo data for testing
* [**scafe.check.dependencies**](#23) ---> check dependencies


### scafe.workflow.sc.subsample [[top]](#0)<a name="1"></a>
   This workflow subsamples a ctss file, defines tCRE and generate a tCRE UMI/cellbarcode count matrix
   Subsampling is useful to investigate the effect of sequencing depth to tCRE definition

```
 Usage:
   scafe.workflow.sc.subsample [options] --UMI_CB_ctss_bed_path --run_cellbarcode_path --subsample_num --genome --run_tag --run_outDir
   
   --UMI_CB_ctss_bed_path <required> [string]  ctss file for subsampling, one line one cellbarcode-UMI combination,
                                               *UMI_CB.ctss.bed.gz from scafe.tool.sc.bam_to_ctss.pl, 
                                               4th column cellbarcode-UMI and 5th column is number of unencoded-G
   --run_cellbarcode_path <required> [string]  tsv file contains a list of cell barcodes,
                                               barcodes.tsv.gz from cellranger
   --subsample_num        <required> [integer] number of UMI to be subsampled
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
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.workflow.sc.subsample \
   --overwrite=yes \
   --UMI_CB_ctss_bed_path=./demo/input/sc.subsample/demo.UMI_CB.ctss.bed.gz \
   --run_cellbarcode_path=./demo/input/sc.subsample/demo.barcodes.tsv.gz \
   --subsample_num=100000 \
   --genome=hg19.gencode_v32lift37 \
   --run_tag=demo \
   --run_outDir=./demo/output/sc.subsample/
```

### scafe.workflow.sc.solo [[top]](#0)<a name="2"></a>
   This workflow process a single sample, from a cellranger bam file to tCRE UMI/cellbarcode count matrix

```
 Usage:
   scafe.workflow.sc.solo [options] --run_bam_path --run_cellbarcode_path --genome --run_tag --run_outDir
   
   --run_bam_path                  <required> [string] bam file from cellranger, can be read 1 only or pair-end
   --run_cellbarcode_path          <required> [string] tsv file contains a list of cell barcodes,
                                                       barcodes.tsv.gz from cellranger
   --genome                        <required> [string] name of genome reference, e.g. hg19.gencode_v32lift37
   --run_tag                       <required> [string] prefix for the output files
   --run_outDir                    <required> [string] directory for the output files
   --training_signal_path          (optional) [string] quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                                       regions (e.g. annotated CRE, in bed format) used for training of logical 
                                                       regression model If null, $usr_glm_model_path must be supplied for 
                                                       pre-built logical regression model. It overrides usr_glm_model_path 
                                                       (default=null)
   --testing_signal_path           (optional) [string] quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                                       regions (e.g. annotated CRE, in bed format) used for testing the performance 
                                                       of the logical regression model. If null, annotated TSS from $genome will be 
                                                       used as binary genomic regions. (default=null)
   --detect_TS_oligo (optional) [match|trim|skip|auto] in bam_to_ctss step, the modes of detecting TS oligo. 1. match: search for 
                                                       TS oligo sequence on the read, identify the TSO/cDNA junction as 5'end of 
                                                       the read. This works only when the error rate of the TS oligo region on 
                                                       the read is low, otherwise a considerable number of read will be invalid. 
                                                       2. trim: assuming the 1st N bases of the reads are TS oligo, without 
                                                       checking the actual sequence. N is determined by the length of TS oligo. 
                                                       3. skip: assuming the TS oligo was not sequenced, the 1st base of the read
                                                       will be treated as the 1st base after the TS oligo. 4. auto: automatically 
                                                       determines the best mode, best of the observed error rate of the TS oligo
                                                       and the frequency of 5'end softclipped bases by the aligner. If softcliped 
                                                       bases is close to the length of TS oligo, mode 1 or 2 will be chosen, 
                                                       depending on the observed error rate of the TS oligo (error rate <= 0.1, 
                                                       mode 1 will be chosen or mode 2 otherwise). If softcliped base os close to 
                                                       zero, mode 3 will be chosen. (default=auto).
   --max_thread                   (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                                       avoid memory overflow (default=5)
   --overwrite                     (optional) [yes/no] erase run_outDir before running (default=no)

 Dependencies:
   R packages: 'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales', 'reshape2'
   bigWigAverageOverBed
   bedGraphToBigWig
   bedtools
   samtools
   paraclu
   paraclu-cut.sh
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.workflow.sc.solo \
   --overwrite=yes \
   --run_bam_path=./demo/input/sc.solo/demo.cellranger.bam \
   --run_cellbarcode_path=./demo/input/sc.solo/demo.barcodes.tsv.gz \
   --genome=hg19.gencode_v32lift37 \
   --run_tag=demo \
   --run_outDir=./demo/output/sc.solo/
```

### scafe.workflow.cm.aggregate [[top]](#0)<a name="3"></a>
   This workflow process a multiple samples from scafe.workflow.sc.solo or scafe.workflow.bk.solo to define CRE using aggregated signal

```
 Usage:
   scafe.workflow.cm.aggregate [options] --lib_list_path --genome --run_tag --run_outDir
   
   --lib_list_path        <required> [string] a list of libraries, in formation of 
                                              <lib_ID><\t><collapse_ctss><\t><unencoded_G_collapse_ctss>
                                              lib_ID = Unique ID of the lib
                                              collapse_ctss = *.collapse.ctss.bed.gz from scafe.tool.sc.bam_to_ctss or scafe.tool.bk.bam_to_ctss
                                                              or *.pass.ctss.bed.gz from scafe.tool.cm.remove_strand_invader
                                              unencoded_G_collapse_ctss = *unencoded_G.collapse.ctss.bed.gz from 
                                                                          scafe.tool.sc.bam_to_ctss or scafe.tool.bk.bam_to_ctss
   --genome               <required> [string] name of genome reference, e.g. hg19.gencode_v32lift37
   --run_tag              <required> [string] prefix for the output files
   --run_outDir           <required> [string] directory for the output files
   --training_signal_path (optional) [string] quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                              regions (e.g. annotated CRE, in bed format) used for training of logical 
                                              regression model If null, $usr_glm_model_path must be supplied for 
                                              pre-built logical regression model. It overrides usr_glm_model_path 
                                              (default=null)
   --testing_signal_path (optional) [string]  quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                              regions (e.g. annotated CRE, in bed format) used for testing the performance 
                                              of the logical regression model. If null, annotated TSS from $genome will be 
                                              used as binary genomic regions. (default=null)
   --max_thread          (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                              avoid memory overflow (default=5)
   --overwrite           (optional) [yes/no]  erase run_outDir before running (default=no)

 Dependencies:
   R packages: 'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales', 'reshape2'
   bigWigAverageOverBed
   bedGraphToBigWig
   bedtools
   samtools
   paraclu
   paraclu-cut.sh
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.workflow.cm.aggregate \
   --overwrite=yes \
   --lib_list_path=./demo/input/cm.aggregate/lib_list_path.txt \
   --genome=hg19.gencode_v32lift37 \
   --run_tag=demo \
   --run_outDir=./demo/output/cm.aggregate/
```

### scafe.workflow.bk.subsample [[top]](#0)<a name="4"></a>
   This workflow subsamples a ctss file, defines tCRE and generate tCRE read count
   Subsampling is useful to investigate the effect of sequencing depth to tCRE definition

```
 Usage:
   scafe.workflow.bk.subsample [options] --long_ctss_bed_path --subsample_num --genome --run_tag --run_outDir
   
   --long_ctss_bed_path    <required> [string] ctss file for subsampling, one line one read
                                               *long.ctss.bed.gz from scafe.tool.bk.bam_to_ctss.pl, 
   --subsample_num        <required> [integer] number of UMI to be subsampled
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

 To demo run, cd to SCAFE dir and run:
   scafe.workflow.bk.subsample \
   --overwrite=yes \
   --long_ctss_bed_path=./demo/input/bk.subsample/demo.long.ctss.bed.gz \
   --subsample_num=100000 \
   --genome=hg19.gencode_v32lift37 \
   --run_tag=demo \
   --run_outDir=./demo/output/bk.subsample/
```

### scafe.workflow.bk.solo [[top]](#0)<a name="5"></a>
   This workflow process a single sample, from a bulk CAGE bam file to read count per tCRE 

```
 Usage:
   scafe.workflow.bk.solo [options] --run_bam_path --genome --run_tag --run_outDir
   
   --run_bam_path         <required> [string]  bam file (of CAGE reads), can be read 1 only or pair-end
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
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.workflow.bk.solo \
   --overwrite=yes \
   --run_bam_path=./demo/input/bk.solo/demo.CAGE.bam \
   --genome=hg19.gencode_v32lift37 \
   --run_tag=demo \
   --run_outDir=./demo/output/bk.solo/
```

### scafe.tool.sc.subsample\_ctss [[top]](#0)<a name="6"></a>
   This tool subsample a ctss bed file and maintains the cellbarcode and UMI information

```
 Usage:
   scafe.tool.sc.subsample_ctss [options] --UMI_CB_ctss_bed_path --subsample_num --outputPrefix --outDir
   
   --UMI_CB_ctss_bed_path <required> [string]  ctss file for subsampling, one line one cellbarcode-UMI combination,
                                               *UMI_CB.ctss.bed.gz from scafe.tool.sc.bam_to_ctss, 
                                               4th column cellbarcode-UMI and 5th column is number of unencoded-G
   --subsample_num        <required> [integer] number of UMI to be subsampled
   --outputPrefix         <required> [string]  prefix for the output files
   --outDir               <required> [string]  directory for the output files
   --overwrite            (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.tool.sc.subsample_ctss \
   --overwrite=yes \
   --UMI_CB_ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.UMI_CB.ctss.bed.gz \
   --subsample_num=100000 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.subsample/subsample_ctss/
```

### scafe.tool.sc.count [[top]](#0)<a name="7"></a>
   This tool counts the UMI within a set of user-defined regions, e.g. tCRE, and returns a UMI/cellbarcode matrix

```
 Usage:
   scafe.tool.sc.count [options] --countRegion_bed_path --cellBarcode_list_path --ctss_bed_path --outputPrefix --outDir
   
   --countRegion_bed_path   <required> [string] bed file contains the regions for counting CTSS, e.g. tCRE ranges, 
                                                *.CRE.coord.bed.gz from scafe.tool.cm.annotate.pl
   --cellBarcode_list_path  <required> [string] tsv file contains a list of cell barcodes,
                                                barcodes.tsv.gz from cellranger
   --ctss_bed_path          <required> [string] ctss file for counting,
                                                *CB.ctss.bed.gz from scafe.tool.sc.bam_to_ctss.pl, 
                                                4th column cellbarcode and 5th column is number UMI
   --genome                 <required> [string] name of genome reference, e.g. hg19.gencode_v32lift37
   --ctss_scope_bed_path    <optional> [string] bed file contains the regions for filtering CTSS, e.g. tssCluster ranges, 
                                                so only the ctss within these ranges (i.e. scope) will be count. This is to 
                                                prevent over permissive counting to ctss in the CRE range by stricting only 
                                                ctss within valid tssClusters to be counted. 
                                                *.tssCluster.default.filtered.bed.gz from scafe.tool.cm.filter.
                                                It will skip filtering if not file was provide (default=null).
   --ctss_scope_slop_bp    <optional> [integer] the length of the boundary extension (in bp) for filter region provided in
                                                option ctss_scope_bed_path. All regions in ctss_scope_bed_path will be 
                                                extended both side by ctss_scope_slop_bp, for controlling the permissiveness
                                                of the counting. A large value of ctss_scope_slop_bp (e.g. 400) will be 
                                                equivalent to no filtering (default=0).
   --outputPrefix           <required> [string] prefix for the output files
   --outDir                 <required> [string] directory for the output files
   --overwrite              (optional) [yes/no] erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools

 To demo run, cd to SCAFE dir and run:
   scafe.tool.sc.count \
   --overwrite=yes \
   --genome=hg19.gencode_v32lift37 \
   --countRegion_bed_path=./demo/output/sc.solo/annotate/demo/bed/demo.CRE.annot.bed.gz \
   --cellBarcode_list_path=./demo/input/sc.solo/demo.barcodes.tsv.gz \
   --ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.CB.ctss.bed.gz \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/count/
```

### scafe.tool.sc.bam\_to\_ctss [[top]](#0)<a name="8"></a>
   This tool converts a bam file to a ctss bed file, identifies read 5'end (capped TSS, i.e. ctss),
   extracts the unencoded G information, pileup ctss, and deduplicate the UMI

```
 Usage:
   scafe.tool.sc.bam_to_ctss [options] --bamPath --genome --outputPrefix --outDir
   
   --bamPath                       <required> [string] bam file from cellranger, can be read 1 only or pair-end
   --genome                        <required> [string] name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix                  <required> [string] prefix for the output files
   --outDir                        <required> [string] directory for the output files
   --include_flag                  (optional) [string] samflag to be included, comma delimited 
                                                       e.g. '64' to include read1, (default=null)
   --exclude_flag                  (optional) [string] samflag to be excluded, comma delimited, 
                                                       e.g. '128,256,4' to exclude read2, secondary alignment 
                                                       and unaligned reads (default=128,256,4)
   --min_MAPQ                     (optional) [integer] minimum MAPQ to include (default=0)
   --max_thread                   (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                                       avoid memory overflow (default=5)
   --TS_oligo_seq                  (optional) [string] Template switching oligo sequence for identification of 
                                                       5'end (default=TTTCTTATATGGG) 
   --detect_TS_oligo (optional) [match/trim/skip/auto] in bam_to_ctss step, the modes of detecting TS oligo. 1. match: search for 
                                                       TS oligo sequence on the read, identify the TSO/cDNA junction as 5'end of 
                                                       the read. This works only when the error rate of the TS oligo region on 
                                                       the read is low, otherwise a considerable number of read will be invalid. 
                                                       2. trim: assuming the 1st N bases of the reads are TS oligo, without 
                                                       checking the actual sequence. N is determined by the length of TS oligo. 
                                                       3. skip: assuming the TS oligo was not sequenced, the 1st base of the read
                                                       will be treated as the 1st base after the TS oligo. 4. auto: automatically 
                                                       determines the best mode, best of the observed error rate of the TS oligo
                                                       and the frequency of 5'end softclipped bases by the aligner. If softcliped 
                                                       bases is close to the length of TS oligo, mode 1 or 2 will be chosen, 
                                                       depending on the observed error rate of the TS oligo (error rate <= 0.1, 
                                                       mode 1 will be chosen or mode 2 otherwise). If softcliped base os close to 
                                                       zero, mode 3 will be chosen. (default=auto).
   --overwrite                    (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   samtools
   tabix
   bgzip
   

 To demo run, cd to SCAFE dir and run:
   scafe.tool.sc.bam_to_ctss \
   --overwrite=yes \
   --bamPath=./demo/input/sc.solo/demo.cellranger.bam \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/bam_to_ctss/
```

### scafe.tool.cm.remove\_strand\_invader [[top]](#0)<a name="9"></a>
   This tool identify and remove strand invader artefact from a ctss bed file, 
   by aligning the sequence immediate upstream of a ctss to TS oligo sequence

```
 Usage:
   scafe.tool.cm.remove_strand_invader [options] --ctss_bed_path --genome --outputPrefix --outDir
   
   --ctss_bed_path      <required> [string]  "collapse" ctss file from scafe.tool.sc.bam_to_ctss.pl, 
                                             4th column is number of cells and 5th column is number UMI
   --genome             <required> [string]  name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix       <required> [string]  prefix for the output files
   --outDir             <required> [string]  directory for the output files
   --min_edit_distance  (optional) [integer] edit distance threshold to define strand invader 
                                             the smaller value, the more stringent defintion of strand invader
                                             (default=5)
   --min_end_non_G_num  (optional) [integer] immediate upstream non-G number threshold to define strand invader
                                             the smaller value, the more stringent defintion of strand invader
                                             (default=2)
   --max_thread         (optional) [integer] maximum number of parallel threads, capped at 
                                             10 to avoid memory overflow (default=5)
   --TS_oligo_seq       (optional) [string]  Template switching oligo sequence for identification 
                                             of 5'end (default=TTTCTTATATGGG) 
   --overwrite          (optional) [yes/no]  [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.remove_strand_invader \
   --overwrite=yes \
   --ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/remove_strand_invader/
```

### scafe.tool.cm.prep\_genome [[top]](#0)<a name="10"></a>
   This tool prepares a reference genome assembly and its gene models for others tools in scafe.

```
 Usage:
   scafe.tool.cm.prep_genome [options] --gtf_path --fasta_path --chrom_list_path --mask_bed_path --outputPrefix --outDir
   
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
   scafe.tool.cm.prep_genome \
   --overwrite=yes \
   --gtf_path=./demo/input/genome/TAIR10.AtRTDv2.gtf.gz \
   --fasta_path=./demo/input/genome/TAIR10.genome.fa.gz \
   --chrom_list_path=./demo/input/genome/TAIR10.chrom_list.txt \
   --mask_bed_path=./demo/input/genome/TAIR10.ATAC.bed.gz \
   --outputPrefix=TAIR10.AtRTDv2 \
   --outDir=./demo/output/genome/
```

### scafe.tool.cm.filter [[top]](#0)<a name="11"></a>

```
 Usage:
   scafe.tool.cm.filter [options] --ctss_bed_path --ung_ctss_bed_path --tssCluster_bed_path --genome --outputPrefix --outDir
   
   --ctss_bed_path           <required> [string]  ctss file contains all ctss,
                                                  *.collapse.ctss.bed.gz from scafe.tool.sc.bam_to_ctss.pl, 
                                                  5th column is number reads/UMI
   --ung_ctss_bed_path       <required> [string]  ctss file contains only ctss with unencoded G,
                                                  *.unencoded_G.collapse.ctss.bed.gz from scafe.tool.sc.bam_to_ctss.pl, 
                                                  5th column is number reads/UMI
   --tssCluster_bed_path     <required> [string]  bed file contains all TSS clusters,
                                                  *.tssCluster.bed.gz from scafe.tool.cm.cluster.pl
   --genome                  <required> [string]  name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix            <required> [string]  prefix for the output files
   --outDir                  <required> [string]  directory for the output files
   --tssCluster_flank_size   (optional) [integer] size of regions (each side) flanking a TSS cluster summit for
                                                  counting UMI/reads for expression levels calculation (default = 75)
   --local_bkgd_extend_size  (optional) [integer] size of regions (each side) flanking a TSS cluster summit for 
                                                  defining the scope for calculating local background (default = 500)
   --min_gold_num            (optional) [integer] minimum number of gold standard regions for training and testing the
                                                  logical regression model (default = 100)
   --training_pct            (optional) [float]   top and bottom percentage of the TSS clusters, ranked by signal in 
                                                  $training_signal_path, used for training of logical regression model
                                                  (default = 5)
   --training_signal_path    (optional) [string]  quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                                  regions (e.g. annotated CRE, in bed format) used for training of logical 
                                                  regression model If null, $usr_glm_model_path must be supplied for 
                                                  pre-built logical regression model. It overrides usr_glm_model_path 
                                                  (default=null)
   --testing_signal_path     (optional) [string]  quantitative signal (e.g. ATAC -logP, in bigwig format), or binary genomic 
                                                  regions (e.g. annotated CRE, in bed format) used for testing the performance 
                                                  of the logical regression model. If null, annotated TSS from $genome will be 
                                                  used as binary genomic regions. (default=null)
   --usr_glm_model_path      (optional) [string]  pre-built logical regression model from the Caret package in R. Used only if 
                                                  training_signal_path is not supplied. Models were pre-built for each genome
                                                  and used as default.
   --Rscript_bin             (optional) [string]  path to the Rscript bin, aim to allow users to supply an R version other the 
                                                  system wide R version. Package Caret must be installed. (Defaul = Rscript)
   --default_cutoff          (optional) [integer] logistic probablity cutoffs for the "default" stringency (Default = 0.5)
   --exclude_chrom_list      (optional) [string]  a list of comma delimited chromosome to be excluded in the training and 
                                                  testing of the logical regression model (Default = chrM)
   --overwrite               (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   R packages: 'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales', 'reshape2'
   bedtools
   bigWigAverageOverBed

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.filter \
   --overwrite=yes \
   --ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --ung_ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.unencoded_G.collapse.ctss.bed.gz \
   --tssCluster_bed_path=./demo/output/sc.solo/cluster/demo/bed/demo.tssCluster.bed.gz \
   --training_signal_path=./demo/input/atac/demo.atac.bw \
   --testing_signal_path=./demo/input/atac/demo.atac.bw \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/filter/
```

### scafe.tool.cm.directionality [[top]](#0)<a name="12"></a>
   This tool counts the ctss within tCREs and calculate their strand bias (i.e. directionality)

```
 Usage:
   scafe.tool.cm.directionality [options] --CRE_bed_path --CRE_info_path --ctss_bed_path --outputPrefix --outDir
   
   --CRE_bed_path           <required> [string] bed file contains the tCRE ranges, 
                                                *.CRE.coord.bed.gz from scafe.tool.cm.annotate
   --CRE_info_path          <required> [string] information table of tCRE, 
                                                *.CRE.info.tsv.gz from scafe.tool.cm.annotate
   --ctss_bed_path          <required> [string] ctss file for counting,
                                                *collapse.ctss.bed.gz from scafe.tool.sc.bam_to_ctss, 
                                                5th column is number of read
   --ctss_scope_bed_path    <optional> [string] bed file contains the regions for filtering CTSS, e.g. tssCluster ranges, 
                                                so only the ctss within these ranges (i.e. scope) will be count. This is to 
                                                prevent over permissive counting to ctss in the CRE range by stricting only 
                                                ctss within valid tssClusters to be counted. 
                                                *.tssCluster.default.filtered.bed.gz from scafe.tool.cm.filter.
                                                It will skip filtering if not file was provide (default=null).
   --outputPrefix           <required> [string] prefix for the output files
   --outDir                 <required> [string] directory for the output files
   --overwrite              (optional) [yes/no] erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.directionality \
   --overwrite=yes \
   --CRE_bed_path=./demo/output/sc.solo/annotate/demo/bed/demo.CRE.coord.bed.gz \
   --CRE_info_path=./demo/output/sc.solo/annotate/demo/log/demo.CRE.info.tsv.gz \
   --ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/directionality/
```

### scafe.tool.cm.ctss\_to\_bigwig [[top]](#0)<a name="13"></a>
   This tool converts a ctss bed file into two bigwig file, one for each strand, for visualization purpose 

```
 Usage:
   scafe.tool.cm.ctss_to_bigwig [options] --ctss_bed_path --genome --outputPrefix --outDir
   
   --ctss_bed_path  <required> [string] "collapse" ctss file from scafe.tool.sc.bam_to_ctss.pl
   --genome         <required> [string] name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix   <required> [string] prefix for the output files
   --outDir         <required> [string] directory for the output files
   --overwrite      (optional) [yes/no] erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedGraphToBigWig

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.ctss_to_bigwig \
   --ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/ctss_to_bigwig/
```

### scafe.tool.cm.cluster [[top]](#0)<a name="14"></a>
   This tool generate TSS cluster from a ctss bed file, using an external tool paraclu with user-defined cutoffs

```
 Usage:
   scafe.tool.cm.cluster [options] --cluster_ctss_bed_path --outputPrefix --outDir
   
   --cluster_ctss_bed_path       <required> [string]  ctss file used for clustering,
                                                      "collapse" ctss file from scafe.tool.sc.bam_to_ctss.pl, 
                                                      4th column is number of cells and 5th column is number UMI
   --outputPrefix                <required> [string]  prefix for the output files
   --outDir                      <required> [string]  directory for the output files
   --count_ctss_bed_path_list    (optional) [string]  comma delimited list of ctss bed file, 
                                                      using for filtering of clusters based signal 
                                                      (default=$cluster_ctss_bed_path) 
   --count_scope_bed_path        (optional) [string]  a bed file specify the scope for counting in $count_ctss_bed_path_list, 
                                                      using for filtering of clusters based signal
                                                      (default=$cluster_ctss_bed_path) 
   --min_pos_count               (optional) [integer] minimum counts per position, used for filtering the raw signal 
                                                      in $cluster_ctss_bed_path before clustering (default = 1)
   --min_cluster_cpm             (optional) [float]   minimum counts per million (cpm) for a cluster (default = 1e-5)
   --min_summit_count            (optional) [integer] minimum counts at the summit of a cluster (default = 3)
   --min_cluster_count           (optional) [integer] minimum counts within a cluster (default = 5)
   --min_num_sample_expr_cluster (optional) [integer] minimum number of samples (or cells) detected at the 
                                                      summit of a cluster (default = 3)
   --min_num_sample_expr_summit  (optional) [integer] minimum number of samples (or cells) detected within 
                                                      of a cluster (default = 5)
   --merge_dist                  (optional) [integer] maximum distance for merging closely located clusters, 
                                                      -1 to turn off merging (default = -1)
   --overwrite                   (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   paraclu
   paraclu-cut.sh
   bedtools

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.cluster \
   --overwrite=yes \
   --cluster_ctss_bed_path=./demo/output/sc.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/cluster/
```

### scafe.tool.cm.annotate [[top]](#0)<a name="15"></a>
   This tool defines tCRE from TSS clusters and annotates them based their overlap with gene models.

```
 Usage:
   scafe.tool.cm.annotate [options] --tssCluster_bed_path --tssCluster_info_path --genome --outputPrefix --outDir
   
   --tssCluster_bed_path       <required> [string]   bed file contains the ranges of filtered TSS clusters,
                                                     *.tssCluster.*.filtered.bed.gz from scafe.tool.cm.filter.pl
   --tssCluster_info_path      <required> [string]   tsv file contains the information of all TSS clusters,
                                                     *.tssCluster.log.tsv from scafe.tool.cm.filter.pl
   --genome                    <required> [string]   name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix              <required> [string]   prefix for the output files
   --outDir                    <required> [string]   directory for the output files
   --up_end5Rng                (optional) [integer]  TSS clusters will be classified as gene TSS, exonic, intron 
                                                     and intergenic. $up_end5Rng determines the range upstream of 
                                                     annotated gene TSS to be used for gene TSS assignment 
                                                     (default = 500)
   --dn_end5Rng                (optional) [integer]  TSS clusters will be classified as gene TSS, exonic, intron 
                                                     and intergenic. $dn_end5Rng determines the range downstream of 
                                                     annotated gene TSS to be used for gene TSS assignment 
                                                     (default = 500)
   --exon_slop_rng             (optional) [integer]  TSS clusters will be classified as gene TSS, exonic, intron 
                                                     and intergenic. $exon_slop_rng determines the range to be extended
                                                     (i.e. slopped) from exon for assignment of exonic class. 
                                                     Used -1 to NOT to extend (default = -1)
   --merge_dist                (optional) [integer]  TSS clusters outside annotated gene promoters are grouped
                                                     as "dummy genes" (for operational uniformity) by merging closely 
                                                     located TSS clusters.  $merge_dist determines the maximum distances 
                                                     between TSS clusters to be merged (default = 500)
   --addon_length              (optional) [integer]  see $merge_dist. add-on "dummy transcrips" will assigned to TSS cluster of 
                                                     "dummy genes" (for operational uniformity).$addon_length determines 
                                                     the length of these add-on "dummy transcrips" (default = 500).
   --proximity_slop_rng        (optional) [integer]  TSS clusters will be assigned to annotated gene TSS are "proximal"
                                                     TSS clusters. $proximity_slop_rng determines the range to be extended
                                                     (i.e. slopped) from gene TSS for assignment of proximal TSS clusters. 
                                                     (default = 500)
   --merge_strandness          (optional) [string]   see $merge_dist. $merge_strandness decides the merge to be 
                                                     strand-aware ("stranded") or strand-agnostic "strandless".
                                                     (default = strandless)
   --proximal_strandness       (optional) [string]   closely located proximal TSS clusters are merged  
                                                     tCREs. $proximal_strandness decides the merge to be 
                                                     strand-aware ("stranded") or strand-agnostic "strandless".
                                                     (default = stranded)
   --CRE_extend_size           (optional) [integer]  tCREs were defined by merging the extended ranges of TSS clusters.
                                                     $CRE_extend_size determine the size of this range (both sides of 
                                                     summit) (default = 500)
   --CRE_extend_upstrm_ratio   (optional) [float]    see $CRE_extend_size. $CRE_extend_upstrm_ratio determines the ratio 
                                                     (X:1) of flanking sizes on the upstream and downstream of summit. 
                                                     e.g. $CRE_extend_upstrm_ratio=4, upstream and downstream size will be 
                                                     taken as 4:1 ratio. $CRE_extend_size=500 and $CRE_extend_upstrm_ratio=4,
                                                     upstream and downstream will be 400 and 100 respectively 
                                                     (default = 4)
   --stitch_distance           (optional) [integer]  distance (nt) for stitching distal tCRE for defining hyperactive distal loci.
                                                     aka superenhancer candidates. If undefined, an optimized value will be 
                                                     determined based on the tangent to rank-distance plot. As a note, the original 
                                                     distance for stitching enhancer in 
                                                     ROSE (http://younglab.wi.mit.edu/super_enhancer_code.html) is 12,500.
                                                     (default = undefined)
   --min_total_exp_frac        (optional) [fraction] minimum fraction of the expression amount (read/UMI) within a gene for an annotated 
                                                     tCRE to be regarded as an unannotated promoter of the gene. The total expression amount 
                                                     of a gene is defined as the total number of UMI/read of all its annotated promoters.
                                                     (default = 0.1)
   --min_gene_strand_read_frac (optional) [fraction] minimum fraction of the expression amount (read/UMI) on the gene strand (in total number 
                                                     of UMI/read both strand) of a tCRE to be regarded as an unannotated promoter of the gene.
                                                     (default = 0.8)
   --min_spreadness            (optional) [fraction] minimum spreadness of the distal CRE to be considered as hyperactive distal loci.
                                                     Spreadness is defined as "num-of-CRE/top-fraction". Top-fraction is defined as the expression 
                                                     amount (read/UMI) on the highest expressed distal CRE in the expression amount of all distal CRE 
                                                     within the locus. Num-of-CRE refers to the number of CRE within the distal CRE locus 
                                                     (default = 4, equivalent to Num-of-CRE=3 and Top-fraction = 0.75, i.e. 3/0.75 = 4)
   --Rscript_bin               (optional) [string]   path to the Rscript bin, aim to allow users to supply an R version other the 
                                                     system wide R version. Package Caret must be installed. (default = Rscript)
   --overwrite                 (optional) [yes/no]   erase outDir/outputPrefix before running (default=no)
 

 Dependencies:
   bedtools,
   R packages: 'ggplot2'

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.annotate \
   --overwrite=yes \
   --tssCluster_bed_path=./demo/output/sc.solo/filter/demo/bed/demo.tssCluster.default.filtered.bed.gz \
   --tssCluster_info_path=./demo/output/sc.solo/filter/demo/log/demo.tssCluster.log.tsv \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/sc.solo/annotate/
```

### scafe.tool.cm.aggregate [[top]](#0)<a name="16"></a>
   This tool aggregate multiple ctss bed files, regardless of single cell or bulk

```
 Usage:
   scafe.tool.cm.aggregate [options] --lib_list_path --genome --outputPrefix --outDir
   
   --lib_list_path  <required> [string] a list of libraries, in formation of 
                                        <lib_ID><\t><collapse_ctss><\t><unencoded_G_collapse_ctss>
                                        lib_ID = Unique ID of the lib
                                        collapse_ctss = *.collapse.ctss.bed.gz from scafe.tool.sc.bam_to_ctss or scafe.tool.bk.bam_to_ctss
                                                        or *.pass.ctss.bed.gz from scafe.tool.cm.remove_strand_invader
                                        unencoded_G_collapse_ctss = *unencoded_G.collapse.ctss.bed.gz from 
                                                                    scafe.tool.sc.bam_to_ctss or scafe.tool.bk.bam_to_ctss
   --genome        <required> [string]  name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix  <required> [string]  prefix for the output files
   --outDir        <required> [string]  directory for the output files
   --max_thread    (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                        avoid memory overflow (default=5)
   --overwrite     (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   tabix
   bgzip
   

 To demo run, cd to SCAFE dir and run:
   scafe.tool.cm.aggregate \
   --overwrite=yes \
   --lib_list_path=./demo/input/cm.aggregate/lib_list_path.txt \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/cm.aggregate/aggregate/
```

### scafe.tool.bk.subsample\_ctss [[top]](#0)<a name="17"></a>
   This tool subsample a ctss bed file from bulk CAGE ctss

```
 Usage:
   scafe.tool.bk.subsample_ctss [options] --UMI_CB_ctss_bed_path --subsample_num --outputPrefix --outDir
   
   --long_ctss_bed_path <required> [string]  ctss file for subsampling, one line read in "long" format,
                                             *long.ctss.bed.gz from scafe.tool.bk.bam_to_ctss.pl, 
   --subsample_num      <required> [integer] number of UMI to be subsampled
   --outputPrefix       <required> [string]  prefix for the output files
   --outDir             <required> [string]  directory for the output files
   --overwrite          (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.tool.bk.subsample_ctss \
   --overwrite=yes \
   --long_ctss_bed_path=./demo/output/bk.solo/bam_to_ctss/demo/bed/demo.long.ctss.bed.gz \
   --subsample_num=100000 \
   --outputPrefix=demo \
   --outDir=./demo/output/bk.subsample/subsample_ctss/
```

### scafe.tool.bk.count [[top]](#0)<a name="18"></a>
   This tool counts the CAGE reads within a set of user-defined regions, e.g. tCRE, and 
   returns the reads per regions

```
 Usage:
   scafe.tool.bk.count [options] --countRegion_bed_path --ctss_bed_path --outputPrefix --outDir
   
   --countRegion_bed_path   <required> [string] bed file contains the regions for counting CTSS, e.g. tCRE ranges, 
                                                *.CRE.coord.bed.gz from scafe.tool.cm.annotate
   --ctss_bed_path          <required> [string] ctss file for counting,
                                                *collapse.ctss.bed.gz from scafe.tool.bk.bam_to_ctss, 
                                                5th column is number of read
   --outputPrefix           <required> [string] prefix for the output files
   --outDir                 <required> [string] directory for the output files
   --overwrite              (optional) [yes/no] erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools

 To demo run, cd to SCAFE dir and run:
   scafe.tool.bk.count \
   --overwrite=yes \
   --countRegion_bed_path=./demo/output/bk.solo/annotate/demo/bed/demo.CRE.annot.bed.gz \
   --ctss_bed_path=./demo/output/bk.solo/bam_to_ctss/demo/bed/demo.collapse.ctss.bed.gz \
   --outputPrefix=demo \
   --outDir=./demo/output/bk.solo/count/
```

### scafe.tool.bk.bam\_to\_ctss [[top]](#0)<a name="19"></a>
   This tool converts a bulk CAGE bam file to a ctss bed file, identifies read 5'end 
   (capped TSS, i.e. ctss), extracts the unencoded G information, pileup ctss, 
   and deduplicate the UMI

```
 Usage:
   scafe.tool.bk.bam_to_ctss [options] --bamPath --genome --outputPrefix --outDir
   
   --bamPath      <required> [string]  bam file (of CAGE reads), can be read 1 only or pair-end
   --genome       <required> [string]  name of genome reference, e.g. hg19.gencode_v32lift37
   --outputPrefix <required> [string]  prefix for the output files
   --outDir       <required> [string]  directory for the output files
   --include_flag (optional) [string]  samflag to be included, comma delimited 
                                       e.g. '64' to include read1, (default=null)
   --exclude_flag (optional) [string]  samflag to be excluded, comma delimited, 
                                       e.g. '128,256,4' to exclude read2, secondary alignment 
                                       and unaligned reads (default=128,256,4)
   --min_MAPQ     (optional) [integer] minimum MAPQ to include (default=0)
   --max_thread   (optional) [integer] maximum number of parallel threads, capped at 10 to 
                                       avoid memory overflow (default=5)
   --overwrite    (optional) [yes/no]  erase outDir/outputPrefix before running (default=no)

 Dependencies:
   bedtools
   samtools
   tabix
   bgzip

 To demo run, cd to SCAFE dir and run:
   scafe.tool.bk.bam_to_ctss \
   --overwrite=yes \
   --bamPath=./demo/input/bk.solo/demo.CAGE.bam \
   --genome=hg19.gencode_v32lift37 \
   --outputPrefix=demo \
   --outDir=./demo/output/bk.solo/bam_to_ctss/
```

### scafe.download.resources.genome [[top]](#0)<a name="20"></a>
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
   scafe.download.resources.genome \
   --genome=hg19.gencode_v32lift37
```

### scafe.download.demo.input [[top]](#0)<a name="21"></a>
   This scripts download demo data and save in ./demo/input dir.

```
 Usage:
   download.demo.input

 Dependencies:
   wget
   tar

 To demo run, cd to SCAFE dir and run:
   scafe.download.demo.input
```

### scafe.demo.test.run [[top]](#0)<a name="22"></a>
   This scripts test run for demo data in the ./demo/input dir. It runs user-selected workflows.
   Demo input data must be downloaded from using ./script/download.demo.input
   Genome reference hg19.gencode_v32lift37 must be downloaded using ./scripts/download.resources.genome

```
 Usage:
   demo.test.run [options] --run_outDir
   
   --run_outDir           <required> [string]  directory for the output test runs
   --workflow             (optional) [string]  comma delimited list of workflows, 
                                               or use 'all' to run all workflows.
                                               Available workflows includes,
                                               scafe.workflow.sc.subsample ---> workflow, single-cell mode, subsample ctss
                                               scafe.workflow.sc.solo ---> workflow, single-cell mode, process a single sample
                                               scafe.workflow.cm.aggregate ---> workflow, common mode, aggregate ctss of multiple samples
                                               scafe.workflow.bk.subsample ---> workflow, bulk mode, subsample ctss
                                               scafe.workflow.bk.solo ---> workflow, bulk mode, process a single sample
                                               (default=all)
   --overwrite            (optional) [yes/no]  erase run_outDir before running (default=no)

 Dependencies:
   R packages: 'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales', 'reshape2'
   bigWigAverageOverBed
   bedGraphToBigWig
   bedtools
   samtools
   paraclu
   paraclu-cut.sh

 To demo run, cd to SCAFE dir and run:
   scafe.demo.test.run \
   --overwrite=yes \
   --run_outDir=./demo/output/
```

### scafe.check.dependencies [[top]](#0)<a name="23"></a>
   This scripts check the integrity of tools and workflow scripts, 3rd executable dependencies and R packages.

```
 Usage:
   check.dependencies

 Dependencies:
   wget
   tar
   Rscript

 To demo run, cd to SCAFE dir and run:
   scafe.check.dependencies
```
