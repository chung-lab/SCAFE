## Demo data and test runs
Please use scafe.download.demo.input to download data and scafe.demo.test.run to test run. Please also make sure genome reference hg19.gencode_v32lift37 is download to ./resources/genome/ using scafe.download.resources.genome

### download.demo.input [[top]](#0)<a name="22"></a>
   This scripts download demo data and save in ./demo/input dir.

```
 Usage:
   scafe.download.demo.input

 Dependencies:
   wget
   tar

 To demo run, cd to SCAFE dir and run:
```

### demo.test.run [[top]](#0)<a name="23"></a>
   This scripts test run for demo data in the ./demo/input dir. It runs all six workflows.
   Demo input data must be downloaded from using ./script/download.demo.input
   Genome reference hg19.gencode_v32lift37 must be downloaded using scafe.download.resources.genome

```
 Usage:
   scafe.demo.test.run [options] --run_outDir
   
   --run_outDir           <required> [string]  directory for the output test runs
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

### download.resources.genome [[top]](#0)<a name="21"></a>
   This script download reference genome data and save in ./resources/genome.

```
 Usage:
   scafe.download.resources.genome --genome
   
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
