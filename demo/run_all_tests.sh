#!/bin/bash
cwd=`dirname $0`
cd $cwd

#---[2021/03/22 20:30] cd to SCAFE
cd ..

scafe.check.dependencies
scafe.download.demo.input
scafe.download.resources.genome --genome=hg19.gencode_v32lift37

scafe.tool.sc.link \
--overwrite=yes \
--max_thread=10 \
--CRE_bed_path=./demo/input/sc.link/demo.CRE.coord.bed.gz \
--CRE_info_path=./demo/input/sc.link/demo.CRE.info.tsv.gz \
--count_dir=./demo/input/sc.link/matrix/ \
--genome=hg19.gencode_v32lift37 \
--outputPrefix=demo \
--outDir=./demo/output/sc.link/

scafe.demo.test.run \
--run_outDir=./demo/output/
