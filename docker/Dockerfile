FROM debian:stretch

RUN echo "deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/" >>/etc/apt/sources.list \
 && apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y perl git wget \
 && apt-get install -y --allow-unauthenticated r-base-core \
 && Rscript -e 'install.packages("ROCR")' \
 && Rscript -e 'install.packages("PRROC")' \
 && Rscript -e 'install.packages("caret")' \
 && Rscript -e 'install.packages("e1071")' \
 && Rscript -e 'install.packages("ggplot2")' \
 && Rscript -e 'install.packages("scales")' \
 && Rscript -e 'install.packages("reshape2")'

RUN apt-get -y install libcurl4-gnutls-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libudunits2-dev libgdal-dev

RUN Rscript -e 'install.packages("Matrix", repos="http://R-Forge.R-project.org")' \
 && Rscript -e 'install.packages("R.utils")' \
 && Rscript -e 'install.packages("data.table")' \
 && Rscript -e 'install.packages("docopt")'
 
RUN Rscript -e 'install.packages("devtools")' \
 && Rscript -e 'install.packages("BiocManager")' \
 && Rscript -e 'devtools::install_version("RcppAnnoy", "0.0.16")' \
 && Rscript -e 'install.packages("BiocManager")' \
 && Rscript -e 'BiocManager::install(c("BiocGenerics", "DelayedArray", \
    "DelayedMatrixStats","limma", "S4Vectors", "SingleCellExperiment", \
    "SummarizedExperiment", "batchelor", "Matrix.utils"))' \
 && Rscript -e 'devtools::install_github("cole-trapnell-lab/leidenbase")' \
 && Rscript -e 'devtools::install_github("cole-trapnell-lab/monocle3")' \
 && Rscript -e 'BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))' \
 && Rscript -e 'devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")'

RUN git clone -v https://github.com/chung-lab/SCAFE 2>2021032201

WORKDIR ./SCAFE

RUN chmod 755 -R ./scripts/ \
 && chmod 755 -R ./resources/bin/ \
 && echo "export PATH=\$PATH:/SCAFE/scripts" >>~/.bashrc
