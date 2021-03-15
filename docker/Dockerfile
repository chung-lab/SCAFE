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
 && Rscript -e 'install.packages("reshape2")' \
 && git clone https://github.com/chung-lab/SCAFE
