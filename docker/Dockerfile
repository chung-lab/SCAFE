FROM debian:stretch
RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y \
    perl r-base git
Run git clone https://github.com/chung-lab/SCAFE
Run Rscript -e 'install.packages("ROCR", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("PRROC", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("caret", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("e1071", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("scales", repos = "http://cran.us.r-project.org")' \
 && Rscript -e 'install.packages("reshape2", repos = "http://cran.us.r-project.org")'
