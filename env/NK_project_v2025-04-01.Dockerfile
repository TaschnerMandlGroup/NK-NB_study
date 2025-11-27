# RUN AS: podman build --build-arg GITHUB_PAT=XXX -t docker.io/swernig/nk_project:v1.0 -f NK_project_v2025-04-01.Dockerfile .

# pull base image
FROM docker.io/rocker/tidyverse:4.4.1

# who maintains this image
LABEL maintainer="Sara Wernig-Zorc <sara.wernig-zorc@ccri.at>"
LABEL version="v1.0"

# make debian noninteractive
ENV DEBIAN_FRONTEND="noninteractive"

# add GitHUb token for authethication
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

# change permissions so that users can dynamically install more libraries within container:
RUN chmod -R a+rw ${R_HOME}/site-library
RUN chmod -R a+rw ${R_HOME}/library

# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf

# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

# install Python dependencies:
RUN apt-get -y update && apt-get -y install \
  fftw3 \
  libfftw3-dev \
  libglpk40 \
  libxt6 \
  libcurl4-openssl-dev \
  libssl-dev \
  libjq-dev \
  libprotobuf-dev \
  protobuf-compiler \
  make \
  libgeos-dev \
  libudunits2-dev \
  libgdal-dev \
  gdal-bin \
  libproj-dev \
  libv8-dev \
  libbz2-dev \
  libhdf5-dev \
  libbz2-dev \
  python3 \
  python3-venv \   
  python3.10-dev \
  python3.10-venv \
  libncurses5-dev \
  libncursesw5-dev \
  dialog \
  apt-utils \
  libgsl-dev \
  libxml2 \
  libxml2-dev \
  libgmp-dev \
  libglpk-dev \
  nano \
  htop \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install general R packages:
RUN R -e "BiocManager::install(c('ggplot2','markdown', 'sparseMatrixStats', 'edgeR', \
 'apeglm', 'DESeq2', 'patchwork', 'glmGamPoi'))"

RUN R -e "BiocManager::install(c('fgsea', 'hypeR', 'msigdbr', 'HDF5Array','terra',  \
'multtest','Rsamtools','GenomicRanges','dplyr','tidyr', 'DelayedArray', 'DelayedMatrixStats', \
'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor'))"


RUN R -e "BiocManager::install(c('DirichletMultinomial','TFBSTools','biomaRt', 'ComplexHeatmap', \
'JASPAR2020','filelock','ProtGenerics','RcppRoll','shinyBS','shinyjs','DT','ensembldb','RcppRoll', \
'Biobase', 'BiocNeighbors', 'BiocGenerics'))"

RUN R -e "install.packages(c('latex2exp','scico','log4r', \
'leiden','leidenAlg','RSQLite'))"
RUN R -e "install.packages('Cairo', repo='https://RForge.net')"

RUN R -e "remotes::install_github(repo = c('sqjin/CellChat', 'VPetukhov/ggrastr', \ 
'cole-trapnell-lab/monocle3'))"

RUN R -e "remotes::install_github(repo = c('bnprks/BPCells', 'immunogenomics/presto', \ 
'cancerbits/canceRbits', 'cancerbits/DElegate'))"

# install Seurat v5 and dependencies:
RUN R -e "remotes::install_github(repo = 'satijalab/seurat', ref = 'v5.1.0')"
RUN R -e "install.packages('Matrix', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('R.utils')"
RUN R -e "remotes::install_github(repo = 'satijalab/seurat-object', ref = 'v5.0.2')"

# install sceasy (for scVI):
RUN R -e "remotes::install_github('cellgeni/sceasy', ref = '0cfc0e39da4b3ce4abf19d2171daa0e4d2acdd03')"

# install scVI:
RUN apt-get -y update && apt-get -y install \ 
  python3-pip \
  && apt-get clean \
  && rm -rf /tmp/* /var/tmp/*
RUN pip install scanpy scvi-tools leidenalg

RUN R -e "remotes::install_github(repo = c('satijalab/seurat-wrappers', \
 'satijalab/seurat-data','mojaveazure/seurat-disk','stuart-lab/signac','satijalab/azimuth'))"   

RUN R -e "remotes::install_github(repo = c('miccec/yaGST','AntonioDeFalco/SCEVAN'))"
RUN R -e "install.packages(c('reticulate', 'SingleCellExperiment', 'cowplot', 'igraph'))"

RUN R -e "BiocManager::install(c('scater', 'zellkonverter','SCpubr', 'scCustomize', \
'ggVennDiagram', 'miloR', 'SingleR', 'clusterProfiler', 'slingshot', 'ggthemes', \
'statmod', 'ggpubr', 'reshape2', 'TSCAN', 'scater', 'org.Hs.eg.db', 'enrichplot', \
'AnnotationDbi', 'pathview', 'aggregateBioVar'))"

# dependencies
RUN R -e "BiocManager::install(c('sva','genefilter'))"

# CCC tools
RUN R -e "remotes::install_github(repo = 'saeyslab/nichenetr')"
RUN R -e "remotes::install_github(repo = 'saeyslab/multinichenetr')"

# Clean up: unset GITHUB_PAT for security
RUN unset GITHUB_PAT