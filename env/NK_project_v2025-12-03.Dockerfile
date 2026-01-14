# RUN AS:
# docker build --build-arg GITHUB_PAT=XXX -t docker.io/swernig/nk_project:v1.0 -f NK_project_v2025-12-03.Dockerfile .

FROM docker.io/rocker/tidyverse:4.4.2

LABEL maintainer="Sara Wernig-Zorc <sara.wernig-zorc@ccri.at>"
LABEL version="v1.0"

# Ensure we run as root during setup (rocker base already does this)
USER root

# Noninteractive apt
ENV DEBIAN_FRONTEND=noninteractive

# Optional GitHub token for private/ratelimited installs
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

# System dependencies (R + Python + libs)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libfftw3-dev \
    libglpk40 \
    libglpk-dev \
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
    python3 \
    python3-dev \
    python3-venv \
    libncurses5-dev \
    libncursesw5-dev \
    dialog \
    apt-utils \
    libgsl-dev \
    libxml2 \
    libxml2-dev \
    libgmp-dev \
    nano \
    htop \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Make sure rstudio user owns its home and libraries
RUN chown -R rstudio:rstudio /home/rstudio \
    && chown -R rstudio:rstudio ${R_HOME}/site-library \
    && chown -R rstudio:rstudio ${R_HOME}/library

# Keep RStudio sessions alive
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf \
    && echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf \
    && echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

# Install R packages
# Use Ncpus for parallel install to speed up the build
RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    BiocManager::install(c( \
        'ggplot2','markdown','sparseMatrixStats','edgeR','apeglm','DESeq2', \
        'patchwork','glmGamPoi', \
        'fgsea','hypeR','msigdbr','HDF5Array','terra','multtest','Rsamtools', \
        'GenomicRanges','dplyr','tidyr','DelayedArray','DelayedMatrixStats', \
        'limma','lme4','S4Vectors','SingleCellExperiment','SummarizedExperiment', \
        'batchelor', \
        'DirichletMultinomial','TFBSTools','biomaRt','ComplexHeatmap', \
        'JASPAR2020','filelock','ProtGenerics','RcppRoll','shinyBS','shinyjs', \
        'DT','ensembldb', \
        'Biobase','BiocNeighbors','BiocGenerics', \
        'scater','zellkonverter','SCpubr','scCustomize','ggVennDiagram','miloR', \
        'SingleR','clusterProfiler','slingshot','ggthemes','statmod','ggpubr', \
        'reshape2','TSCAN','org.Hs.eg.db','enrichplot','AnnotationDbi', \
        'pathview','aggregateBioVar', \
        'sva','genefilter' \
    ))"

RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    install.packages(c( \
        'latex2exp','scico','log4r','leiden','leidenAlg','RSQLite', \
        'Matrix','R.utils','reticulate','cowplot','igraph' \
    ), repos = 'https://cloud.r-project.org'); \
    install.packages('Cairo', repos = 'https://RForge.net')"

# GitHub packages (split in a couple of calls to reduce failure impact)
RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    remotes::install_github(c( \
    'sqjin/CellChat', \
    'VPetukhov/ggrastr', \
    'cole-trapnell-lab/monocle3' \
    ))"

RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    remotes::install_github(c( \
        'bnprks/BPCells', \
        'immunogenomics/presto', \
        'cancerbits/canceRbits', \
        'cancerbits/DElegate'))"

# Seurat v5 stack
RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    remotes::install_github('satijalab/seurat', ref = 'v5.2.0'); \
    remotes::install_github('satijalab/seurat-object', ref = 'v5.2.0'); \
    remotes::install_github(c( \
        'satijalab/seurat-wrappers', \
        'satijalab/seurat-data', \
        'mojaveazure/seurat-disk', \
        'stuart-lab/signac', \
        'satijalab/azimuth'))"

# CCC tools
RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    remotes::install_github(c( \
        'saeyslab/nichenetr', \
        'saeyslab/multinichenetr', \
        'miccec/yaGST', \
        'AntonioDeFalco/SCEVAN'))"

# Clear GitHub token from environment for security
ENV GITHUB_PAT=""

# Switch back to rstudio user for runtime
WORKDIR /home/rstudio