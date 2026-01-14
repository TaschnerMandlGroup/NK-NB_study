# RUN AS:
# docker build --build-arg GITHUB_PAT=XXX -t docker.io/swernig/nk_project:v1.1 -f NK_project_v2026-01-14.Dockerfile .

FROM docker.io/rocker/tidyverse:4.4.2

LABEL maintainer="Sara Wernig-Zorc <sara.wernig-zorc@ccri.at>"
LABEL version="v1.1"

# Ensure we run as root during setup 
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
    python3-pip \
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
    wget \
    bzip2 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Miniconda
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniconda.sh \
    && ${CONDA_DIR}/bin/conda clean -afy \
    && ln -s ${CONDA_DIR}/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". ${CONDA_DIR}/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate base" >> ~/.bashrc

# Initialize conda for bash
RUN ${CONDA_DIR}/bin/conda init bash

# Install MACS2 via conda
RUN ${CONDA_DIR}/bin/conda install -c bioconda -y macs2 \
    && ${CONDA_DIR}/bin/conda clean -afy

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

# Install ArchR dependencies and ArchR
RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    BiocManager::install(c( \
        'GenomicRanges', \
        'SummarizedExperiment', \
        'GenomeInfoDb', \
        'BiocGenerics', \
        'S4Vectors', \
        'IRanges', \
        'matrixStats', \
        'Biostrings', \
        'BiocParallel', \
        'rtracklayer', \
        'BSgenome', \
        'BSgenome.Hsapiens.UCSC.hg38', \
        'chromVAR', \
        'motifmatchr' \
    ))"

RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"

# Configure reticulate to use conda environment
RUN echo 'Sys.setenv(PATH = paste("/opt/conda/bin", Sys.getenv("PATH"), sep = ":"))' >> ${R_HOME}/etc/Rprofile.site

# Clear GitHub token from environment for security
ENV GITHUB_PAT=""

# Switch back to rstudio user for runtime
WORKDIR /home/rstudio