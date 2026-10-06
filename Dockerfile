# Dockerfile -------------------------------------------------------------
# One environment for two uses:
#   1. Binder  : https://mybinder.org builds this file and opens RStudio in the browser
#   2. Local   : docker build -t airway-deseq2 .
#
# rocker/binder = R 4.4.2 + RStudio Server + Jupyter, pinned to a dated CRAN snapshot.
FROM rocker/binder:4.4.2

USER root

# System libraries needed by tidyverse, Bioconductor and plotting packages
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages as root so they land in the shared site library.
# Copied on its own first so Docker can cache this slow layer.
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R

# Copy the project into the home directory (required by Binder)
USER ${NB_USER}
WORKDIR ${HOME}
COPY --chown=${NB_USER}:${NB_USER} . ${HOME}
