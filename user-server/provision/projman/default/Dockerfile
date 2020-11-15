################################################################################
# build using build.sh script to ensure files are copied to the right location before building 
################################################################################

# Jupyer - RStudio Server 1.2.5019 - R 4.0.3 - Seurat 3.1.2
FROM quay.io/cellgeni/projman:latest

LABEL maintainer="my4@sanger.ac.uk"
LABEL name="MYBase"

RUN sudo apt-get update
RUN sudo apt-get install -y \
        htop \
        tmux

#Create the symlink to home.  Assumes volumes mounted from higher level.
RUN ln -s /nfs/users/nfs_m/my4/ /home/jovyan/trueHome

# R Packages
RUN Rscript -e 'devtools::install_github("jalvesaq/colorout")'
#RUN Rscript -e 'BiocManager::install("multtest")'
#RUN Rscript -e 'devtools::install_version(package = "Seurat", version = package_version("2.3.4"))'

# Copy everything from the init folder into home directory
COPY --chown=jovyan:users files/initJupyter/ /home/jovyan/
#COPY files/.Rprofile /home/jovyan/


