FROM jupyter/base-notebook

ARG rstudio_version="1.2.5019"
ARG r_cran="cran40"
USER root

# general OS packages
RUN apt-get update && apt-get install -yq --no-install-recommends \
    # basic packages 
    htop vim emacs unzip git wget axel tmux nano rsync lsb-release sshpass \
    # sshfs dependencies
    fuse libfuse2 sshfs s3fs \
    # build dependencies for singularity
    build-essential uuid-dev libgpgme-dev squashfs-tools libseccomp-dev \
    pkg-config gcc g++ cryptsetup-bin \
    # base dependencies for general stuff
    libhdf5-dev hdf5-tools libigraph0-dev \
    # R dependecies
    libgsl0-dev libxml2-dev libboost-all-dev libssl-dev libhdf5-dev unzip curl libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    build-essential xorg-dev libreadline-dev libc6-dev libclang-8-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libcairo2-dev \
    libpango1.0-dev tcl-dev tk-dev openjdk-8-jdk gfortran \ 
    # RStudio dependencies
    gdebi-core ant clang cmake debsigs dpkg-sig expect fakeroot gnupg1 libacl1-dev libattr1-dev libbz2-dev libcap-dev libclang-6.0-dev \
    libclang-dev libcurl4-openssl-dev libegl1-mesa libfuse2 libgl1-mesa-dev libgtk-3-0 libpam-dev libpango1.0-dev libpq-dev libsqlite3-dev \
    libssl-dev libuser1-dev libxslt1-dev lsof ninja-build patchelf pkg-config psmisc rrdtool software-properties-common uuid-dev zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# conda packages
RUN conda install -c conda-forge mamba
RUN mamba install -y \
      jupyter-server-proxy \
      nbresuse \ 
      qtconsole \
      ipywidgets \
      ipykernel \
      -c conda-forge

RUN jupyter labextension install @jupyterlab/server-proxy

# clean conda cache
RUN  mamba clean --index-cache --tarballs --yes

# fix permissions
RUN fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# rclone
RUN RCLONE_DEB=rclone-current-linux-amd64.deb && \
    cd /tmp && \
    wget --quiet https://downloads.rclone.org/$RCLONE_DEB && \
    apt-get install -y /tmp/$RCLONE_DEB && \
    rm /tmp/$RCLONE_DEB

# go
RUN export VERSION=1.15 OS=linux ARCH=amd64 && \
    cd /tmp && \
    wget --quiet https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

ENV PATH=/usr/local/go/bin:$PATH

# singularity
RUN export VERSION=3.6.1 && \
    cd /tmp && \
    wget --quiet https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    rm singularity-${VERSION}.tar.gz && \
    cd /tmp/singularity && \
    ./mconfig && make -C builddir && make -C builddir install && \
    rm -rf /tmp/singularity/

# install R
# https://cran.r-project.org/bin/linux/ubuntu/README.html
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -c | awk '{print $2}')-$r_cran/" | sudo tee -a /etc/apt/sources.list && \
    add-apt-repository ppa:c2d4u.team/c2d4u4.0+ && \
    apt-get update && apt-get install -yq --no-install-recommends \
        r-base \
        r-base-dev \
    && apt-get clean \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    rm -rf /var/lib/apt/lists/*

# install RStudio
RUN RSTUDIO_PKG=rstudio-server-${rstudio_version}-amd64.deb && \
    cd /tmp && \
    wget -q https://download2.rstudio.org/server/bionic/amd64/${RSTUDIO_PKG} && \
    gdebi -n /tmp/${RSTUDIO_PKG} && \
    rm /tmp/${RSTUDIO_PKG}
ENV PATH="${PATH}:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="/usr/lib/R/lib:/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-8-openjdk-amd64/jre/lib/amd64/server"


# pip packages
RUN pip install --upgrade --no-cache \
      jupyter-rsession-proxy \ 
      rpy2

# Use this URL to work with the binary packages available 
# https://packagemanager.rstudio.com/client/#/repos/1/overview
RUN echo 'options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"))' > ~/.Rprofile

## R packages
RUN Rscript -e 'install.packages(c("IRkernel", "Rmagic", "rJava", "BiocManager", "ggplot2", "devtools", "Seurat"), dependencies = TRUE)'
RUN Rscript -e 'IRkernel::installspec(prefix="/opt/conda")'
RUN Rscript -e 'BiocManager::install()'
RUN Rscript -e 'devtools::install_github("jalvesaq/colorout")'

# clean conda cache
RUN  mamba clean --index-cache --tarballs --yes

# fix permissions
RUN fix-permissions /opt/conda/share/jupyter/kernels/ && \
    fix-permissions /usr/lib/R/ && \
    fix-permissions /usr/local/lib/R/site-library

# give jovyan sudo permissions
RUN sed -i -e "s/Defaults    requiretty.*/ #Defaults    requiretty/g" /etc/sudoers && \
    echo "jovyan ALL= (ALL) NOPASSWD: ALL" >> /etc/sudoers.d/jovyan

RUN mkdir -p /sanger/ && \
    conda list > /sanger/conda.info && \
    lsb_release -a > /sanger/ubuntu.info && \
    singularity --version > /sanger/singularity.info && \
    rclone --version > /sanger/rclone.info && \
    R --version > /sanger/r.info && \
    echo "$rstudio_version" > /sanger/rstudio.info
