FROM continuumio/miniconda3:4.9.2

# Make RUN commands use `bash --login` (always source ~/.bashrc on each RUN)
SHELL ["/bin/bash", "--login", "-c"]

# install apt dependencies and update conda
RUN apt-get update --allow-releaseinfo-change && apt-get install git -y \
    && apt-get install -y apt-transport-https ca-certificates wget unzip bzip2 libfontconfig1 \
    && update-ca-certificates \
    && apt-get -qq -y remove curl \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && apt update -y \
    && apt-get install -y build-essential zlib1g-dev libbz2-dev  libncurses5 liblzma-dev make gcc g++  libopenblas-dev gnuplot

RUN conda config --add channels biobuilds; conda config --add channels bioconda
# COPY environment_artic.yml /opt/environment_artic.yml
# RUN conda env create -f /opt/environment_artic.yml

# RUN conda create --name medaka --yes python=3.6 && conda clean --all --yes
# RUN conda activate medaka && conda install -c bioconda --yes medaka=1.0.3
RUN git clone --recurse-submodules https://github.com/artic-network/artic-ncov2019 \
    && rm -rf artic-ncov2019/.git 
RUN conda env create -f artic-ncov2019/environment.yml 

RUN wget --no-check-certificate https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_5.1.15_linux64.tar.gz \
    && tar -xzf ont-guppy-cpu_5.1.15_linux64.tar.gz \
    && rm ont-guppy-cpu_5.1.15_linux64.tar.gz

COPY environment.yml /opt/environment.yml
RUN conda env create  -f /opt/environment.yml


ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "sandbox"]