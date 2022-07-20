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
RUN conda create --name consensus && conda activate consensus && conda install --yes -c bioconda artic

COPY environment.yml /opt/environment.yml
RUN conda env create  -f /opt/environment.yml



ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "sandbox"]