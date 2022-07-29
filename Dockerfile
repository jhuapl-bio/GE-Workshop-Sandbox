FROM condaforge/mambaforge:4.12.0-2

# Make RUN commands use `bash --login` (always source ~/.bashrc on each RUN)
SHELL ["/bin/bash", "--login", "-c"]

#nstall apt dependencies and update conda
RUN apt-get update --allow-releaseinfo-change && apt-get install git -y \
    && apt-get install -y apt-transport-https ca-certificates wget unzip bzip2 libfontconfig1 
    # && apt-get -qq -y remove curl \
    # && apt-get -qq -y autoremove \
    # && apt-get autoclean \
    # && apt update -y \
    # && apt-get install -y build-essential zlib1g-dev libbz2-dev  libncurses5 liblzma-dev make gcc g++  libopenblas-dev gnuplot

RUN mamba create --name artic --yes python=3.6
RUN source /opt/conda/etc/profile.d/conda.sh && conda activate artic && mamba install -c bioconda -c conda-forge --yes  artic

COPY environment.yml /opt/environment.yml
RUN mamba env create -f /opt/environment.yml

RUN apt-get install -y build-essential wget
RUN apt-get install -y zlib1g-dev libbz2-dev  libncurses5 liblzma-dev libopenblas-dev
# ENV PATH /opt/conda/envs/sandbox/bin:$PATH
# RUN /bin/bash -c "source activate sandbox"


# ENTRYPOINT ["mamba", "run", "--no-capture-output", "-n", "sandbox"]