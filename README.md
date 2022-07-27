# Sandbox for JHUAPL/Fogarty Workshops

## Introduction 

This Repo is intended to be run with `docker` on Windows, Mac OS, or most Linux distros. Please see [here](https://www.docker.com/) for setting up Docker on your systems.

## 1a. Project setup pull

```
docker pull jhuaplbio/sandbox
```

## 1b. Project setup build

```
docker build . -t jhuaplbio/sandbox
```

## 2. Test Data Retrieval 

<!-- If you have gdown on your system, download it here with: `gdown 1zrgwheJxhMTvd7zu0fuRhVYYM0aGY5XS -O ./test-data.zip` 

else  -->

Get the test data from [here](https://drive.google.com/file/d/1zrgwheJxhMTvd7zu0fuRhVYYM0aGY5XS/view?usp=sharing). Place it where you have this repository to make it easier to track. 

## 3. Brief Introduction to your Docker environment



This will simply start up your environment in the `sandbox` conda env which contains all of your necessary binaries to run downstream processes

A list of installed binaries for bioinformatics teaching is as follows:

1. [minimap2](https://github.com/lh3/minimap2)  - Alignment
2. [bowtie2](https://github.com/BenLangmead/bowtie2)  - Alignment
3. [trimmomatic](https://github.com/usadellab/Trimmomatic)  - Pre Processing
4. [NanoPlot](https://github.com/wdecoster/NanoPlot)  - PLotting 
5. [FastQC](https://github.com/s-andrews/FastQC)  - Plotting
6. [bamstats](http://bamstats.sourceforge.net/)  - Plotting 
7. [samtools](http://www.htslib.org/)  - Parsing and Variant Calling
8. [bedtools](https://bedtools.readthedocs.io/en/latest/)  - Parsing
9. [medaka](https://github.com/nanoporetech/medaka) - Consensus

## 3a. Running an interactive container 

```
docker container run --rm -it -v $pwd/test-data:/data jhuaplbio/sandbox  bash 
```

To Exit once done, type: `exit`

You've successfully run your first interactive docker container!


If we wanted to start up an interactive shell, we would simply override the entry command (the command run at start) with `bash` like so:


### Windows Powershell

```
docker container run -w /data -v $pwd/test-data:/data --rm -it jhuaplbio/sandbox bash
```

### Unix Terminal

```
docker container run -w /data -v $PWD/test-data:/data --rm -it jhuaplbio/sandbox bash

```


## 3b. Running without interactive container

You can of course run all things without being in the interactive shell, which simulates the linux environment you'd be working in. It would look like this for each and all commands going forwards (e.g. `minimap2`, `samtools`, etc.)

### Windows Powershell

```
docker container run -w /data -v $pwd/test-data:/data --rm -it jhuaplbio/sandbox <command>
```

### Unix Terminal

```
docker container run -w /data -v $PWD/test-data:/data --rm -it jhuaplbio/sandbox <command>
```


:warning:

**All commands going forward will be assuming that we're working in the interactive environment that was just mentioned so make sure to run that before continuing. Your terminal once in should look something like:**

```
(base) root@addcf87af2d3:/# 
```

## 5. Activating Your Conda Environment

Let's take a quick look at the envrionment we are in at startup. It is `base` which does not have any of our dependencies. First, we must activate `sandbox` with 


```
conda activate sandbox
```

Your terminal screen should look like: 

```
(sandbox) root@addcf87af2d3:/# 
```

Let's check our current packages installed in this environment with: `conda list` 

```
(sandbox) root@addcf87af2d3:/# conda list
# packages in environment at /opt/conda/envs/sandbox:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
biopython                 1.79                     pypi_0    pypi
blast                     2.10.1          pl526he19e7b1_3    bioconda
bowtie2                   2.4.1            py38he513fc3_0    bioconda
brotlipy                  0.7.0           py38h27cfd23_1003  
bzip2                     1.0.8                h7b6447c_0  
...
```

If you ever forget the name of your environments, you can run `conda env list` where the `*` indicates which one you're in currently.

```
(sandbox) root@addcf87af2d3:/# conda env list
# conda environments:
#
base                     /opt/conda
consensus                /opt/conda/envs/consensus
sandbox               *  /opt/conda/envs/sandbox
```

## Running Test Data throught the bioinformatic commands

### a. Running Trimming on some Illumina Data 

```
trimmomatic PE viruses/sars_cov_2/ERR6913101_1.fastq.gz \
    viruses/sars_cov_2/ERR6913101_2.fastq.gz \
    viruses_trimmed/ERR6913101_1.trim.fastq.gz \
    viruses_trimmed/ERR6913101_1.untrim.fastq.gz \
    viruses_trimmed/ERR6913101_2.trim.fastq.gz \
    viruses_trimmed/ERR6913101_2.untrim.fastq.gz \
    SLIDINGWINDOW:4:20 \
    MINLEN:25
```

Your output files should be (since this is paired-end reads) as the `1.trim.fastq.gz` and `2.trim.fastq.gz` in `viruses_trimmed`

### b. Running Alignment on Trimmed Illumina Data

First, let's index our FASTA genome file. Since we are aligning against SARS-CoV-2, we will use the reference genome in the test data folder `reference/nCoV-2019.reference.fasta`. We first need to head into the reference folder and run `bowtie2-build`

```
cd reference; bowtie2-build \
    -f /data/reference/nCoV-2019.reference.fasta nCoV-2019; \
    cd /data
```

Next, run alignment against the indexed reference FASTA file

```
mkdir -p alignments;
bowtie2 -x /data/reference/nCoV-2019 \
    -1 viruses_trimmed/ERR6913101_1.trim.fastq.gz \
    -2 viruses_trimmed/ERR6913101_2.trim.fastq.gz \
    -S alignments/ERR6913101_alignments.sam
```

Let's now try minimap2, which is recommended for long reads (Oxford Nanopore runs)

```

mkdir -p /data/alignments/minimap2 && \
minimap2 \
    -x map-ont \
    -a /data/reference/nCoV-2019.reference.fasta \
    -o /data/alignments/minimap2/alignment.sam demux-fastq_pass/NB03.fastq

```



### c. Converting SAM to BAM to reduce filesize

For our short read alignments, let's work on that first...

```
samtools view \
    -S -b alignments/ERR6913101_alignments.sam > alignments/ERR6913101_alignments.bam && \
    rm alignments/ERR6913101_alignments.sam

```

Next, let's take a look at the equivalent with minimap2 for the long reads


```
samtools view \
    -S \
    -b /data/alignments/minimap2/alignment.sam > /data/alignments/minimap2/alignment.bam

```

### d1. Prepping for Variant Calling - Sorting BAM

```
samtools sort alignments/ERR6913101_alignments.bam > alignments/ERR6913101_alignments_sorted.bam
```

### d2. Prepping for Variant Calling - Indexing

```
samtools index alignments/ERR6913101_alignments_sorted.bam 
```

### e. Variant Calling

```
bcftools mpileup \
    -f reference/nCoV-2019.reference.fasta \
    alignments/ERR6913101_alignments_sorted.bam \
        | bcftools call \
            -mv \
            -Ov \
            -o alignments/ERR6913101_alignment.vcf
```


### f.2 Plotting with Bamstats

First, like always, let's focus on the short reads

```
samtools stats alignments/ERR6913101_alignments_sorted.bam > alignments/ERR6913101_alignments_sorted.stats

plot-bamstats \
    -p alignments/plots_bamstats alignments/ERR6913101_alignments_sorted.stats 

```

Next, the long reads doing the same process except on the minimap2 alignment output file (.bam)

```
samtools stats /data/alignments/minimap2/alignment.bam > /data/alignments/minimap2/alignment.stats

plot-bamstats \
    -p alignments/minimap2/plots_bamstats alignments/minimap2/alignment.stats

```


### g. Plotting a sample Oxford Nanopore run with NanoPlot

```
NanoPlot \
    --summary /data/20200514_2000_X3_FAN44250_e97e74b4/sequencing_summary_FAN44250_77d58da2.txt \
    -o nanoplots
```

### h. Plotting a sample Illumina Paired end read set with FastQC

```
mkdir /data/fastqc_plots_fastq;

fastqc \
    viruses_trimmed/ERR6913101_1.trim.fastq.gz viruses_trimmed/ERR6913101_2.trim.fastq.gz \
    -o fastqc_plots

mkdir /data/fastqc_plots_bam;

fastqc \
    -o fastqc_plots_bam \
    -f bam alignments/ERR6913101_alignments_sorted.bam


```


### j. Running Classification with Kraken2

:warning:Requires internet to download the minikraken database. You can also get it from [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)

#### j1. Getting the minikraken database for Kraken2

```
mkdir -p /data/databases

wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz -O /data/databases/minikraken2.tar.gz; 

tar -xvzf --directory /data/databases/minikraken2

```

#### j2. Running Kraken2 using the minikraken2 database (Illumina paired end)

```
mkdir /data/classifications;

kraken2 --db /data/databases/minikraken2 --gzip-compressed --paired --classified-out ERR6913101#.fq viruses_trimmed/ERR6913101_1.trim.fastq.gz viruses_trimmed/ERR6913101_2.trim.fastq.gz --report classifications/ERR6913101.kraken.report
```

#### j3. Running Kraken2 using the minikraken2 database (Nanopore)

```
mkdir /data/classifications;

kraken2 --db /data/databases/minikraken2 --classified-out sample_metagenome#.fq metagenome/sample_metagenome.fastq --report classifications/sample_metagenome.kraken.report
```

#### j4. Prepping the database for Krona Plots

:warning:Requires internet connectivity

```
ktUpdateTaxonomy.sh /data/databases/minikraken2
```

#### j5. Creating a Krona Plot of your Tax calls

```
ktImportTaxonomy -i classifications/ERR6913101.kraken.report -o classifications/ERR6913101.kraken.html -tax /data/databases/minikraken2/

ktImportTaxonomy -i classifications/sample_metagenome.kraken.report -o classifications/sample_metagenome.kraken.html -tax /data/databases/minikraken2/

```

#### k. Creating a consensus genome from Nanopore Reads using Medaka


You may have already noticed that we have a demultiplexed set of files in the `/data/demux-fastq_pass` folder. However, this is not likely to happen when you have a traditional run from MinKNOW. Usually it will spit out a set of fastq files rather than just one single one per sample. To gather them all up, we can run one of the subcommands called `artic gather`. This will push all fastq files into a single one like what we see in `/data/demux-fastq_pass`

```
artic gather --directory /data/20200514_2000_X3_FAN44250_e97e74b4/fastq_pass

```

[Medaka](https://github.com/nanoporetech/medaka) is a good tool for generating consensues and variant calling for most organisms. 

```

conda activate artic-ncov2019 && \
    mkdir -p /data/consensus/medaka && \
    cd /data/consensus/medaka

medaka_consensus \
    -i /data/demux-fastq_pass/NB11.fastq \
    -d /data/reference/nCoV-2019.reference.fasta \
    -o /data/consensus/medaka \
    -m r941_min_high_g360

```

However, if we're working with a select few organisms utilizing a primer scheme, we should opt for the artic bioinformatics pipeline. This pipeline can use both nanopolish and medaka as it's primary variant and consensus calling set of scripts. In addition, it will perform the necessary trimming/pre-processing for you on your demultiplexed fastq files. It also comes with several subcommands that can help prep your data, like concatenated all of your fastq files into a single one to prepare for consensus building, no matter the purpose. 

First, we need to get the primer schemes for SARS-CoV-2 (our organism of interest for this example). These primers will be used for the artic minion process

```

wget --no-check-certificate https://github.com/artic-network/primer-schemes/archive/refs/heads/master.zip && unzip master.zip && mv primer-schemes-master /data/primer_schemes; rm master.zip

```

This is pulling in the set of primer schemes directly from the artic pipeline toolkit's set of primer schemes for EBOLA, SARS-CoV-2, and Nipah

```

artic minion --medaka \
    --medaka-model r941_min_high_g360 \
    --normalise 1000000 \
    --read-file /data/20200514_2000_X3_FAN44250_e97e74b4/fastq_pass \
    --scheme-directory /data/primer_schemes \
    --scheme-version V3 \
    nCoV-2019/3 NB11;

cd /data 

```


## Creating your Own Analysis Runs using 3rd party Docker Images

So now we've used our single sandbox Docker Image. But what if, for example, you have a script that isn't available in this pre-made image? For example, `medaka`, which is as we described previously is available in sandbox but also as a standlone Docker Image from Docker hub [here](https://hub.docker.com/r/staphb/medaka)

Let's run a consensus run, like we did in the consensus pipeline without pre-installing it AND automatically running our command once it is downloaded. 



| Command | Platform |
| ------- | -------- |
| `docker container run -it -v $pwd/test-data:/data -w /data staphb/medaka bash -c "medaka_consensus -i /data/demux-fastq_pass/NB03.fastq -d /data/reference/nCoV-2019.reference.fasta -o /data/output/medaka_consensus -m r941_min_high_g303 -f"` | Windows Powershell |
| `docker container run -it -v $PWD/test-data:/data -w /data staphb/medaka bash -c "medaka_consensus -i /data/demux-fastq_pass/NB03.fastq -d /data/reference/nCoV-2019.reference.fasta -o /data/output/medaka_consensus -m r941_min_high_g303 -f"` | Unix |

Notice something interesting....

We can immediately run a command directly from our Shell, and have all dependencies automatically installed and our command run. This is one of the many strengths of using Docker as you have a fully runnable command with all required dependencies right out of the box!

That's it for the copy+paste section of the Docker portion of this workshop. Now that you've hopefully had a grasp on the inner-workings of Docker, we can move on to making and installation for your own pipeline(s)!


