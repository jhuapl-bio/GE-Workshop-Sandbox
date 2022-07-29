# Sandbox for JHUAPL/Fogarty Workshops

## Introduction 

This Repo is intended to be run with `docker` on Windows, Mac OS, or most Linux distros. Please see [here](https://www.docker.com/) for setting up Docker on your systems.

## 1a. Project setup pull

```
docker pull jhuaplbio/sandbox
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
10. [artic](https://github.com/nanoporetech/medaka) - Consensus

## 3a. Running an interactive container 


If we wanted to start up an interactive shell, we would simply override the entry command (the command run at start) with `bash` like so:

To download and run our primary sandbox container do

### Windows Powershell

```
docker container run -w /data -v $pwd/test-data:/data --rm -it jhuaplbio/sandbox bash
```

### Unix Terminal

```
docker container run -w /data -v $PWD/test-data:/data --rm -it jhuaplbio/sandbox bash
```

To Exit a contauner once you're done, type: `exit`

You've successfully run your first interactive docker container!

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


(base) root@addcf87af2d3:/# 


## 5. Activating Your Conda Environment

Let's take a quick look at the envrionment we are in at startup. It is `base` which does not have any of our dependencies. First, we must activate `sandbox` with 


```
conda activate sandbox
```

Your terminal screen should look like: 


(sandbox) root@addcf87af2d3:/# Notice the (sandbox) on the left side


## Running Test Data throught the bioinformatic commands

### a. Running Trimming on some Illumina Data 

```
trimmomatic PE viruses/sars_cov_2/ERR6913101_1.fastq.gz \
    viruses/sars_cov_2/ERR6913101_2.fastq.gz \
    viruses_trimmed/ERR6913101_1.trim.fastq.gz \
    viruses_trimmed/ERR6913101_1.untrim.fastq.gz \
    viruses_trimmed/ERR6913101_2.trim.fastq.gz \
    viruses_trimmed/ERR6913101_2.untrim.fastq.gz \
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
    -o /data/alignments/minimap2/NB03.sam demux-fastq_pass/NB03.fastq

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
    -b /data/alignments/minimap2/NB03.sam > /data/alignments/minimap2/NB03.bam

```

### d1. Prepping for Variant Calling - Sorting BAM


With Bowtie2 Short Reads output

```
samtools sort alignments/ERR6913101_alignments.bam > alignments/ERR6913101_alignments_sorted.bam
```

With Minimap2 Long Reads output

```
samtools sort /data/alignments/minimap2/NB03.bam > /data/alignments/minimap2/NB03_sorted.bam
```

### d2. Prepping for Variant Calling - Indexing

```
samtools index alignments/ERR6913101_alignments_sorted.bam 
```

```
samtools index alignments/minimap2/NB03_sorted.bam 
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

Again, lets do the same for our long read example fastq file

```
bcftools mpileup \
    -f reference/nCoV-2019.reference.fasta \
    alignments/minimap2/NB03_sorted.bam \
        | bcftools call \
            -mv \
            -Ov \
            -o alignments/minimap2/NB03_alignment.vcf
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
samtools stats /data/alignments/minimap2/NB03.bam > /data/alignments/minimap2/NB03.stats

plot-bamstats \
    -p alignments/minimap2/plots_bamstats alignments/minimap2/NB03.stats

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

#### j1. Running Kraken2 using the flukraken2 database (Nanopore)

```
mkdir /data/classifications;

kraken2 --db /data/databases/flukraken2 /data/metagenome/flu_ont/BC01.fastq --report /data/classifications/BC01.kraken.report


```

#### j2. Prepping the database for Krona Plots


Try this out in pavian. First exit from the docker container with `exit`, then run 

```
docker container run -it --rm -p 8004:80  florianbw/pavian

```

:warning: If you receive an error about another container using the same port (like 8004) check it with `docker container ls`. You can also remove all containers that might be running with `docker container prune -f`

:warning: This next section is optional for the workshop. Feel free to try at home. It will require substantial filesize(s) to be downloaded so please be patient!

:warning:Requires internet to download the minikraken database. You can also get it from [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)


#### j3. Getting the minikraken database for Kraken2

```
mkdir -p /data/databases

wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz -O /data/databases/minikraken2.tar.gz; 

tar -xvzf --directory /data/databases/minikraken2 /data/databases/minikraken2.tar.gz

```


#### j4. Running Kraken2 using the flukraken2 database (Nanopore)

```
mkdir /data/classifications;

kraken2 --db /data/databases/minikraken2 --classified-out sample_metagenome#.fq metagenome/sample_metagenome.fastq --report classifications/sample_metagenome.kraken.report


```

#### j5. Prepping the database for Krona Plots

:warning:Requires internet connectivity

```
ktUpdateTaxonomy.sh /data/databases/minikraken2
```

#### j6. Creating a Krona Plot of your Tax calls

```
ktImportTaxonomy -i /data/classifications/ERR6913101.kraken.report -o /data/classifications/ERR6913101.kraken.html -tax /data/databases/minikraken2/

ktImportTaxonomy -i /data/classifications/sample_metagenome.kraken.report -o /data/classifications/sample_metagenome.kraken.html -tax /data/databases/minikraken2/

```

Now take a look at the resulting `classifications/sample_metagenome.kraken.html` file by double-clicking it to open in your browser automatically

#### k. Preparing to create a consensus 

*You need to use a different Docker image than sandbox for this. If you're in a Docker container already you need to type: `exit`*

#### k.1. Creating a consensus genome from Nanopore Reads using Medaka


:warning: This step is optional for your own data! All MinKNOW runs can be demultiplexed automatically after basecalling through the MinKNOW UI options


You may have already noticed that we have a demultiplexed set of files in the `/data/demux-fastq_pass` folder. However, this is not likely to happen when you have a traditional run from MinKNOW. Usually it will spit out a set of fastq files rather than just one single one per sample. To gather them all up, we can run one of the subcommands called `artic gather`. This will push all fastq files into a single one like what we see in `/data/demux-fastq_pass`

Let's first enter our new image by running. Remember that the name for the current directory is `$PWD` for Unix (Mac or Linux) and `$pwd` for Windows Powershell

To Demultiplex a run, use `guppy barcoder` on the `fastq_pass` folder of interest

### Windows Powershell

```
docker container run -w /data -v $pwd/test-data:/data  --rm --name articdemux genomicpariscentre/guppy guppy_barcoder --require_barcodes_both_ends -i /data/20200514_2000_X3_FAN44250_e97e74b4/fastq_pass -s /data/20200514_2000_X3_FAN44250_e97e74b4/demux --recursive
```

### Unix

```
docker container run -w /data -v $PWD/test-data:/data --rm --name articdemux genomicpariscentre/guppy guppy_barcoder --require_barcodes_both_ends -i /data/20200514_2000_X3_FAN44250_e97e74b4/fastq_pass -s /data/20200514_2000_X3_FAN44250_e97e74b4/demux --recursive
```

#### k.2 Gathering multiple fastqs 

*You need to use a different Docker image than sandbox for this. If you're in a Docker container already you need to type: `exit`*

Once you have a demultiplexed directory for your sample(s) you then need to gather them all into a single fastq file for the next step. 

**You have 2 options**

You only need to select one of the options (Option A or B) below

#### Option A


The easier method is to simply use the `artic` command from `staphb/artic-ncov2019` to gather up all fastq files in your directory (demultiplexed result) and gather all into one file. 


#### Windows Powershell

```
docker container run -w /data/20200514_2000_X3_FAN44250_e97e74b4/demux -v $pwd/test-data:/data --rm --name articgather staphb/artic-ncov2019 artic guppyplex --directory /data/20200514_2000_X3_FAN44250_e97e74b4/demux/barcode03 --output /data/demultiplexed/barcode03.fastq
```

#### Unix

```
docker container run -w /data/20200514_2000_X3_FAN44250_e97e74b4/demux -v $PWD/test-data:/data --rm --name articgather staphb/artic-ncov2019 artic guppyplex --directory /data/20200514_2000_X3_FAN44250_e97e74b4/demux/barcode03 --output /data/demultiplexed/barcode03.fastq
```

IF you look at the `demux-fastq_pass` folder, you should see that the `NB03.fastq` is the same length as the `demultiplexed/barcode03.fastq` file you made just now. 


```

mkdir -p test-data/demultiplexed

for fastq in $( ls test-data/20200514_2000_X3_FAN44250_e97e74b4/demux/barcode03 ); do \
    cat test-data/20200514_2000_X3_FAN44250_e97e74b4/demux/barcode03/$fastq ; \
done > test-data/demultiplexed/barcode03.fastq

```

You now have created a merged barcode from all the fastqs that are attributed to it! Congrats!

## l. Consensus Generation with medaka and artic

Let's start easy and run artic. This is not in your sandbox Docker Image so first, we need to `exit`

We first need to run 


#### Windows Powershell 

```

docker container run -w /data/consensus/artic -v $pwd/test-data:/data  --rm --name articconsensus staphb/artic-ncov2019 artic minion --medaka --medaka-model r941_min_high_g360  --strict --normalise 1000000 --scheme-directory /data/primer-schemes  --scheme-version V3 --read-file /data/demultiplexed/barcode03.fastq nCoV-2019/V3 barcode03;

```

#### Unix

```

docker container run -w /data/consensus/artic -v $PWD/test-data:/data  --rm --name articconsensus staphb/artic-ncov2019 artic minion --medaka --medaka-model r941_min_high_g360  --strict --normalise 1000000 --scheme-directory /data/primer-schemes  --scheme-version V3 --read-file /data/demultiplexed/barcode03.fastq nCoV-2019/V3 barcode03;

```

Finally, to make a report, we can run `multiqc .` in the /data/consensus folder to make a helpful report for variant information

#### Windows Powershell

```
docker container run -w /data/consensus -v $pwd/test-data:/data --rm --name articreport staphb/artic-ncov2019 bash -c "export LC_ALL=C.UTF-8; export LANG=C.UTF-8; multiqc . "

```

#### Unix

```

docker container run -w /data/consensus -v $PWD/test-data:/data  --rm --name articreport staphb/artic-ncov2019 bash -c "export LC_ALL=C.UTF-8; export LANG=C.UTF-8; multiqc . "

```

Your output files will be in `test-data/consensus/artic`, primarily you want the `<barcode_name>.consensus.fasta` and `multiqc_report.html` if you made it


[Medaka](https://github.com/nanoporetech/medaka) is a good tool for generating consensues and variant calling for most organisms. 

We have some demultiplexed SARS-CoV-2 test data available in the test directory. Lets try to make a consensus out of a sample's full fastq file. Let's ASSUME that we did not use any amplicons in the same fastq file we've been working on that is long reads

```

conda activate consensus && \
    mkdir -p /data/consensus/medaka && \
    cd /data/consensus/medaka

medaka_consensus \
    -i /data/demux-fastq_pass/NB03.fastq \
    -d /data/reference/nCoV-2019.reference.fasta \
    -o /data/consensus/medaka \
    -m r941_min_high_g360

```

 
## Creating your Own Ivar Consensus Runs using 3rd party Docker Images

But, what if we wanted to run our own images one by one NOT in an interactive environment. 


If we're working with a select few organisms utilizing a primer scheme for Ilumina, we should opt for the ivar bioinformatics pipeline. This pipeline can use both nanopolish and medaka as it's primary variant and consensus calling set of scripts. In addition, it will perform the necessary trimming/pre-processing for you on your demultiplexed fastq files. It also comes with several subcommands that can help prep your data

First, we need to get the primer schemes for SARS-CoV-2 (our organism of interest for this example). These primers will be used for the artic minion process. Luckily, we have this present in your test-data/primer-schemes


This is pulling in the set of primer schemes directly from the artic pipeline toolkit's set of primer schemes for EBOLA, SARS-CoV-2, and Nipah



1. 

| Command | Platform |
| ------- | -------- |
| `docker container run -w /data -v $pwd/test-data:/data --rm --name artic staphb/ivar bash -c "ivar trim -i /data/alignments/minimap2/NB03_sorted.bam -b /data/primer-schemes/nCoV-2019/V3/nCoV-2019.primer.bed -p /data/alignments/minimap2/NB03_trimmed_unsorted.bam"` | Windows Powershell |
| `docker container run -w /data -v $PWD/test-data:/data  --rm --name artic staphb/ivar bash -c "ivar trim -i /data/alignments/ERR6913101_alignments_sorted.bam -b /data/primer-schemes/nCoV-2019/V3/nCoV-2019.primer.bed -p /data/alignments/ERR6913101_trimmed_unsorted.bam"` | Unix |


2. 


| Command | Platform |
| ------- | -------- |
| `docker container run -w /data -v  $pwd/test-data:/data --rm --name artic staphb/ivar bash -c "samtools sort -o /data/alignments/minimap2/NB03_trimmed_unsorted.bam > /data/alignments/minimap2/NB03_trimmed_sorted.bam"` | Windows Powershell |
| `docker container run -w /data -v $PWD/test-data:/data --rm --name artic staphb/ivar bash -c "samtools sort -o /data/alignments/ERR6913101_trimmed_unsorted.bam > /data/alignments/ERR6913101_trimmed_sorted.bam"` | Unix |

3. 


| Command | Platform |
| ------- | -------- |
| `docker container run -w /data -v $pwd/test-data:/data --rm --name artic staphb/ivar bash -c "mkdir /data/variants/; samtools mpileup -A -aa -d 0 -Q 0 --reference /data/reference/nCoV-2019.reference.fasta  /data/alignments/ERR6913101_trimmed_sorted.bam > /data/variants/ERR6913101_pileup.txt"` | Windows Powershell |
| `docker container run -w /data -v $PWD/test-data:/data --rm --name artic staphb/ivar bash -c "mkdir /data/variants/; samtools mpileup -A -aa -d 0 -Q 0 --reference /data/reference/nCoV-2019.reference.fasta  /data/alignments/ERR6913101_trimmed_sorted.bam > /data/variants/ERR6913101_pileup.txt"` | Unix |


4. 

| Command | Platform |
| ------- | -------- |
| ```docker container run -w /data -v $pwd/test-data:/data --rm --name artic staphb/ivar bash -c "cat /data/variants/NB03_pileup.txt \|  ivar consensus  -p /data/consensus/consensus.fa  -m 10 -t 0.5 -n N" ``` | Windows Powershell |
| `docker container run -w /data -v $PWD/test-data:/data --rm --name artic staphb/ivar bash -c "cat /data/variants/NB03_pileup.txt \|  ivar consensus  -p /data/consensus/consensus.fa  -m 10 -t 0.5 -n N" ` | Unix |

## Creating your Own Medaka Consensus Runs using 3rd party Docker Images

So now we've used our single sandbox Docker Image. But what if, for example, you have a script that isn't available in this pre-made image? For example, `medaka`, which is as we described previously is available in sandbox but also as a standlone Docker Image from Docker hub [here](https://hub.docker.com/r/staphb/medaka)

Let's run a consensus run, like we did in the consensus pipeline without pre-installing it AND automatically running our command once it is downloaded. 


| Command | Platform |
| ------- | -------- |
| `docker container run -v $pwd/test-data:/data -w /data staphb/medaka bash -c "medaka_consensus -i /data/demux-fastq_pass/NB03.fastq -d /data/reference/nCoV-2019.reference.fasta -o /data/output/medaka_consensus -m r941_min_high_g303 -f"` | Windows Powershell |
| `docker container run -v $PWD/test-data:/data -w /data staphb/medaka bash -c "medaka_consensus -i /data/demux-fastq_pass/NB03.fastq -d /data/reference/nCoV-2019.reference.fasta -o /data/output/medaka_consensus -m r941_min_high_g303 -f"` | Unix |

Notice something interesting....

We can immediately run a command directly from our Shell, and have all dependencies automatically installed and our command run. This is one of the many strengths of using Docker as you have a fully runnable command with all required dependencies right out of the box!

That's it for the copy+paste section of the Docker portion of this workshop. Now that you've hopefully had a grasp on the inner-workings of Docker, we can move on to making and installation for your own pipeline(s)!

