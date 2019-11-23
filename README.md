# qdgedtu
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[General information about the project]

## Snakemake pipeline execution
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires Python 3 and can be most easily installed via the bioconda package from the anaconda cloud service.

### Step 1: Download and installation of Miniconda3
Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  source .bashrc
  ```

macOS:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  source .bashrc
  ```

### Step 2: Snakemake installation

Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest version:
  ```bash
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml 

Write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are four scripts to start the pipeline, depending on whether you want to run locally/on a SLURM computational cluster and if you prefer to use conda virtual environment or a specific container from dockerhub. In order to use singularity containers you need to have `singularity` installed. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  bash snakemake_local_run_singularity_container.sh
  bash snakemake_cluster_run_singularity_container.sh
  ```

## License

GPL v3.0

--------------------


# DGE_DTU

Transcript quantification with salmon; Differential Gene Expression and Differential Transcripts Usage according to: https://f1000research.com/articles/7-952/v3

# Dependencies
- bash
- git
- wget
- conda (instructions below)

# Download and install miniconda 3 (if not already installed)

This pipeline has been tested with conda version 4.6.7.

## Download installation file

for Linux:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

for Mac OSX:
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

## Installation

Make sure that you run the 'bash' shell and execute:

for Linux:
```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

for Mac:
```bash
bash Miniconda3-latest-MacOSX-x86_64.sh
```

# Clone the repository and install snakemake in a virtual environment

```bash
git clone https://git.scicore.unibas.ch/AnnotationPipelines/dge_dtu.git
cd dge_dtu
conda create -n dge_dtu_snakemake -c conda-forge -c bioconda snakemake=5.4
conda activate dge_dtu_snakemake
```

# Download resources

Please download the annotation files as following:

Download the genome from ENSEMBL

```bash
mkdir -p RESOURCES
wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
zcat Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > RESOURCES/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
rm Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

Download annotation file
```bash
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.chr.gtf.gz
zcat Homo_sapiens.GRCh38.90.chr.gtf.gz > RESOURCES/Homo_sapiens.GRCh38.90.chr.gtf
rm Homo_sapiens.GRCh38.90.chr.gtf.gz
```

# User's input

The input consist of two files and only these two files should be adjusted per any workflow execution:

* Main config file for snakemake: `config.yaml`
which contains all the paths needed for successful workflow execution as well as all the parameters for
gene/transcript expression analysis. Please go carefully through all the entries and adjust all of them according to your desires.

* A design table with detailed information about sequencing data. A sample design table is available in the repository as `sample_design_table.tsv`.
The design table is  a TSV file with 4 columns; column names are fixed; first column stands as sample ID; last column represents the condition ('treated' or 'untreated' only);
2nd and 3rd column are paths to the fastq/fasta files and their order does not matter.
In case the sequencing data are single-end please provide paths only in the 2nd columnd and leave the 3rd column with empty strings.


# Contributions

Contributions are more than welcome, but please follow the following instructions:

Create a new branch
```bash
git checkout -b your_feature_branch
```

Add your changes and push them remotely

```bash
git push origin your_feature_branch
```

Create a merge request via the gitlab website: https://git.scicore.unibas.ch/AnnotationPipelines/dge_dtu/merge_requests
