# Snakemake pipeline for transcripts expression analyses
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

This is a small snakemake pipeline I have put together for quantification of transcripts expression, differential gene expression and differential transcripts usage from RNA-Seq data.  
Transcript quantification is performed with [salmon](https://combine-lab.github.io/salmon/).  
Differential Gene Expression ([DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)) and Differential Transcripts Usage ([DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) and [DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html)) according to: https://f1000research.com/articles/7-952/v3

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
### Step 2: Pandas and Snakemake installation

To execute the workflow one would require pandas python library and snakemake workflow menager.  
Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest versions:
  ```bash
  conda install -c conda-forge pandas
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml. Apart from the genome in FASTA and genomic annotation in GTF the user needs to provide a TSV design table with detailed information about sequencing data. A sample design table is available in the repository as `sample_design_table.tsv`.
Note that:
* the design table requires 4 columns; column names are fixed;
* first column stands as sample ID; last column represents the condition ('treated' or 'untreated' only);
* 2nd and 3rd column are paths to the fastq/fasta files and their order does not matter.
* in case the sequencing data are single-end please provide paths only in the 2nd columnd and leave the 3rd column with empty strings.

Once the metadata are ready write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are two scripts to start the pipeline, depending on whether you want to run locally or on a SLURM computational cluster. In order to execute the workflow snakemake automatically creates internal conda virtual environments and installs software from anaconda cloud service. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  ```

## License

GPL v3.0
