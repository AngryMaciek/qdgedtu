##############################################################################
#
#   Snakemake pipeline:
#   Gene/transcripts quantification, Differential Gene Expression
#   and Differential Transcript Usage.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 22-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################

# imports
import sys
import os
import pandas as pd




'''

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_random_samples = os.path.join("{output_dir}", "random_samples"),
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_random_samples}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Sample some random data
##############################################################################

rule generate_files:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        SCRIPT = \
            os.path.join(config["src_dir"], "mb_random_sample.py")
    output:
        TXT_random_sample = \
            os.path.join("{output_dir}", "random_samples", "{file}")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "generate_files_{file}.log"),
        queue = "30min",
        time = "0:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "generate_files_{file}.log"),
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "generate_files_{file}_benchmark.log")
    conda:
        "packages.yaml"
    singularity:
        ""
    shell:
        """
        python {input.SCRIPT} \
        --outfile {output.TXT_random_sample} \
        &> {log.LOG_local_log};
        """

##############################################################################
### Merge the results
##############################################################################

rule merge_results:
    input:
        TXT_result_files = \
            lambda wildcards: [os.path.join(wildcards.output_dir,
                "random_samples", f) for f in config["samples_filenames"]]
    output:
        TXT_final_results = os.path.join("{output_dir}", "results.txt")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/merge_results.log"),
        queue = "30min",
        time = "00:05:00"
    resources:
        threads = 1,
        mem = 5000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_results.log")
    run:
        # read all the sampled numbers:
        numbers = []
        for i in input.TXT_result_files:
            with open(i) as f:
                numbers.append(f.read().splitlines()[0])
        # save into one file:
        with open(output.TXT_final_results, "w") as outfile:
                outfile.write("\n".join(numbers))

'''









# local rules
localrules: all, create_output_dir, merge_TPM_tables, prepare_table_for_R_workflow

#localrules: copy_configfile, merge_TPM_tables, prepare_table_for_R_workflow, final

def get_samples():
    design_table = pd.read_csv(config["design_table"], sep="\t")
    return list(design_table["sample"])

def get_mate_1(wildcards):
    design_table = pd.read_csv(config["design_table"], sep="\t", index_col=0)
    return str(design_table.loc[wildcards.sample]["fq1"])

def get_mate_2(wildcards):
    design_table = pd.read_csv(config["design_table"], sep="\t", index_col=0)
    return str(design_table.loc[wildcards.sample]["fq2"])

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        all_genes = config["output_dir"] + "/genes.tsv",
        all_transcripts = config["output_dir"] + "/transcripts.tsv",
        edgeR = config["output_dir"] + "/edgeR_DGE.tsv",
        DESeq = config["output_dir"] + "/DESeq_DGE.tsv",
        DRIMSeq = config["output_dir"] + "/StageR_DRIMSeq.tsv",
        DEXSeq = config["output_dir"] + "/StageR_DEXSeq.tsv",
        #TXT_final_results = \
        #    expand(os.path.join("{output_dir}", "results.txt"),
        #        output_dir=config["output_dir"])

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """











#################################################################################
### Extract transcriptome & biuld index for it
#################################################################################

rule extract_transcriptome:
    input:
        gtf = config["gtf"],
        genome = config["genome"],
        TMP_output = os.path.join(config["output_dir"], "dir_created"),
    output:
        transcriptome = "{output_dir}" + "/transcriptome.fasta"
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/extract_transcriptome.log"),
        queue = "30min",
        time = "00:10:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/extract_transcriptome.log"),
    resources:
        threads = 8,
        mem = 5000
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        gffread {input.gtf} \
        -g {input.genome} \
        -w {output.transcriptome} \
        &> {log.LOG_local_log};
        """











rule index_transcriptome:
    input:
        transcriptome = config["output_dir"] + "/transcriptome.fasta"
    output:
        transcriptome_index = directory("{output_dir}" + "/transcriptome_index")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/index_transcriptome.log"),
        queue = "30min",
        time = "00:15:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/index_transcriptome.log"),
    resources:
        threads = 4,
        mem = 10000
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        salmon index \
        -t {input.transcriptome} \
        -i {output.transcriptome_index} \
        --type quasi \
        -k 31 \
        &> {log.LOG_local_log};
        """

#################################################################################
### Quantify transcripts expression
#################################################################################

rule salmon_quantify:
    input:
        index = config["output_dir"] + "/transcriptome_index",
        gtf = config["gtf"]
    output:
        salmon_out = "{output_dir}" + "/{sample}/quant.sf",
    params:
        mate_1 = lambda wildcards: get_mate_1(wildcards),
        mate_2 = lambda wildcards: get_mate_2(wildcards),
        libType = "A",
        seqtype = config["seqtype"],
        salmon_dir = config["output_dir"] + "/{sample}",
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/salmon_quantify_{sample}.log"),
        queue = "6hours",
        time = "02:00:00",
    resources:
        threads = 8,
        mem = 10000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/salmon_quantify_{sample}.log"),
    threads:    6
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        if [ {params.seqtype} == paired_end ]
        then
            salmon quant \
            --index {input.index} \
            --geneMap {input.gtf} \
            --libType {params.libType} \
            -1 {params.mate_1} \
            -2 {params.mate_2} \
            --seqBias \
            --threads {threads} \
            --output {params.salmon_dir} \
            &> {log.LOG_local_log};
        else
            salmon quant \
            --index {input.index} \
            --geneMap {input.gtf} \
            --libType {params.libType} \
            -r {params.mate_1} \
            --seqBias \
            --threads {threads} \
            --output {params.salmon_dir} \
            &> {log.LOG_local_log};
        fi        
       """

rule merge_TPM_tables:
    input:
        salmon_out = expand(config["output_dir"] + "/{sample}/quant.sf", output_dir=config["output_dir"], sample=get_samples())
    output:
        all_genes = "{output_dir}" + "/genes.tsv",
        all_transcripts = "{output_dir}" + "/transcripts.tsv",
    params:
        queue = "30min",
        time = "00:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/merge_TPM_tables.log"),
    resources:
        threads = 1,
        mem = 5000
    run:
        quant_genes_list = []
        quant_transcripts_list = []
        design_table = pd.read_csv(config["design_table"],sep="\t")
        for i,row in design_table.iterrows():
            path = os.path.join(config["output_dir"],row["sample"],"quant.sf")
            df = pd.read_csv(path,sep="\t",index_col=0)[["TPM"]]
            df.columns = [row["sample"]]
            quant_transcripts_list.append(df)
            path = os.path.join(config["output_dir"],row["sample"],"quant.genes.sf")
            df = pd.read_csv(path,sep="\t",index_col=0)[["TPM"]]
            df.columns = [row["sample"]]
            quant_genes_list.append(df)
        x = pd.concat(quant_transcripts_list,axis=1)
        x.index.name = "ID"
        x.to_csv(output.all_transcripts,sep="\t")
        x = pd.concat(quant_genes_list,axis=1)
        x.index.name = "ID"
        x.to_csv(output.all_genes,sep="\t")

#################################################################################
### DTU + DGE workflow
#################################################################################

rule prepare_table_for_R_workflow:
    input:
        salmon_out = expand(config["output_dir"] + "/{sample}/quant.sf", output_dir=config["output_dir"], sample=get_samples())
    output:
        table = "{output_dir}"+ "/workflow_table.csv"
    params:
        queue = "30min",
        time = "00:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/prepare_table.log"),
    resources:
        threads = 1,
        mem = 5000
    run:
        design_table = pd.read_csv(config["design_table"],sep="\t")
        with open(output.table,"w") as csv_table:
            csv_table.write("sample_id,condition\n")
            for i,row in design_table.iterrows():
                condition = 1 if row["condition"]=="untreated" else 2
                csv_table.write(row["sample"]+","+str(condition)+"\n")

rule DTU_and_DGE_workflow:
    input:
        table = config["output_dir"] + "/workflow_table.csv",
        SCRIPT = os.path.join(config["scripts_directory"], "mb_DTU_and_DGE_workflow.R")
    output:
        edgeR = "{output_dir}" + "/edgeR_DGE.tsv",
        DESeq = "{output_dir}" + "/DESeq_DGE.tsv",
        DRIMSeq = "{output_dir}" + "/StageR_DRIMSeq.tsv",
        DEXSeq = "{output_dir}" + "/StageR_DEXSeq.tsv",
    params:
        out_dir = config["output_dir"],
        minimal_gene_expression = config["minimal_gene_expression"],
        minimal_transcript_expression = config["minimal_transcript_expression"],
        minimal_transcripts_proportion = config["minimal_transcripts_proportion"],
        alpha = config["statistical_alpha"],
        gtf = config["gtf"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/DTU_and_DGE_workflow.log"),
        queue = "6hours",
        time = "06:00:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log/DTU_and_DGE_workflow.log"),
    resources:
        threads = 1,
        mem = 10000
    conda:
        "env_yaml/dgedtu.yaml"
    shell:
        """
        Rscript {input.SCRIPT} \
        --gtf {params.gtf} \
        --design_table {input.table} \
        --output_dir {params.out_dir} \
        --alpha {params.alpha} \
        --minimal_gene_expression {params.minimal_gene_expression} \
        --minimal_transcript_expression {params.minimal_transcript_expression} \
        --minimal_proportion {params.minimal_transcripts_proportion} \
        &> {log.LOG_local_log};
        """
