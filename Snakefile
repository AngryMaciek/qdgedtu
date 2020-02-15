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

# local rules
localrules: all, create_output_dir, create_index_dir, extract_decoys, \
    concatenate_transcriptome_and_genome, merge_TPM_tables, \
    prepare_table_for_R_workflow

# get all sample names from the design table
def get_samples():
    design_table = pd.read_csv(config["design_table"], sep="\t")
    return list(design_table["sample"])

# get mate1 of the RNA-Seq data for a given sample
def get_mate_1(wildcards):
    design_table = pd.read_csv(config["design_table"], sep="\t", index_col=0)
    return str(design_table.loc[wildcards.sample]["fq1"])

# get mate2 of the RNA-Seq data for a given sample
def get_mate_2(wildcards):
    design_table = pd.read_csv(config["design_table"], sep="\t", index_col=0)
    return str(design_table.loc[wildcards.sample]["fq2"])

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TSV_all_genes = \
            expand(os.path.join("{output_dir}", "genes.tsv"), \
                output_dir=config["output_dir"]),
        TSV_all_transcripts = \
            expand(os.path.join("{output_dir}", "transcripts.tsv"), \
                output_dir=config["output_dir"]),
        TSV_edgeR = \
            expand(os.path.join("{output_dir}", "edgeR_DGE.tsv"), \
                output_dir=config["output_dir"]),
        TSV_DESeq = \
            expand(os.path.join("{output_dir}", "DESeq_DGE.tsv"), \
                output_dir=config["output_dir"]),
        TSV_DRIMSeq = \
            expand(os.path.join("{output_dir}", "StageR_DRIMSeq.tsv"), \
                output_dir=config["output_dir"]),
        TSV_DEXSeq = \
            expand(os.path.join("{output_dir}", "StageR_DEXSeq.tsv"), \
                output_dir=config["output_dir"])

##############################################################################
### Create directories for the results
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "output_dir_created"))
    params:
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log")
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log")
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Create directories for the transcriptome index
##############################################################################

rule create_index_dir:
    output:
        TMP_index = temp(os.path.join("{index_dir}", "index_dir_created"))
    params:
        DIR_index_dir = "{index_dir}",
        DIR_cluster_log = os.path.join("{index_dir}", "cluster_log")
    log:
        DIR_local_log = os.path.join("{index_dir}", "local_log")
    shell:
        """
        mkdir -p {params.DIR_index_dir}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_index}
        """

##############################################################################
### Extract transcriptome & biuld index for it
##############################################################################

rule extract_transcriptome:
    input:
        TMP_index = os.path.join("{index_dir}", "index_dir_created")
    output:
        FASTA_transcriptome = \
            os.path.join("{index_dir}", "transcriptome.fasta")
    params:
        GTF_gtf = config["gtf"],
        FASTA_genome = config["genome"],
        LOG_cluster_log = \
            os.path.join("{index_dir}", "cluster_log", \
                "extract_transcriptome.log"),
        queue = "30min",
        time = "00:10:00"
    resources:
        threads = 1,
        mem = 5000
    log:
        LOG_local_log = \
            os.path.join("{index_dir}", "local_log", \
                "extract_transcriptome.log")
    benchmark:
        os.path.join("{index_dir}", "local_log", \
            "extract_transcriptome_benchmark.log")
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        gffread {params.GTF_gtf} \
        -g {params.FASTA_genome} \
        -w {output.FASTA_transcriptome} \
        &> {log.LOG_local_log};
        """

rule extract_decoys:
    input:
        TMP_index = os.path.join("{index_dir}", "index_dir_created")
    output:
        TXT_decoys = \
            os.path.join("{index_dir}", "decoys.txt")
    params:
        FASTA_genome = config["genome"]
    log:
        LOG_local_log = \
            os.path.join("{index_dir}", "local_log", \
                "extract_decoys.log")
    benchmark:
        os.path.join("{index_dir}", "local_log", \
            "extract_decoys_benchmark.log")
    shell:
        """
        (grep "^>" <{params.FASTA_genome} \
        | cut -d " " -f 1 > {output.TXT_decoys} && \
        sed -i.bak -e 's/>//g' {output.TXT_decoys}) \
        2> {log.LOG_local_log};        
        """

rule concatenate_transcriptome_and_genome:
    input:
        FASTA_transcriptome = \
            os.path.join("{index_dir}", "transcriptome.fasta")
    output:
        FASTA_merged = \
            os.path.join("{index_dir}", "genome_transcriptome.fasta")
    params:
        FASTA_genome = config["genome"]
    log:
        LOG_local_log = \
            os.path.join("{index_dir}", "local_log", \
                "concatenate_transcriptome_and_genome.log")
    benchmark:
        os.path.join("{index_dir}", "local_log", \
            "concatenate_transcriptome_and_genome_benchmark.log")
    shell:
        """
        cat {input.FASTA_transcriptome} {params.FASTA_genome} \
        1> {output.FASTA_merged} \
        2> {log.LOG_local_log};   
        """

rule index_transcriptome:
    input:
        TXT_decoys = \
            os.path.join("{index_dir}", "decoys.txt"),
        FASTA_merged = \
            os.path.join("{index_dir}", "genome_transcriptome.fasta")
    output:
        DIR_transcriptome = \
            directory(os.path.join("{index_dir}", "index"))
    params:
        LOG_cluster_log = \
            os.path.join("{index_dir}", "cluster_log", \
                "index_transcriptome.log"),
        queue = "6hours",
        time = "06:00:00",
    resources:
        threads = 4,
        mem = 50000
    log:
        LOG_local_log = \
            os.path.join("{index_dir}", "local_log", \
                "index_transcriptome.log")
    benchmark:
        os.path.join("{index_dir}", "local_log", \
            "index_transcriptome_benchmark.log")
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        salmon index \
        -t {input.FASTA_merged} \
        -d {input.TXT_decoys} \
        --threads {resources.threads} \
        -i {output.DIR_transcriptome} \
        &> {log.LOG_local_log};
        """

##############################################################################
### Quantify transcripts expression
##############################################################################

rule salmon_quantify:
    input:
        TMP_output = os.path.join("{output_dir}", "output_dir_created"),
        DIR_transcriptome = \
            expand(os.path.join("{index_dir}", "index"), index_dir=config["transcriptome_index"])

        #DIR_transcriptome = \
        #    os.path.join("{output_dir}", "transcriptome_index")
    output:
        TSV_salmon_out = \
            os.path.join("{output_dir}", "{sample}", "quant.sf")
    params:
        GTF_gtf = config["gtf"],
        STRING_mate_1 = lambda wildcards: get_mate_1(wildcards),
        STRING_mate_2 = lambda wildcards: get_mate_2(wildcards),
        libType = "A",
        seqtype = config["seqtype"],
        DIR_salmon_dir = os.path.join("{output_dir}", "{sample}"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "salmon_quantify_{sample}.log"),
        queue = "6hours",
        time = "02:00:00",
    resources:
        threads = 4,
        mem = 50000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "salmon_quantify_{sample}.log")
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "salmon_quantify_{sample}_benchmark.log")
    conda:
        "env_yaml/quantification.yaml"
    shell:
        """
        if [ {params.seqtype} == paired_end ]
        then
            salmon quant \
            --index {input.DIR_transcriptome} \
            --geneMap {params.GTF_gtf} \
            --libType {params.libType} \
            -1 {params.STRING_mate_1} \
            -2 {params.STRING_mate_2} \
            --seqBias \
            --gcBias \
            --posBias \
            --validateMappings \
            --threads {threads} \
            --output {params.DIR_salmon_dir} \
            &> {log.LOG_local_log};
        else
            salmon quant \
            --index {input.DIR_transcriptome} \
            --geneMap {params.GTF_gtf} \
            --libType {params.libType} \
            -r {params.STRING_mate_1} \
            --seqBias \
            --posBias \
            --validateMappings \
            --threads {threads} \
            --output {params.DIR_salmon_dir} \
            &> {log.LOG_local_log};
        fi
       """

##############################################################################
### Merged expressed data
##############################################################################

rule merge_TPM_tables:
    input:
        TSV_salmon_out = \
            expand(os.path.join("{output_dir}", "{sample}", "quant.sf"), \
                output_dir=config["output_dir"], sample=get_samples())
    output:
        TSV_all_genes = os.path.join("{output_dir}", "genes.tsv"),
        TSV_all_transcripts = os.path.join("{output_dir}", "transcripts.tsv"),
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "merge_TPM_tables.log")
    run:
        # merge TPM values for all genes/transcripts for all samples
        quant_genes_list = []
        quant_transcripts_list = []
        design_table = pd.read_csv(config["design_table"], sep="\t")
        for i,row in design_table.iterrows():
            path = \
                os.path.join(config["output_dir"], row["sample"], "quant.sf")
            df = pd.read_csv(path, sep="\t", index_col=0)[["TPM"]]
            df.columns = [row["sample"]]
            quant_transcripts_list.append(df)
            path = \
                os.path.join(\
                    config["output_dir"], row["sample"], "quant.genes.sf")
            df = pd.read_csv(path, sep="\t", index_col=0)[["TPM"]]
            df.columns = [row["sample"]]
            quant_genes_list.append(df)
        t_df = pd.concat(quant_transcripts_list, axis=1)
        t_df.index.name = "ID"
        t_df.to_csv(output.TSV_all_transcripts, sep="\t")
        g_df = pd.concat(quant_genes_list, axis=1)
        g_df.index.name = "ID"
        g_df.to_csv(output.TSV_all_genes, sep="\t")

##############################################################################
### DTU + DGE workflow
##############################################################################

rule prepare_table_for_R_workflow:
    input:
        TSV_salmon_out = \
            expand(os.path.join("{output_dir}", "{sample}", "quant.sf"), \
                output_dir=config["output_dir"], sample=get_samples())
    output:
        CSV_samples_table = \
            os.path.join("{output_dir}", "workflow_table.csv")
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "prepare_table_for_R_workflow.log")
    run:
        # parse design table and create a simplified version for the
        # workflow below
        design_table = pd.read_csv(config["design_table"], sep="\t")
        with open(output.CSV_samples_table, "w") as csv_table:
            csv_table.write("sample_id,condition\n")
            for i,row in design_table.iterrows():
                condition = 1 if row["condition"]=="untreated" else 2
                csv_table.write(row["sample"]+","+str(condition)+"\n")

rule DTU_and_DGE_workflow:
    input:
        CSV_samples_table = \
            os.path.join("{output_dir}", "workflow_table.csv")
    output:
        TSV_edgeR = os.path.join("{output_dir}", "edgeR_DGE.tsv"),
        TSV_DESeq = os.path.join("{output_dir}", "DESeq_DGE.tsv"),
        TSV_DRIMSeq = os.path.join("{output_dir}", "StageR_DRIMSeq.tsv"),
        TSV_DEXSeq = os.path.join("{output_dir}", "StageR_DEXSeq.tsv")
    params:
        SCRIPT = \
            os.path.join(\
                config["scripts_directory"], "mb_DGE_DTU_workflow.r"),
        DIR_out_dir = "{output_dir}",
        minimal_gene_expression = \
            config["minimal_gene_expression"],
        minimal_transcript_expression = \
            config["minimal_transcript_expression"],
        minimal_transcripts_proportion = \
            config["minimal_transcripts_proportion"],
        alpha = config["statistical_alpha"],
        GTF_gtf = config["gtf"],
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "DTU_and_DGE_workflow.log"),
        queue = "6hours",
        time = "06:00:00"
    resources:
        threads = 1,
        mem = 10000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "DTU_and_DGE_workflow.log")
    benchmark:
        os.path.join("{output_dir}", "local_log", \
            "DTU_and_DGE_workflow_benchmark.log")
    conda:
        "env_yaml/dgedtu.yaml"
    shell:
        """
        Rscript {params.SCRIPT} \
        --gtf {params.GTF_gtf} \
        --design_table {input.CSV_samples_table} \
        --output_dir {params.DIR_out_dir} \
        --alpha {params.alpha} \
        --minimal_gene_expression \
        {params.minimal_gene_expression} \
        --minimal_transcript_expression \
        {params.minimal_transcript_expression} \
        --minimal_proportion \
        {params.minimal_transcripts_proportion} \
        &> {log.LOG_local_log};
        """
