#==================#
#   HEADER START   #
#==================#
### Created: Oct 29, 2018
### Author: Maciej Bak
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Description: Calculates the Differential Transcript Usage and Differential Gene Expression.
### Based on: https://f1000research.com/articles/7-952/v3
#==================#
#    HEADER END    #
#==================#


########################
### LOAD PACKAGES	 ###
########################
library("optparse")
library("tximport")
library("GenomicFeatures")
library("DRIMSeq")
library("DEXSeq")
library("DESeq2")
library("stageR")
library("edgeR")


#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option(c("--gtf"), action="store", type="character", default="", help="REQUIRED: GTF file with annotated transcripts", metavar="file"),
	make_option(c("--design_table"), action="store", type="character", default="", help="REQUIRED: BED file with regions of interest", metavar="file"),
	make_option(c("--output_dir"), action="store", type="character", default="", help="REQUIRED: Output directory", metavar="directory"),
	make_option(c("--alpha"), action="store", type="numeric", default=0.05, help="Alpha level for statistical testing."),
	make_option(c("--minimal_gene_expression"), action="store", type="numeric", default=10, help="Minimal gene expression for DRIMSeq filter."),
	make_option(c("--minimal_transcript_expression"), action="store", type="numeric", default=10, help="Minimal transcripts expression for DRIMSeq filter."),
	make_option(c("--minimal_proportion"), action="store", type="numeric", default=0.1, help="Minimal transcripts proportion for DRIMSeq filter."),
	make_option(c("--help"), action="store_true", default=FALSE, help="Show this information and die"),
	make_option(c("--verbose"), action="store_true", default=FALSE, help="Be Verbose")
)
## Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list, add_help_option=FALSE, description="")
opt <- parse_args(opt_parser)
stageR_alpha = opt$alpha
isoform_switch_alpha = opt$alpha
min_gene_expr = opt$minimal_gene_expression
min_feature_expr = opt$minimal_transcript_expression
min_feature_prop = opt$minimal_proportion

## Die if any required arguments are missing...
if ( opt$gtf=="" || opt$design_table=="" || opt$output_dir=="") {
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
if ( opt$verbose ){
	options(warn=0)
}

# CREATE THE OUTPUT DIRECTORY
outdir = strsplit(opt$output_dir,"/")[[1]][1]
plot_dir = paste(opt$output_dir,"DRIMSeq_plots",sep="/")
dir.create(opt$output_dir)
dir.create(plot_dir)


############################################
### Importing counts into R/Bioconductor ###
############################################

# We can use tximport to import the estimated counts, abundances and effective transcript lengths into R.
# The tximport package offers three methods for producing count matrices from transcript-level quantification files:
# (1) original estimated counts with effective transcript length information used as a statistical offset,
# (2) generating “counts from abundance” by scaling TPM abundance estimates per sample such that they sum to the total number of mapped reads (scaledTPM),
# (3) generating “counts from abundance” by scaling TPM abundance estimates first by the average effective transcript length over samples,
# and then per sample such that they sum to the total number of mapped reads (lengthScaledTPM).
# We will use scaledTPM for DTU detection, with the statistical motivation described below, and the original estimated counts with length offset for DGE detection.

 # read the design table...
samps <- read.csv(opt$design_table)
samps$condition <- factor(samps$condition)
# ...and the quantification files
files <- file.path(outdir, samps$sample_id, "quant.sf")
names(files) <- samps$sample_id

# We can then import transcript-level counts using tximport.
# For DTU analysis, we suggest generating counts from abundance, using the scaledTPM method.
# The countsFromAbundance option of tximport uses estimated abundances to generate roughly count-scaled data, such that each column will sum to the number of reads mapped for that library.
# By using scaledTPM counts, the estimated proportions fit by DRIMSeq, which are generated from counts, will be equivalent to proportions of the abundance of the isoforms.
# If instead of scaledTPM, we used the original estimated transcript counts (countsFromAbundance="no"), or if we used lengthScaledTPM transcript counts,
# then a change in transcript usage among transcripts of different length could result in a changed total count for the gene, even if there is no change in total gene expression.
# Briefly, this is because the original transcript counts and lengthScaledTPM transcript counts scale with transcript length, while scaledTPM transcript counts do not.
# A change in the total count for the gene, not corrected by an offset, could then bias the calculation of proportions and therefore confound DTU analysis.
# For testing DTU using DRIMSeq and DEXSeq, it is convenient if the count-scale data do not scale with transcript length within a gene.
# This could be corrected by an offset, but this is not easily implemented in the current DTU analysis packages.
# For gene-level analysis (DGE), we can use gene-level original counts with an average transcript length offset, gene-level scaledTPM counts,
# or gene-level lengthScaledTPM counts, as all of these three methods control for the length bias.
# A final note is that, the motivation for using scaledTPM counts hinges on the fact that estimated fragment counts scale with transcript length in fragmented RNA-seq data.
# If a different experiment is performed and a different quantification method used to produce counts per transcript which do not scale with transcript length,
# then the recommendation would be to use these counts per transcript directly.
# Examples of experiments producing counts per transcript that would potentially not scale with transcript length include counts of full-transcript-length
# or nearly-full-transcript-length reads, or counts of 3’ tagged RNA-seq reads aggregated to transcript groups.
# In either case, the statistical methods for DTU could be provided directly with the transcript counts.

txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
cts <- txi$counts
# keep only those transcripts that are expressed in at least one sample
cts <- cts[rowSums(cts) > 0,]

############################################
###     Transcript-to-gene mapping       ###
############################################

# Make a dataframe with mapping:
# GENE ID : TRANSCRIPT ID : NUMBER OF TRANSCRIPTS
txdb <- makeTxDbFromGFF(opt$gtf, format="auto")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

############################################
###       DTU analysis by DRIMSeq        ###
############################################

# Make sure that the database is complete
if (!all(rownames(cts) %in% txdf$TXNAME)){
	stop("Quantified transcripts not in the databse. Aborted.")
}
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
if (!all(rownames(cts) == txdf$TXNAME)){
	stop("Quantified transcripts not in the databse. Aborted.")
}

#In order to run DRIMSeq, we build a data.frame with the gene ID, the feature (transcript) ID, and then columns for each of the samples:
counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts, check.names = FALSE)
# We can now create a dmDSdata object from DRIMSeq package, with our counts and samps dataframe.
d <- dmDSdata(counts=counts, samples=samps)

# It will be useful to first filter the object, before running procedures to estimate model parameters.
# This greatly speeds up the fitting and removes transcripts that may be troublesome for parameter estimation,
# e.g. estimating the proportion of expression among the transcripts of a gene when the total count is very low.
# We first define n to be the total number of samples, and n.small to be the sample size of the smallest group.
# We use all three of the possible filters: for a transcript to be retained in the dataset, we require that
# (1) it has a count of at least min_feature_expr in at least n.small samples,
# (2) it has a relative abundance proportion of at least min_feature_prop in at least n.small samples,
# (3) the total count of the corresponding gene is at least min_gene_expr in all n samples.
# We used all three possible filters, whereas only the two count filters are used in the DRIMSeq vignette example code.

# It is important to consider what types of transcripts may be removed by the filters, and potentially adjust depending on the dataset.
# If n was large, it would make sense to allow perhaps a few samples to have very low counts,
# so lowering min_samps_gene_expr to some factor multiple (< 1) of n, and likewise for the first two filters for n.small.
# The second filter means that if a transcript does not make up more than x% of the gene’s expression for at least n.small samples, it will be removed.
# If this proportion seems too high, for example, if very lowly expressed isoforms are of particular interest, then the filter can be omitted or the min_feature_prop lowered.
# For a concrete example, if a transcript goes from a proportion of 0% in the control group to a proportion of 9% in the treatment group, this would be removed by the above 10% filter.

n <- nrow(samps)
n.small <- min(aggregate(sample_id ~ condition, samps, function(x) length(unique(x)))$sample_id)
d <- dmFilter(d,
                min_samps_feature_expr=n.small, min_feature_expr=min_feature_expr,
                min_samps_feature_prop=n.small, min_feature_prop=min_feature_prop,
                min_samps_gene_expr=n, min_gene_expr=min_gene_expr)

# The dmDSdata object only contains genes that have more that one isoform, which makes sense as we are testing for differential transcript usage.
# We can find out how many of the remaining genes have N isoforms by tabulating the number of times we see a gene ID, then tabulating the output again:
binned_by_isoform_number = table(table(counts(d)$gene_id))

# We create a design matrix, using a design formula and the sample information contained in the object, accessed via samples.
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))

# We then use the following three functions to estimate the model parameters and test for DTU.
# We first estimate the precision, which is related to the dispersion in the Dirichlet Multinomial model as dispersion = 1/(1+precision)
# Higher dispersion – counts more variable around their expected value – is associated with lower precision.
# For full details about the DRIMSeq model, one should read both the detailed software vignette and the publication.
# After estimating the precision, we fit regression coefficients and perform null hypothesis testing on the coefficient of interest.
# Because we have a simple two-group model, we test the coefficient associated with the difference between condition 2 and condition 1, called condition2.

set.seed(1)
d <- dmPrecision(d, design=design_full)
d <- dmFit(d, design=design_full)
d <- dmTest(d, coef="condition2")

#To build a results table, we run the results function.
#We can generate a single p-value per gene, which tests whether there is any differential transcript usage within the gene,
#or a single p-value per transcript, which tests whether the proportions for this transcript changed within the gene:
res <- DRIMSeq::results(d)
res.txp <- DRIMSeq::results(d, level="feature")

# Because the pvalue column may contain NA values, we use the following function to turn these into 1’s.
# The NA values would otherwise cause problems for the stage-wise analysis.
# From investigating these NA p-value cases for DRIMSeq, they all occur when one condition group has all zero counts for a transcript,
# but sufficient counts from the other condition group, and sufficient counts for the gene.
# DRIMSeq will not estimate a precision for such a gene.
# These all happen to be true positive genes for DTU in the simulation, where the isoform switch is total or nearly total.
# DEXSeq, shown in a later section, does not produce NA p-values for these genes.
# A potential fix would be to use a plug-in common or trended precision for such genes, but this is not implemented in the current version of DRIMSeq.
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

# Save the result tables
write.table(res, file = paste(opt$output_dir,"DRIMSeq_genes_DTU.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")
write.table(res.txp, file = paste(opt$output_dir,"DRIMSeq_transcripts_proportions.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")

# We can plot the estimated proportions for one of the significant genes, where we can see evidence of switching
signif_idx <- which(res$adj_pvalue < isoform_switch_alpha)

for (i in 1:length(signif_idx)){
	idx = signif_idx[i]
	pdf(paste(plot_dir,paste(res$gene_id[idx],"pdf",sep="."),sep="/"))
	print(plotProportions(d, res$gene_id[ which(res$adj_pvalue < isoform_switch_alpha)[i] ], "condition"))
	dev.off()
}

############################################
###    StageR analysis after DRIMSeq     ###
############################################

# A typical analysis of differential transcript usage would involve asking first: “which genes contain any evidence of DTU?”,
# and secondly, “which transcripts in the genes that contain some evidence may be participating in the DTU?”.
# Note that a gene may pass the first stage without exhibiting enough evidence to identify one or more transcripts that are participating in the DTU.
# The stageR package is designed to allow for such two-stage testing procedures, where the first stage is called a screening stage and the second stage a confirmation stage.
# We show below how stageR is used to detect DTU and how to interpret its output.
# We first construct a vector of p-values for the screening stage.
# Because of how the stageR package will combine transcript and gene names, we need to strip the gene and transcript version numbers from their Ensembl IDs
# (this is done by keeping only the first 15 characters of the gene and transcript IDs).

pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)

# We construct a one column matrix of the confirmation p-values:
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)

# We arrange a two column data.frame with the transcript and gene identifiers.
tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# The following functions then perform the stageR analysis.
# We must specify an alpha, which will be the overall false discovery rate target for the analysis, defined below.
# Unlike typical adjusted p-values or q-values, we cannot choose an arbitrary threshold later: after specifying alpha=0.05, we need to use 5% as the target in downstream steps.
# There are also convenience functions getSignificantGenes and getSignificantTx, which are demonstrated in the stageR vignette.

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                         pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=stageR_alpha)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                      onlySignificantGenes=FALSE)
})

# save the table
write.table(drim.padj, file = paste(opt$output_dir,"StageR_DRIMSeq.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")

# The final table with adjusted p-values summarizes the information from the two-stage analysis.
# Only genes that passed the filter are included in the table, so the table already represents screened genes.
# The transcripts with values in the column, transcript, less than stageR_alpha pass the confirmation stage on a target stageR_alpha overall false discovery rate, or OFDR.
# This means that, in expectation, no more than stageR_alpha of the genes that pass screening will either:
# (1) not contain any DTU, so be falsely screened genes, or
# (2) contain a falsely confirmed transcript.
# A falsely confirmed transcript is a transcript with an adjusted p-value less than stageR_alpha which does not exhibit differential usage across conditions.
# The stageR procedure allows us to look at both the genes that passed the screening stage and the transcripts with adjusted p-values less than our target alpha,
# and understand what kind of overall error rate this procedure entails.
# This cannot be said for an arbitrary procedure of looking at standard gene adjusted p-values and transcript adjusted p-values, where the adjustment was performed independently.

############################################
### Post-hoc filter on SD in proportions ###
############################################

# We found that DRIMSeq was sensitive to detect DTU, but could exceed its false discovery rate (FDR) bounds,
# particularly on the transcript-level tests, and that a post-hoc, non-specific filtering of the DRIMSeq transcript p-values and adjusted p-values improved the FDR and OFDR control.
# We considered the standard deviation (SD) of the per-sample proportions as a filtering statistic.
# This statistic does not use the information about which samples belong to which condition group.
# We set the p-values and adjusted p-values for transcripts with small per-sample proportion SD to 1.
# We do not recompute adjusted p-values, although we will provide the filtered p-values to the stageR procedure.
# We note that the p-values are no longer necessarily uniform after filtering out small effect size transcripts and genes,
# although we find that in our simulation at least, the filtering made the procedure more conservative:
# excluding transcripts with small SD of the per-sample proportions brought both the observed FDR and the observed OFDR closer to their nominal targets.

res.txp.filt <- DRIMSeq::results(d, level="feature")

smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}

filt <- smallProportionSD(d)
res.txp.filt$pvalue[filt] <- 1
res.txp.filt$adj_pvalue[filt] <- 1
write.table(res.txp.filt, file = paste(opt$output_dir,"DRIMSeq_transcripts_proportions_post_hoc_SD_adjusted.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")

# The above post-hoc filter is not part of the DRIMSeq modeling steps, and to avoid interfering with the modeling, we run it after DRIMSeq.
# The other three filters used before have been tested by the DRIMSeq package authors, and are therefore a recommended part of an analysis before the modeling begins.
# We do not apply this post-hoc filter to DEXSeq in this workflow, as DEXSeq seemed to be closer to controlling its FDR and OFDR in the evaluations,
# after using the DRIMSeq filters recommended in this workflow.

############################################
###       DTU analysis by DEXSeq         ###
############################################

# The DEXSeq package was originally designed for detecting differential exon usage, but can also be adapted to run on estimated transcript counts, in order to detect DTU.
# Using DEXSeq on transcript counts was evaluated previously, showing the benefits in FDR control from filtering lowly expressed transcripts for a transcript-level analysis.
# We benchmarked DEXSeq here, beginning with the DRIMSeq filtered object, as these filters are intuitive, they greatly speed up the analysis,
# and such filtering was shown to be beneficial in FDR control.
# The two factors of:
# (1) working on isoform counts rather than individual exons and
# (2) using the DRIMSeq filtering procedure
# dramatically increase the speed of DEXSeq, compared to running an exon-level analysis.
# Another advantage is that we benefit from the sophisticated bias models of Salmon,
# which account for drops in coverage on alternative exons that can otherwise throw off estimates of transcript abundance.
# A disadvantage over the exon-level analysis is that we must know in advance all of the possible isoforms that can be generated from a gene locus,
# all of which are assumed to be contained in the annotation files (FASTA and GTF).

# We first load the DEXSeq package and then build a DEXSeqDataSet from the data contained in the dmDStest object
# (the class of the DRIMSeq object changes as the results are added).
# The design formula of the DEXSeq-DataSet here uses the language “exon” but this should be read as “transcript” for our analysis.
# DEXSeq will test – after accounting for total gene expression for this sample and for the proportion of this transcript relative to the others –
# whether there is a condition-specific difference in the transcript proportion relative to the others.

sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                       sampleData=sample.data,
                       design=~sample + exon + condition:exon,
                       featureID=counts(d)$feature_id,
                       groupID=counts(d)$gene_id)

# The following functions run the DEXSeq analysis.
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)

# We then extract the results table, not filtering on mean counts (as we have already conducted filtering via DRIMSeq functions).
# We compute a per-gene adjusted p-value, using the perGeneQValue function, which aggregates evidence from multiple tests within a gene to a single p-value for the gene
# and then corrects for multiple testing across genes.
# Other methods for aggregative evidence from the multiple tests within genes have been discussed in a recent publication and may be substituted at this step.
# Finally, we build a simple results table with the per-gene adjusted p-values.

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])

write.table(dxr, file = paste(opt$output_dir,"DEXSeq_DTU.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")

############################################
###    StageR analysis after DEXSeq      ###
############################################

# we run code very similar to that used above for DRIMSeq two-stage testing, with a target alpha=stageR_alpha.
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# The following three functions provide a table with the OFDR control described above.
# To repeat, the set of genes passing screening should not have more than stageR_alpha of either:
# - genes which have in fact no DTU or
# - genes which contain a transcript with an adjusted p-value less than stageR_alpha which do not participate in DTU.

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=stageR_alpha)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
})

write.table(dex.padj, file = paste(opt$output_dir,"StageR_DEXSeq.tsv",sep="/"), row.names=FALSE, na="",col.names=TRUE, sep="\t")

############################################
###        DGE analysis by DESeq         ###
############################################

# In the final section we demonstrate how differential transcript usage, summarized to the gene level, can be visualized with respect to differential gene expression analysis results.
# We use tximport and summarize counts to the gene level and compute an average transcript length offset for count-based methods.
# We will then show code for using DESeq2 and edgeR to assess differential gene expression.

# The following line of code generate an object txi.g which contains the gene-level counts, abundances and average transcript lengths.
txi.g <- tximport(files, type="salmon", tx2gene=txdf[,2:1])

# We load the DESeq2 package and build a DESeqDataSet from txi.g, providing also the sample information and a design formula.
dds <- DESeqDataSetFromTximport(txi.g, samps, ~condition)

#The following two lines of code run the DESeq2 analysis.
dds <- DESeq(dds)
dres <- DESeq2::results(dds)

write.table(dres, file = paste(opt$output_dir,"DESeq_DGE.tsv",sep="/"), row.names=TRUE, na="",col.names=TRUE, sep="\t")

############################################
###        DGE analysis by edgeR         ###
############################################

# We can also perform differential gene expression analysis using edgeR as the inference engine.
# The following code incorporates the average transcript length matrix as an offset for an edgeR analysis.
cts.g <- txi.g$counts
normMat <- txi.g$length
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts.g/normMat)) + log(colSums(cts.g/normMat))
y <- DGEList(cts.g)
y <- scaleOffset(y, t(t(log(normMat)) + o))
keep <- filterByExpr(y)
y <- y[keep,]

# The basic edgeR model fitting and results extraction can be accomplished with the following lines:
y <- estimateDisp(y, design_full)
fit <- glmFit(y, design_full)
lrt <- glmLRT(fit)
tt <- topTags(lrt, n=nrow(y), sort="none")[[1]]
write.table(tt, file = paste(opt$output_dir,"edgeR_DGE.tsv",sep="/"), row.names=TRUE, na="",col.names=TRUE, sep="\t")
