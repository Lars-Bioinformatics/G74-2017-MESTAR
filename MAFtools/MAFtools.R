if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("maftools")
install.packages("mclust")
install.packages("NMF")
install.packages("pheatmap")
install.packages("barplot3d")


library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
# library(tidyverse)

setwd("/work/Data/Connor/mutect2_vcf/")

maf.files = list.files("mutect2_maf_PASSonly", full.names = T)
laml = merge_mafs(maf.files)
# laml = read.maf(maf = "G74-NP1_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY_vs_G74-NN1_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY_somatic_mutect2_filterFlag.maf")

#Shows sample summry - count of mutations
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Shows all fields in MAF
getFields(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')


plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = laml, top = 50, bgCol = "lightgrey")
oncoplot(maf = laml, top = 50)

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

lollipopPlot(maf = laml, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'TP53', showDomainLabel = FALSE, labelPos = 882)

for (maf.file in maf.files[9:10]){ # Only files with changes: maf.files[9:10]
  maf.single = read.maf(maf.file)
  rainfallPlot(maf = maf.single, detectChangePoints = TRUE, pointSize = 0.4)
}

laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

laml@data$t_vaf = laml@data$t_alt_count/laml@data$t_depth
plotVaf(maf = laml, vafCol = 't_vaf', top = 20)
plotVaf(maf = laml, vafCol = 't_vaf')



# 9 Analysis

# 9.1 Somatic interactions
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))

# 9.2 Detecting cancer driver genes based on positional clustering
laml.sig = oncodrive(maf = laml, AACol = 'HGVSp', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# 9.3 Adding and summarizing pfam domains
laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp_Short', top = 10)
laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp', top = 10)

#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]

#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]

# 9.7 Drug-Gene Interactions
dgi = drugInteractions(maf = laml, fontSize = 0.75)

# 9.8 Oncogenic Signaling Pathways
OncogenicPathways(maf = laml)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")
PlotOncogenicPathways(maf = laml, pathways = "TP53")


# 9.9 Tumor heterogeneity and MATH scores
# 9.9.1 Heterogeneity in tumor samples
tum_sample_barcodes = unique(laml@data$Tumor_Sample_Barcode)
for (tsb in tum_sample_barcodes){
  inferred.het = inferHeterogeneity(maf = laml, tsb = tsb, vafCol = 't_vaf')
  plotClusters(clusters = inferred.het)
}

d = density(laml@data$t_vaf)
stats::plot(d)


# 9.10 Mutational Signatures
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

laml.tnm = trinucleotideMatrix(maf = laml, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# 9.10.1 APOBEC Enrichment estimation
# 9.10.2 Differences between APOBEC enriched and non-enriched samples
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)


# 9.10.3 Signature analysis
# Note: If either extractSignatures or estimateSignatures stops in between, 
#       its possibly due to low mutation counts in the matrix. In that case rerun
#       the functions with pConstant argument set to small positive value (e.g, 0.1).
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
dev.off()
plotCophenetic(res = laml.sign)

# Choose n where the correlation value (plot above) on the y-axis drops significantly.
laml.sig = extractSignatures(mat = laml.tnm, n = 4)

#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures")

#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")
pheatmap::pheatmap(mat = laml.v3.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures")

maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")

# 
#Visualize first signature
sig1 = laml.sig$signatures[,1]
barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)
