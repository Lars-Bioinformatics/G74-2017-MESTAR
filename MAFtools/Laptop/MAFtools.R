# Link to guide: https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html

# if (!require("BiocManager")) install.packages("BiocManager") 
# BiocManager::install("maftools")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# install.packages("mclust")
# install.packages("NMF")
# install.packages("pheatmap")
# install.packages("barplot3d")


library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
library(tidyverse)
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

workdir = "~/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/output/connor_pipeline/MAFtools_output/"
caller = "mutect2"

fix_sample_names_files <- function(){
    input = "../maf_files/mutect2_maf_PASSonly/"
    output = "../maf_files/mutect2_maf_PASSonly_fixedSampleNames/"
    dir.create(output, showWarnings = F)
    maf.files = list.files(input,full.names = T)
    maf.file = maf.files[1]
    # maf.simple = do.call(rbind, lapply(maf.files, read.table, header=T, sep="\t", fill = NA, quote = ""))
    for (maf.file in maf.files){
        maf = read.maf(maf = maf.file)
        maf.simple = maf@data
        # maf.simple = read.table(file = f, header = T, sep = "\t", file = NA, quote = "")
        maf.simple$Tumor_Sample_Barcode = sapply(maf.simple$Tumor_Sample_Barcode, function(s) str_split(s, pattern = "_")[[1]][1])
        maf.simple$Matched_Norm_Sample_Barcode = sapply(maf.simple$Matched_Norm_Sample_Barcode, function(s) str_split(s, pattern = "_")[[1]][1])
        tumor = unique(maf.simple$Tumor_Sample_Barcode)
        normal = unique(maf.simple$Matched_Norm_Sample_Barcode)
        out_file = paste0(tumor,"_vs_",normal,"_somatic_mutect2_filterFlag_PASSonly.maf")
        write.table(maf.simple, file = paste0(output,out_file), sep = "\t", quote = F, row.names = F)
    }
}
# fix_sample_names_files()

fix_sample_names_maf <- function(maf){
    maf@data$Tumor_Sample_Barcode = sapply(maf@data$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@data$Matched_Norm_Sample_Barcode = sapply(maf@data$Matched_Norm_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variants.per.sample$Tumor_Sample_Barcode = sapply(maf@variants.per.sample$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variant.type.summary$Tumor_Sample_Barcode = sapply(maf@variant.type.summary$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variant.classification.summary$Tumor_Sample_Barcode = sapply(maf@variant.classification.summary$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@maf.silent$Tumor_Sample_Barcode = sapply(maf@maf.silent$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@maf.silent$Matched_Norm_Sample_Barcode = sapply(maf@maf.silent$Matched_Norm_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@clinical.data$Tumor_Sample_Barcode = sapply(maf@clinical.data$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    return(maf)
}

# Sample combinations
sample_combinations = c(#"All samples" = ".maf",
                        #"Tumor samples" = "G74-[A-Z][1-9]",
                        # "Plasma samples" = "G74-[A-Z]P[1-9]",
                        # "G74-C" = "G74-C[1-9]",
                        # "G74-D" = "G74-D[1-9]",
                        # "G74-E" = "G74-E[1-9]",
                        # "G74-F" = "G74-F[1-9]",
                        # "G74-G" = "G74-G[1-9]",
                        # "G74-K" = "G74-K[1-9]",
                        # "G74-M" = "G74-M[1-9]",
                        # "G74-N" = "G74-N[1-9]",
                        "manuscript1" =  
                        )

# s_c = list(
# "All samples" = list.files("../maf_files/mutect2_maf_PASSonly", pattern = ".maf", full.names = T),
# "Tumor samples" = list.files("../maf_files/mutect2_maf_PASSonly", pattern = "G74-[A-Z][1-9]", full.names = T),
# "Plasma samples" = list.files("../maf_files/mutect2_maf_PASSonly", pattern = "G74-[A-Z]P[1-9]_", full.names = T)
# )
# 
# 
# list.files("../maf_files/mutect2_maf_PASSonly", pattern = "G74-[A-Z][1-9]", full.names = T)

for (sc in 1:length(sample_combinations)){
    # Reset working dir
    setwd(workdir)
    
    # Read MAF files
    # maf.files = list.files("../maf_files/mutect2_maf_PASSonly", pattern = sample_combinations[sc], full.names = T)
    # maf.files = list.files("../../../data/connor_pipeline/SNV/maf_files/mutect2_force_calling_mergedMatchedAndJointCalling/", pattern = sample_combinations[sc], full.names = T)
    maf.files = list.files("../../../data/connor_pipeline/SNV/maf_files/mutect2_ctDNA_pipeline_maf-no_twist", pattern = sample_combinations[sc], full.names = T)
    mutect2.maf = merge_mafs(maf.files)
    mutect2.maf@data$t_vaf = mutect2.maf@data$t_alt_count/mutect2.maf@data$t_depth
    mutect2.maf@maf.silent$t_vaf = mutect2.maf@maf.silent$t_alt_count/mutect2.maf@maf.silent$t_depth
    mutect2.maf@data$n_vaf = mutect2.maf@data$n_alt_count/mutect2.maf@data$n_depth
    mutect2.maf@maf.silent$n_vaf = mutect2.maf@maf.silent$n_alt_count/mutect2.maf@maf.silent$n_depth
    dim(mutect2.maf@data)
    
    t_alt_minReads=7
    n_maxVAF = 0.01
    n_ref_minReads=20
    mutect2.maf = subsetMaf(mutect2.maf, query = paste0("t_alt_count >= ",t_alt_minReads))
    mutect2.maf = subsetMaf(mutect2.maf, query = paste0("n_vaf <= ",n_maxVAF))
    mutect2.maf = subsetMaf(mutect2.maf, query = paste0("n_ref_count >= ",n_ref_minReads))
    dim(mutect2.maf@data)
    # mutect2.maf = fix_sample_names_maf(mutect2.maf)
    
    # Get sample group
    sample_group = names(sample_combinations[sc])
    
    # Create output dir and update work dir
    output = paste0("mutect2 - ", sample_group, " - t_alt_count_above",t_alt_minReads,"/")
    dir.create(output, showWarnings = F)
    setwd(paste0(workdir, output))

    # extract sample name
    tum_sample_barcodes = unique(mutect2.maf@data$Tumor_Sample_Barcode)
    
    #Shows sample summry - count of mutations
    getSampleSummary(mutect2.maf)
    #Shows gene summary.
    getGeneSummary(mutect2.maf)
    #shows clinical data associated with samples
    getClinicalData(mutect2.maf)
    #Shows all fields in MAF
    getFields(mutect2.maf)
    #Writes maf summary to an output file with basename mutect2.maf.
    write.mafSummary(maf = mutect2.maf, basename = '0. mutect2')
    
    pdf("1. Variant Summary.pdf", width = 16, height = 9)
    plotmafSummary(maf = mutect2.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    top_genes = 50
    pdf(paste0("2. oncoplot - Top ",top_genes," Genes.pdf"), width = 10, height = 12)
    # oncoplot(maf = mutect2.maf, top = 50, bgCol = "lightgrey")
    oncoplot(maf = mutect2.maf, top = top_genes, showTumorSampleBarcodes = T)
    dev.off()
    
    mutect2.maf.titv = titv(maf = mutect2.maf, plot = FALSE, useSyn = TRUE)
    pdf("3. Transitions and Transversions.pdf", width = 12, height = 10)
    #plot titv summary
    plotTiTv(res = mutect2.maf.titv, showBarcodes = T)
    dev.off()
    
    if ("TP53" %in% getGeneSummary(mutect2.maf)$Hugo_Symbol){
        pdf(paste0("4. Lollipop Plot of TP53.pdf"), width = 12, height = 8)
        # lollipopPlot(maf = mutect2.maf, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE)
        # lollipopPlot(maf = mutect2.maf, gene = 'TP53', showDomainLabel = FALSE, labelPos = 273)
        lollipopPlot(maf = mutect2.maf, gene = 'TP53', showDomainLabel = T, labelPos = "all", AACol = "HGVSp_Short")
        dev.off()
    }
    
    # Lollipop plot of top mutated gene
    top_mutated_gene = getGeneSummary(mutect2.maf) %>% arrange(desc(total)) %>% .[1,1]
    pdf(paste0("4. Lollipop Plot of Top Mutated Gene - ",top_mutated_gene,".pdf"), width = 12, height = 8)
    lollipopPlot(maf = mutect2.maf, gene = top_mutated_gene, showDomainLabel = T, labelPos = "all", AACol = "HGVSp_Short")
    dev.off()
    
    # for (maf.file in maf.files[9:10]){ # Only files with changes: maf.files[9:10]
    dir.create("5. Rainfall Plots", showWarnings = F)
    setwd(paste0(workdir,output,"5. Rainfall Plots"))
    for (sample in tum_sample_barcodes){ # Only files with changes: maf.files[9:10]
        try({
            pdf(paste0(sample, " rainfall plot.pdf"), width = 9, height = 6)
            rainfallPlot(maf = mutect2.maf, detectChangePoints = TRUE, pointSize = 0.4, tsb = sample)
            dev.off()
        })
    }
    setwd(paste0(workdir,output))
    
    pdf("6. Mutational load compared to TCGA cohorts.pdf", width = 9, height = 7)
    mutect2.maf.mutload = tcgaCompare(maf = mutect2.maf, cohortName = "MESTAR" , logscale = TRUE, capture_size = 50)
    dev.off()
    
    # mutect2.maf@data$t_vaf = mutect2.maf@data$t_alt_count/mutect2.maf@data$t_depth
    pdf(paste0("7. VAF per Gene - Top ", top_genes, " Genes.pdf"), width = 16, height = 7)
    plotVaf(maf = mutect2.maf, vafCol = 't_vaf', top = top_genes)
    # plotVaf(maf = mutect2.maf, vafCol = 't_vaf')
    dev.off()
    
    
    
    # 9 Analysis
    
    # 9.1 Somatic interactions
    # Mutually exclusive or co-occurring set of genes can be detected using somaticInteractions function, which performs pair-wise Fisher’s Exact test to detect such significant pair of genes.
    try({
        pdf("8. Somatic Interactions.pdf", width = 16, height = 16)
        interactions = somaticInteractions(maf = mutect2.maf, top = top_genes, pvalue = c(0.05, 0.1))
        dev.off()
        write.table(interactions, file = "8. Somatic Interactions List.txt", sep = "\t", quote = F, row.names = F)
    })
    
    # 9.2 Detecting cancer driver genes based on positional clustering
    # If less than 5 mutations in any gene, then set value to max seen, just to include at least one gene
    minMutRequired = min(max(getGeneSummary(mutect2.maf)$total), 4)
    mutect2.maf.sig = oncodrive(maf = mutect2.maf, AACol = 'HGVSp_Short', minMut = minMutRequired, pvalMethod = 'zscore')
    pdf("9. Cancer Driver Genes Analysis - Positional Clustering.pdf", width = 10, height = 8)
    plotOncodrive(res = mutect2.maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
    dev.off()
    write.table(mutect2.maf.sig, file = "9. Cancer Driver Genes List.txt", sep = "\t", quote = F, row.names = F)
    
    # 9.3 Adding and summarizing pfam domains
    top_genes_label = 25
    pdf("10. Pfam Domain Summary.pdf", width = 12, height = 8)
    mutect2.maf.pfam = pfamDomains(maf = mutect2.maf, AACol = 'HGVSp_Short', top = top_genes_label)
                                   # width = 12, height = 8, baseName = "10. Pfam")
    dev.off()
    # mutect2.maf.pfam = pfamDomains(maf = mutect2.maf, AACol = 'HGVSp', top = 10)
    
    #Protein summary (Printing first 7 columns for display convenience)
    mutect2.maf.pfam$proteinSummary[,1:7, with = FALSE]
    write.table(mutect2.maf.pfam$proteinSummary, file = "10. Protein Summary.txt", sep = "\t", quote = F, row.names = F)
    
    #Domain summary (Printing first 3 columns for display convenience)
    # mutect2.maf.pfam$domainSummary[,1:3, with = FALSE]
    write.table(mutect2.maf.pfam$domainSummary, file = "10. Domain Summary.txt", sep = "\t", quote = F, row.names = F)
    
    # 9.7 Drug-Gene Interactions
    pdf("11. Drug-Gene Interactions.pdf", width = 12, height = 8)
    dgi = drugInteractions(maf = mutect2.maf, fontSize = 0.75)
    dev.off()
    write.table(dgi, file = "11. Drug-Gene Interactions.txt", sep = "\t", quote = F, row.names = F)
    
    # 9.8 Oncogenic Signaling Pathways
    pdf("12. Oncogenic Signaling Pathways.pdf", width=10, height=8)
    pathways = OncogenicPathways(maf = mutect2.maf)
    dev.off()
    top_pathways = pathways %>% arrange(desc(Mutated_samples)) %>% slice_head(n = min(nrow(pathways), 2))
    for (pathway in top_pathways$Pathway){
        pdf(paste("12.",pathway,"- Oncogenic Signaling Pathways.pdf"), width=10, height=8)
        PlotOncogenicPathways(maf = mutect2.maf, pathways = pathway)
        dev.off()
    }
    # pdf("12. RTK-RAS - Oncogenic Signaling Pathways.pdf", width=10, height=8)
    # PlotOncogenicPathways(maf = mutect2.maf, pathways = "RTK-RAS")
    # dev.off()
    # pdf("12. TP53 - Oncogenic Signaling Pathways.pdf", width=10, height=8)
    # PlotOncogenicPathways(maf = mutect2.maf, pathways = "TP53")
    # dev.off()
    
    
    # 9.9 Tumor heterogeneity and MATH scores
    # 9.9.1 Heterogeneity in tumor samples
    dir.create("13. Tumor Heterogeneity and MATH scores", showWarnings = F)
    setwd(paste0(workdir,output,"13. Tumor Heterogeneity and MATH scores"))
    for (tsb in tum_sample_barcodes){
        pdf(paste("13. ", tsb, "- Tumor Heterogeneity and MATH score.pdf"), width=8, height=8)
        inferred.het = inferHeterogeneity(maf = mutect2.maf, tsb = tsb, vafCol = 't_vaf')
        plotClusters(clusters = inferred.het)
        dev.off()
    }
    setwd(paste0(workdir,output))
    
    # d = density(mutect2.maf@data$t_vaf)
    # stats::plot(d)
    
    
    # 9.10 Mutational Signatures
    # Only perform signature analysis if more than 10 samples
    if (length(tum_sample_barcodes)>15){
        mutect2.maf.tnm = trinucleotideMatrix(maf = mutect2.maf, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
        write.table(mutect2.maf.tnm, file = "14. Mutational Signature Analysis - mutation counts.txt", quote = F, sep = "\t", row.names = F)
        
        # 9.10.1 APOBEC Enrichment estimation
        # 9.10.2 Differences between APOBEC enriched and non-enriched samples
        # plotApobecDiff(tnm = mutect2.maf.tnm, maf = mutect2.maf, pVal = 0.2)
        
        
        # 9.10.3 Signature analysis
        # Note: If either extractSignatures or estimateSignatures stops in between, 
        #       its possibly due to low mutation counts in the matrix. In that case rerun
        #       the functions with pConstant argument set to small positive value (e.g, 0.1).
        mutect2.maf.sign = estimateSignatures(mat = mutect2.maf.tnm, nTry = 6)
        # mutect2.maf.sign = estimateSignatures(mat = mutect2.maf.tnm, nTry = 6, pConstant = 0.1)
        pdf("14. Mutational Signature Analysis - cophenetic metric.pdf", width = 8, height = 8)
        plotCophenetic(res = mutect2.maf.sign)
        dev.off()
        
        # Choose n where the correlation value (plot above) on the y-axis drops significantly.
        mutect2.maf.sig = extractSignatures(mat = mutect2.maf.tnm, n = 2)
        write.table(tibble::rownames_to_column(as.data.frame(mutect2.maf.sig$signatures), "Mutation"), 
                    file = "14. Mutational Signature Analysis - Novel identified signatures.txt", 
                    sep = "\t", quote = F, row.names = F)
        write.table(tibble::rownames_to_column(as.data.frame(mutect2.maf.sig$contributions), "Signature"), 
                    file = "14. Mutational Signature Analysis - signature contributions.txt", 
                    sep = "\t", quote = F, row.names = F)
        
        #Compate against original 30 signatures 
        mutect2.maf.og30.cosm = compareSignatures(nmfRes = mutect2.maf.sig, sig_db = "legacy")
        pdf("14. Mutational Signature Analysis - Cosine Similarity Against COSMIC Signatures V2.pdf", width = 14, height = 6)
        pheatmap::pheatmap(mat = mutect2.maf.og30.cosm$cosine_similarities, 
                           cluster_rows = FALSE, 
                           main = "Cosine Similarity Against COSMIC Signatures V2")
        dev.off()
        
        #Compate against updated version3 60 signatures 
        mutect2.maf.v3.cosm = compareSignatures(nmfRes = mutect2.maf.sig, sig_db = "SBS")
        pdf("14. Mutational Signature Analysis - Cosine Similarity Against COSMIC Signatures V3.pdf", width = 14, height = 6)
        pheatmap::pheatmap(mat = mutect2.maf.v3.cosm$cosine_similarities, 
                           cluster_rows = FALSE, 
                           main = "cosine similarity against validated signatures")
        dev.off()
        
        pdf("14. Mutational Signature Analysis - Found COSMIC Signatures V2.pdf", width = 10, height = 10)
        maftools::plotSignatures(nmfRes = mutect2.maf.sig, title_size = 1.2, sig_db = "legacy")
        dev.off()
        pdf("14. Mutational Signature Analysis - Found COSMIC Signatures V3.pdf", width = 10, height = 10)
        maftools::plotSignatures(nmfRes = mutect2.maf.sig, title_size = 1.2, sig_db = "SBS")
        dev.off()
        
        # 3d plot of signatures - Not too pretty
        #Visualize first signature
        # sig1 = mutect2.maf.sig$signatures[,1]
        # barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)
    }
    
    
    # 10.3 Prepare MAF file for MutSigCV analysis
    mutect2.maf.mutsig.corrected = prepareMutSig(maf = mutect2.maf)
    write.table(mutect2.maf.mutsig.corrected, file = "15. MutSigCV_Analysis_table.txt", sep = "\t", quote = F, row.names = F)

}
