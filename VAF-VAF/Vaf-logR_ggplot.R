# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("IRanges")

library(VariantAnnotation)
library(tidyverse)
library(RColorBrewer)
# library(parallel)
# library(doParallel)
# library(foreach)
library(IRanges)
library(data.table)
library(tictoc)

setwd("~/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/")

snv_caller = "mutect2"
snv_caller = "mutect2-PASS"
# snv_caller = "varscan2"

cnv_caller = "sequenza"

output = "VAF_LogR_plots/"
dir.create(output, showWarnings = F)

# SNV callers
get_snv_data <- function(snv_caller, index){
    if (snv_caller=="mutect2"){
        # get vcf files
        vcfs = list.files("data/connor_pipeline/SNV/vcf_files/mutect2_matched_vcf", full.names = T, pattern = "filterFlag.vcf.gz$") # $ is end of string
        
        snv_data = get_mutect2_data(index, vcfs)
        
    } else if (snv_caller=="mutect2-PASS"){
        # get vcf files
        vcfs = list.files("data/connor_pipeline/SNV/vcf_files/mutect2_matched_vcf", full.names = T, pattern = "filterFlag.vcf.gz$") # $ is end of string
        
        snv_data = get_mutect2_data(index, vcfs)
        snv_data = filter(snv_data, pass=="Yes")
        
    } else if (snv_caller=="varscan2"){
        # get vcf files
        vcfs = list.files("vcf_files/varscan_somatic_joint_calling", full.names = T, pattern = "pon-filtered_strict.vcf.gz$") # $ is end of string
        
        snv_data = get_varscan2_data(index, vcfs)
    }
    
    return(snv_data)
}

get_mutect2_data <- function(index, vcfs){
    # get sample name
    sample = str_split(vcfs[index], pattern = "/|_")[[1]][10]
    print(sample)
    
    vcf = readVcf(vcfs[index])
    vcf = VariantAnnotation::expand(vcf) # expand ALT to multiple rows
    
    data = data.table(chr = as.character(seqnames(vcf)),
                      start = start(ranges(vcf)),
                      end = end(ranges(vcf)),
                      VAF = geno(vcf)$AF[,1],
                      REF_DP = geno(vcf)$AD[,2,1],
                      ALT_DP = geno(vcf)$AD[,1,1],
                      Total_DP = geno(vcf)$DP[,2],
                      pass = ifelse(rowRanges(vcf)$FILTER=="PASS", "Yes", "No"),
                      sample=sample)
    return(data)
}

get_varscan2_data <- function(index, vcfs){
    # get sample name
    sample = str_split(vcfs[index], pattern = "/|_")[[1]][7]
    print(sample)
    
    vcf = readVcf(vcfs[index])
    vcf = VariantAnnotation::expand(vcf) # expand ALT to multiple rows
    
    data = data.table(chr = as.character(seqnames(vcf)),
                      start = start(ranges(vcf)),
                      end = end(ranges(vcf)),
                      VAF = sapply(str_split(geno(vcf)$FREQ[,2],"%"), function(e) as.numeric(e[1])/100),
                      REF_DP = geno(vcf)$RD[,2],
                      ALT_DP = geno(vcf)$AD[,2],
                      Total_DP = geno(vcf)$DP[,2],
                      pass = ifelse(rowRanges(vcf)$FILTER=="PASS", "Yes", "No"),
                      sample=sample)
    return(data)
}

# CNV callers

get_cnv_data <- function(cnv_caller, index){
    if (cnv_caller=="ngcgh"){
        cnv_files = list.files("CNV/ngcgh", full.names = T, pattern = ".txt$")
        
        cnv_data = fread(cnv_files[index], header = F)
        names(cnv_data) = c("chr","start2","end2","normal_reads","tumor_reads","logR")
        cnv_data$end2 = cnv_data$end2 - 1
        
    } else if (cnv_caller=="sequenza"){
        cnv_files = list.files("data/connor_pipeline/CNV/Sequenza/sequenza_segments", full.names = T, pattern = ".txt$")
        cnv_data = fread(cnv_files[index], header = T)
        cnv_data$logR = log2(2 * cnv_data$depth.ratio) - 1
        cnv_data = dplyr::select(cnv_data, chr=chromosome, start2=start.pos, end2=end.pos, N.BAF, logR)
    }
    
    return(cnv_data)
}


get_combined_snv_cnv_data <- function(snv_caller, cnv_caller, index){
    # snv data
    snv_data = get_snv_data(snv_caller, index)
    # cnv data
    cnv_data = get_cnv_data(cnv_caller, index)
    
    # Set keys used for finding overlaps
    setkey(snv_data, chr, start, end)
    setkey(cnv_data, chr, start2, end2)
    
    data = foverlaps(snv_data,cnv_data) %>%
        dplyr::select(chr, start, end, VAF, logR, REF_DP, ALT_DP, Total_DP, pass, sample) %>% 
        distinct()
    
    return(data)
}

plot_vaf_logR <- function(data){
    sample = unique(data$sample)
    p = data %>%
        ggplot(aes(VAF, logR, color=ALT_DP)) +
        geom_point(size=2) +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(-1,1)) +
        labs(title=paste("Sample:", sample), color = "ALT Depth", ylab = "LogR") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    print(p)
    p
}


for(index in 1:length(vcfs)){
    data = get_combined_snv_cnv_data(snv_caller, cnv_caller, index)
    sample = unique(data$sample)
    
    TDP = 5; ADP = 3; AF = 0.01
    p = data %>% filter(Total_DP>=TDP, ALT_DP>=ADP, VAF>=AF) %>% plot_vaf_logR()
    output2 = paste0(output, "plot_VAF_LogR_",snv_caller,"-PASS_",cnv_caller,"-TumorDP_",TDP,"-TumorAD_",ADP,"-TumorAF_",AF,"/")
    dir.create(output2, showWarnings = F)
    ggsave(paste0(output2, sample, "-VAF_LogR_plot-",snv_caller,"_",cnv_caller,"-TumorDP_",TDP,"-TumorAD_",ADP,"-TumorAF_",AF,".png"), plot = p, width=8, height = 6)
}




# TESTING OVERLAP METHODS
# tic()
# data2 = inner_join(data, select(ngcgh, chr, start, end, logR), by = "chr") %>%
#     filter(pos >= start, pos < end) %>%
#     select(chr, pos, VAF, logR, pass)
# toc()
# 
# tic()
# ngcgh_df = as.data.frame(ngcgh)
# data3 = data %>% mutate(
#     logR=map_dbl(
#         .x = 1:nrow(.),
#         .f = ~ ngcgh_df[chr[.x] == ngcgh_df$chr & pos[.x]>=ngcgh_df$start & pos[.x]<=ngcgh_df$end,"logR"][1]
#     )
# )
# data9 = distinct(data3)
# toc()
# 
# tic()
# data4 = merge(data, select(ngcgh, chr, start, end, logR), by = "chr") %>%
#     filter(pos > start, pos < end)
# toc()
# 
# tic()
# d_range = GRanges(seqnames = data$chr, ranges=IRanges(data$pos, data$pos))
# n_range = GRanges(seqnames = ngcgh$chr, ranges = IRanges(ngcgh$start, ngcgh$end))
# ov = findOverlaps(d_range, n_range, type = "within")
# ov[["chr1"]]
# toc()
# 
# tic()
# data5 = setDT(data)
# data5[, pos2 := data$pos]
# setkey(data5, chr, pos, pos2)
# ngcgh_dt = setDT(ngcgh)
# setkey(ngcgh_dt, chr, start, end)
# data6 = foverlaps(data5,ngcgh_dt, nomatch = 0)
# toc()
