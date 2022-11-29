# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("IRanges")

library(VariantAnnotation)
library(tidyverse)
library(maftools)
library(ggiraph)
library(htmlwidgets)
library(data.table)
# library(RColorBrewer)
library(parallel)
library(writexl)
library(vcfR)
library(hrbrthemes)
# library(doParallel)
# library(foreach)
# library(IRanges)
# library(tictoc)

setwd("~/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/")

# Select variant caller
snv_caller = "mutect2_joint_calling" # same as *_pass as they are filtered
snv_caller = "mutect2_joint_calling_pass"
snv_caller = "mutect2_joint_calling_no-PON"
snv_caller = "mutect2_force_calling_mergedMatchedAndJointCalling"
snv_caller = "mutect2_force_calling_mergedMatchedAndJointCalling_twist"
snv_caller = "mutect2_force_calling_mergedMatched"
snv_caller = "new_mutect2_joint_calling"
snv_caller = "mutect2_ctDNA_pipeline"
snv_caller = "mutect2_ctDNA_pipeline_noFFPE"
snv_caller = "mutect2"
snv_caller = "varscan2"
snv_caller = "varscan2_1000g-pon"
snv_caller = "varscan2_1000g-pon_q30"
snv_caller = "varscan2_1000g-pon_default-settings"
snv_caller = "varscan2_1000g-pon_gnomad"

# Select CNV caller
cnv_caller = "sequenza"

# Select pyclone run
pyclone_run = "pyclone_noPlasma_noFFPE" # Mutect2 analysis with FFPE (mutect2_ctDNA_pipeline), pyclone without FFPE input
pyclone_run = "pyclone" # Mutect2 analysis with FFPE (mutect2_ctDNA_pipeline), pyclone with FFPE input
pyclone_run = "pyclone_noPlasma_noFFPE_alsoNotInJointCalling" # Mutect2 analysis without FFPE (mutect2_ctDNA_noFFPE_pipeline), pyclone without FFPE input

# Choose output folder
# output = "output/connor_pipeline/VAF_VAF_plots2/"
output = "VAF_VAF_plots2/"

# run_vaf_vaf(snv_caller, cnv_caller)


# Run several snv_callers and cnv_callers
# for (snv_caller in c("mutect2_joint_calling_pass")){
#     for (cnv_caller in c("sequenza")){
#         # snv_path = get_snv_path(snv_caller)
#         # cnv_path = get_cnv_path(cnv_caller)
#         run_vaf_vaf(snv_caller, cnv_caller)
#     }
# }


###############################################################################
### Get data paths
###############################################################################
get_snv_path <- function(snv_caller){
    if (snv_caller=="mutect2") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_maf_PASSonly/"
    if (snv_caller=="mutect2_joint_calling") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_joint_calling_maf/"
    if (snv_caller=="mutect2_joint_calling_pass") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_joint_calling_PASS_maf/"
    if (snv_caller=="mutect2_joint_calling_no-PON") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_joint_calling_maf-no-Pon/"
    if (snv_caller=="mutect2_force_calling_mergedMatchedAndJointCalling") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_force_calling_mergedMatchedAndJointCalling/"
    if (snv_caller=="mutect2_force_calling_mergedMatchedAndJointCalling_twist") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_force_calling_mergedMatchedAndJointCalling_twist/"
    if (snv_caller=="mutect2_force_calling_mergedMatched") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_force_calling_mergedMatched/"
    if (snv_caller=="new_mutect2_joint_calling") snv_path = "data/connor_pipeline/SNV/maf_files/new_mutect2_joint_calling_maf/"
    if (snv_caller=="mutect2_ctDNA_pipeline") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_ctDNA_pipeline_maf/"
    # if (snv_caller=="mutect2_ctDNA_pipeline") snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_ctDNA_pipeline_noFFPE_maf/"
    
    if (snv_caller=="varscan2") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf/"
    if (snv_caller=="varscan2_1000g-pon") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered/"
    if (snv_caller=="varscan2_1000g-pon_q30") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered_q30/"
    if (snv_caller=="varscan2_1000g-pon_default-settings") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered_default-settings/"
    if (snv_caller=="varscan2_1000g-pon_gnomad") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g-pon_gnomad-filtered/"
    
    return(snv_path)
}

get_vcf_path <- function(snv_caller){
    if (snv_caller=="mutect2_ctDNA_pipeline") vcf_path = "data/connor_pipeline/SNV/vcf_files/mutect2_ctDNA_pipeline/vcf_matched/"
}

get_cnv_path <- function(cnv_caller){
    if (cnv_caller == "sequenza") cnv_path = "data/connor_pipeline/CNV/Sequenza/sequenza_segments/"
    return(cnv_path)
}

get_ploidy_path <- function(cnv_caller){
    if (cnv_caller == "sequenza") cnv_path = "data/connor_pipeline/CNV/Sequenza/sequenza_cellularity_ploidy/"
    return(cnv_path)
}

###############################################################################
### SNV data - MAF format
###############################################################################

# MAF Utils
fix_maf_sampleNames_and_add_vaf <- function(maf){
    maf@data$Tumor_Sample_Barcode = sapply(maf@data$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@data$Matched_Norm_Sample_Barcode = sapply(maf@data$Matched_Norm_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variants.per.sample$Tumor_Sample_Barcode = sapply(maf@variants.per.sample$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variant.type.summary$Tumor_Sample_Barcode = sapply(maf@variant.type.summary$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@variant.classification.summary$Tumor_Sample_Barcode = sapply(maf@variant.classification.summary$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@maf.silent$Tumor_Sample_Barcode = sapply(maf@maf.silent$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@maf.silent$Matched_Norm_Sample_Barcode = sapply(maf@maf.silent$Matched_Norm_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    maf@clinical.data$Tumor_Sample_Barcode = sapply(maf@clinical.data$Tumor_Sample_Barcode, function(s) str_split(s, "_")[[1]][1])
    
    maf@data$t_vaf = maf@data$t_alt_count/maf@data$t_depth
    maf@data$n_vaf = maf@data$n_alt_count/maf@data$n_depth
    
    maf@data$material = ifelse(str_detect(maf@data$Tumor_Sample_Barcode,"P"), "plasma", "tumor")
    
    return(maf)
}

# Read MAF data
read_mafs <- function(path = "maf_files/mutect2_maf_PASSonly", sample = ""){
    maf.files = list.files(path, pattern = ".maf", full.names = T)
    if (sample != "") maf.files = maf.files[str_detect(maf.files, sample)]
    maf = merge_mafs(maf.files)
    maf = fix_maf_sampleNames_and_add_vaf(maf)
}

read_mafs_dt <- function(path, skip_E=F){
    if (skip_E){ # Varscan
        maf.files = setdiff(list.files(path, full.names = T), list.files(path, pattern = "G74-E", full.names = T))
    } else {
        maf.files = list.files(path, pattern = ".maf", full.names = T)
    }
    data = rbindlist(lapply(maf.files, fread, sep="\t"))
    # data = parLapply(cl, maf.files, fread, sep="\t", select=c("Hugo_Symbol", "Chromosome", "Start_Position","End_Position", "Variant_Type","Tumor_Sample_Barcode",
    #                                                       "t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count"))
    data$t_vaf = data$t_alt_count/data$t_depth
    data$n_vaf = data$n_alt_count/data$n_depth
    
    data$material = ifelse(str_detect(data$Tumor_Sample_Barcode,"P"), "plasma", "tumor")
    return(data)
}

read_vaf_from_vcf <- function(vcf_path){
    vcf_files = list.files(vcf_path, pattern = "gz$", full.names = T)
    vcf_af = map_dfr(
        .x = vcf_files,
        .f = function(file){
            if (str_detect(basename(file), "_vs_")){
                vcf_sample = str_split(basename(file), "_vs_", simplify = T)[1]
            } else {
                vcf_sample = str_split(basename(file), "__", simplify = T)[1]
            }
            vcf_sample = str_remove(vcf_sample, "-twist")
            print(vcf_sample)
            sample = str_split(basename(file), "_", simplify = T)[1]
            print(sample)
            vcf = read.vcfR(file)
            # head(vcf)
            chr_pos = getFIX(vcf)[,1:2]
            # head(chr_pos)
            mutect2_vaf = extract.gt(vcf, element = "AF", as.numeric = T)
            # head(mutect2_vaf)
            setNames(data.frame(chr_pos[,1],as.numeric(chr_pos[,2]), unname(mutect2_vaf[,vcf_sample]), sample), c("chr","start","vcf_vaf","sample"))
        }
    )
}

# read_pyclone <- function(path){
#     path = "output/connor_pipeline/pyclone/G74-N_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY_pyclone/tables/"
#     
#     fread(input = )
# }

###############################################################################
### sequenza CNV data
###############################################################################
manually_calculate_ploidy <- function(data){
    # Tumor ploidy
    data2 = data %>% mutate(segment.size=stop-start)
    genome_lengths = data2 %>% group_by(sample) %>% summarise(genome_length = sum(segment.size))
    data2 = merge(data2, genome_lengths, by = "sample")
    # ploidy = length_of_segment * total_cnv_in_that segment / length_of_genome
    ploidy = data2 %>% group_by(sample) %>% summarise(ploidy = sum(segment.size*tumor.total/genome_length),
                                                      ploidy_rounded = round(sum(segment.size*tumor.total/genome_length))
                                                      )
    data = merge(data, ploidy, by = "sample")
    return(data)
}

read_sequenza <- function(path, ploidy_path="", manual_ploidy=F){
    cnv_files <- list.files(pattern = "segments", path = path)
    data <- map_dfr(
        .x = cnv_files,
        .f = function(file){
            sample = str_split(file, "_", simplify = T)[2]
            print(sample)
            data = read.table(paste0(path,"/",file), header = T) %>% dplyr::select(chr=chromosome, start=start.pos, stop=end.pos, tumor.total=CNt, tumor.major=A, tumor.minor=B)
            data$sample = sample
            data$sample_num = which(cnv_files==file)
            
            # Add ploidy
            if (manual_ploidy || ploidy_path==""){
                data = manually_calculate_ploidy(data)
            } else {
                ploidy_file = list.files(path = ploidy_path, pattern = paste0(sample,"_"), full.names = T)
                stats_data =  read.table(ploidy_file, header = T)
                ploidy = stats_data[2,2]
                purity = stats_data[2,1]
                if (is.na(ploidy)){
                    data = manually_calculate_ploidy(data)
                    purity = as.numeric(purity)
                    data$purity = purity
                } else {
                    ploidy = as.numeric(ploidy)
                    data$ploidy = ploidy
                    data$ploidy_rounded = round(ploidy)
                    purity = as.numeric(purity)
                    data$purity = purity
                }
                
            }
            data
        }
    )
    
    return(data)
}

###############################################################################
### Pyclone data
###############################################################################

pyclone_info <- function(pyclone_run = ""){
    if (pyclone_run == "pyclone"){
        pyclone_path = "output/connor_pipeline/pyclone_manual"
        pyclone_output = "pyclone"
    } else if(pyclone_run =="pyclone_noPlasma_noFFPE"){
        pyclone_path = "output/connor_pipeline/pyclone_noPlasma_noFFPE"
        pyclone_output = "pyclone_noPlasma_noFFPE"
    } else if(pyclone_run =="pyclone_noPlasma_noFFPE_alsoNotInJointCalling"){
        pyclone_path = "output/connor_pipeline/pyclone_noPlasma_noFFPE_alsoNotInJointCalling/"
        pyclone_output = "pyclone_noPlasma_noFFPE_alsoNotInJointCalling"
    }
    return(c(pyclone_path,pyclone_output))
}
read_pyclone_data <- function(path=""){
    # path = "output/connor_pipeline/pyclone_manual"
    files = list.files(path, "loci.tsv", recursive = T, full.names = T)
    data = rbindlist(lapply(files, fread, header=T))
    pyclone_data = data %>% 
        separate(mutation_id, c("chr","start","end"),sep="_", convert = T)
    return(pyclone_data)
}

manual_compute_ccf <- function(merged,vaf_col="t_vaf",ccf_col="ccf", cnv_col="tumor.total", purity_col="purity"){
    # VAF
    f = ifelse(is.na(merged[[vaf_col]]), 0, merged[[vaf_col]])
    # purity
    p = merged[[purity_col]]
    # Total copynumber
    Nt = ifelse(is.na(merged[[cnv_col]]), 2, merged[[cnv_col]])
    # multiplicity
    m = sapply(round( f/p * ( p*Nt + 2*(1-p)) ), function(v) max(1, v))
    # compute CCF
    ccf = f/m*p * (p*Nt + 2*(1-p))
    # return data frame
    merged[[ccf_col]] = ccf
    return(merged)
}

###############################################################################
### TitanCNA
###############################################################################
read_titan <- function(path="data/connor_pipeline/CNV/TitanCNA/optimalClusterSolution/"){
    params_files = list.files(path, pattern = "params")
    segs_files = list.files(path, pattern = "segs.txt")
    titan_df = map_dfr(
        .x = segs_files,
        .f = function(file){
            # file = segs_files[8]
            sample = str_split(file, "_", simplify = T)[1]
            print(sample)
            segs = read.table(paste0(path,"/",file), header = T)
            params_file = list.files(path = path, pattern = paste0(sample,"_.+params.txt$"), full.names = T)
            params = fread(params_file, nrows = 1)
            segs$titan_purity = as.numeric(params[1,2])
            segs$titan_ploidy = as.numeric(params[2,2])
            return(segs)
        }
    )
    titan_df = titan_df %>% dplyr::rename(chr=Chromosome, titan_start=Start_Position.bp., titan_end=End_Position.bp., sample=Sample,
                        titan_cnv=Corrected_Copy_Number, titan_majorCN=Corrected_MajorCN, titan_minorCN=Corrected_MinorCN,
                        titan_call=Corrected_Call)
    return(titan_df)
}
# titan_df = read_titan()

combine_snv_titan_data <- function(snv_data, titan_data){
    
    # snv_data = snv_data %>% dplyr::rename(sample=Tumor_Sample_Barcode, chr=Chromosome, start=Start_Position, end=End_Position)
    titan_data = titan_data %>% as.data.table()
    
    # Set keys used for finding overlaps
    setkey(snv_data, sample, chr, start, end)
    setkey(titan_data, sample, chr, titan_start, titan_end)
    
    combined = foverlaps(snv_data,titan_data) #%>%
        # mutate(sample_group=str_sub(sample, 1, 5))
    # c2=filter(combined, Hugo_Symbol=="DGKE")
    
    return(combined)
}

###############################################################################
### Combined functions for snv-cnv data object
###############################################################################
combine_snv_cnv_data <- function(snv_data, cnv_data){
    
    # snv_data = snv_data %>% dplyr::rename(sample=Tumor_Sample_Barcode, chr=Chromosome, start=Start_Position, end=End_Position)
    cnv_data = cnv_data %>% dplyr::rename(start2=start, end2=stop) %>% as.data.table()
    
    # Set keys used for finding overlaps
    setkey(snv_data, sample, chr, start, end)
    setkey(cnv_data, sample, chr, start2, end2)
    
    combined = foverlaps(snv_data,cnv_data) %>%
        mutate(sample_group=str_sub(sample, 1, 5))
    # c2=filter(combined, Hugo_Symbol=="DGKE")
    
    return(combined)
}

join_vaf_vs_vaf <- function(merged){
    data = merged %>%
        dplyr::select(sample, chr, start, end, t_vaf, t_alt_count, t_ref_count, n_ref_count, Hugo_Symbol, tumor.total, tumor.major, tumor.minor, sample_group, 
                      HGVSc, HGVSp, cluster_id, cellular_prevalence, ccf, vcf_vaf, vcf_ccf, ccf_titan, vcf_ccf_titan) %>%
        # mutate(sample_group=str_sub(sample, 1, 5)) %>%
        inner_join(., ., by=c("chr", "start", "end", "sample_group", "Hugo_Symbol", "HGVSc", "HGVSp","cluster_id")) %>%
        filter(sample.x < sample.y) %>%
        # mutate(sample_name_length_x=nchar(sample.x),sample_name_length_y=nchar(sample.y)) %>%
        rowwise() %>%
        mutate(tooltip=paste0("Gene: ", Hugo_Symbol,
                              "\nPosition: ", chr, ":", start,
                              "\nHGVSc: ", HGVSc,
                              "\nHGVSp: ", HGVSp,
                              # "\n", str_pad(sample.x,max(nchar(sample.y)-nchar(sample.x),0),),"&nbsp;&nbsp;CNV: ", tumor.total.x, "-", tumor.major.x, "-", tumor.minor.x, #"&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              # "\n", str_pad("", 12-sample_name_length_x), "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              "\n", sample.x, ifelse(nchar(sample.y)-nchar(sample.x)>0,"&nbsp;&nbsp;",""),"&nbsp;&nbsp;CNV: ", tumor.total.x, "-", tumor.major.x, "-", tumor.minor.x, "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              # "\n", paste0(rep("&nbsp;",max(nchar(sample.x),nchar(sample.y))), collapse = ""), "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              "\n", sample.y, ifelse(nchar(sample.x)-nchar(sample.y)>0,"&nbsp;&nbsp;",""),"&nbsp;&nbsp;CNV: ", tumor.total.y, "-", tumor.major.y, "-", tumor.minor.y, "&nbsp;&nbsp;SNV: ", t_alt_count.y, "-", t_ref_count.y
                              # "\n", paste0(rep("&nbsp;",max(nchar(sample.x),nchar(sample.y))), collapse = ""), "&nbsp;&nbsp;SNV: ", t_alt_count.y, "-", t_ref_count.y
               ))
    return(data)
}
# joined_data = join_vaf_vs_vaf(filtered)
# joined_data$tooltip

# add_pyclone_info <- function(joined_data, sample){
#     pyclone_info = read.table("output/connor_pipeline/pyclone/G74-C_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY_pyclone/plots", header = T)
# }

join_vaf_vs_vaf_keepAll <- function(merged){
    data = merged %>%
        dplyr::select(sample, chr, start, end, t_vaf, t_alt_count, t_ref_count, n_ref_count, Hugo_Symbol, tumor.total, tumor.major, tumor.minor, sample_group, 
                      HGVSc, HGVSp, cluster_id, cellular_prevalence, ccf, vcf_vaf, vcf_ccf, ccf_titan, vcf_ccf_titan) %>%
        # mutate(sample_group=str_sub(sample, 1, 5)) %>%
        inner_join(., ., by=c("chr", "start", "end", "sample_group", "Hugo_Symbol", "HGVSc", "HGVSp","cluster_id")) %>%
        # filter(sample.x < sample.y) %>%
        # mutate(sample_name_length_x=nchar(sample.x),sample_name_length_y=nchar(sample.y)) %>%
        rowwise() %>%
        mutate(tooltip=paste0("Gene: ", Hugo_Symbol,
                              "\nPosition: ", chr, ":", start,
                              "\nHGVSc: ", HGVSc,
                              "\nHGVSp: ", HGVSp,
                              # "\n", str_pad(sample.x,max(nchar(sample.y)-nchar(sample.x),0),),"&nbsp;&nbsp;CNV: ", tumor.total.x, "-", tumor.major.x, "-", tumor.minor.x, #"&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              # "\n", str_pad("", 12-sample_name_length_x), "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              "\n", sample.x, ifelse(nchar(sample.y)-nchar(sample.x)>0,"&nbsp;&nbsp;",""),"&nbsp;&nbsp;CNV: ", tumor.total.x, "-", tumor.major.x, "-", tumor.minor.x, "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              # "\n", paste0(rep("&nbsp;",max(nchar(sample.x),nchar(sample.y))), collapse = ""), "&nbsp;&nbsp;SNV: ", t_alt_count.x, "-", t_ref_count.x,
                              "\n", sample.y, ifelse(nchar(sample.x)-nchar(sample.y)>0,"&nbsp;&nbsp;",""),"&nbsp;&nbsp;CNV: ", tumor.total.y, "-", tumor.major.y, "-", tumor.minor.y, "&nbsp;&nbsp;SNV: ", t_alt_count.y, "-", t_ref_count.y
                              # "\n", paste0(rep("&nbsp;",max(nchar(sample.x),nchar(sample.y))), collapse = ""), "&nbsp;&nbsp;SNV: ", t_alt_count.y, "-", t_ref_count.y
        ))
    return(data)
}

filter_variants <- function(merged, filters){
    print(filters)
    # filters=c(DP=0, AD=0, RD=10, AF=0)
    # DP = as.integer(filters["DP"]); AD = as.integer(filters["AD"]); RD = as.integer(filters["RD"]); AF = as.integer(filters["AF"])
    t_minVAF=as.numeric(filters[1]); t_alt_minReads=as.numeric(filters[2]); n_maxVAF=as.numeric(filters[3]); n_ref_minReads=as.numeric(filters[4])
    
    
    extra = ""
    # if (DP>0) extra = paste0(extra, "-TumorTotalDepth_", DP)
    if (t_minVAF>0)        extra = paste0(extra, "-TumorAlleleFrequency_", t_minVAF)
    if (t_alt_minReads>0)  extra = paste0(extra, "-TumorAltDepth_", t_alt_minReads)
    if (n_maxVAF<1)        extra = paste0(extra, "-NormalAlleleFrequency_", n_maxVAF)
    if (n_ref_minReads>0)  extra = paste0(extra, "-NormalRefDepth_", n_ref_minReads)
    if (extra=="")         extra = "-No_Filters"
    extra <<- extra
    
    # t_minVAF = 0.1; t_alt_minReads = 10; n_maxVAF = 0.01; n_ref_minReads = 10
    # t_minVAF = 0; t_alt_minReads = 0; n_maxVAF = 1; n_ref_minReads = 0
    # Filter
    filtered_variants = merged %>%
        filter(material=="tumor") %>%
        group_by(chr, start, end, sample_group) %>% 
        summarise(t_alt_count=max(t_alt_count), t_vaf=max(t_vaf), n_vaf=max(n_vaf), n_ref_count=max(n_ref_count)) %>%
        filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads, n_vaf <= n_maxVAF, n_ref_count >= n_ref_minReads)
        # filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads)
    
    data = semi_join(merged, filtered_variants, by = c("chr", "start", "end", "sample_group"))
    return(data)
}



###############################################################################
### VAF-VAF plots
###############################################################################
vaf_vaf_per_sample <- function(sub_data, sample1, sample2, output_dir){
    # sub_data = data %>% filter(sample.x == sample1, sample.y==sample2)
    
    sample = sample_group=str_sub(sample1, 1, 5)
    output3 = paste0(output_dir, "VAF-VAF_",snv_caller,"-per_sample/patient_", sample,"/")
    dir.create(output3, showWarnings = F, recursive = T)
    
    p = sub_data %>%
        ggplot(aes(t_vaf.x, t_vaf.y)) +
        geom_point(size=3, color="blue",alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x=sample1, y=sample2) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    # p
    filename = paste0(output3, sample1, "_vs_", sample2, "_VAF-VAF_",snv_caller,".png")
    ggsave(filename, p, width = 8, height = 8)
}

vaf_vaf_per_sample_interactive <- function(sub_data, sample1, sample2, output_dir){
    # sub_data = data %>% filter(sample.x == sample1, sample.y==sample2)
    
    sample = sample_group=str_sub(sample1, 1, 5)
    output3 = paste0(output_dir, "VAF-VAF_",snv_caller,"-per_sample_interactive/patient_", sample,"/")
    dir.create(output3, showWarnings = F, recursive = T)
    
    p = sub_data %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes(x = t_vaf.x, y = t_vaf.y, tooltip = tooltip), size=3, color="blue",alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x=sample1, y=sample2) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    # girafe(ggobj = p)
    htmlwidgets::saveWidget(widget = girafe(ggobj = p), 
                            file = paste0(output3, sample1, "_vs_", sample2, "_VAF-VAF_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sample1, "_vs_", sample2, "_VAF-VAF_",snv_caller,"_interactive_files"), recursive = T)
}

vaf_vaf_per_sample_general <- function(sub_data, sample1, sample2, output_dir, points="VAF", col.x, col.y, color="cluster_id"){
    # sub_data = data %>% filter(sample.x == sample1, sample.y==sample2)
    
    sample = sample_group=str_sub(sample1, 1, 5)
    output3 = paste0(output_dir, "Per_sample/patient_", sample,"/")
    dir.create(output3, showWarnings = F, recursive = T)
    
    p = sub_data %>%
        ggplot(aes_string(col.x, col.y)) +
        geom_point(aes_string(color=color), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x=sample1, y=sample2, color="Pyclone Clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    # p
    filename = paste0(output3, sample1, "_vs_", sample2, "_",points,"-",points,"_",snv_caller,".pdf")
    ggsave(filename, p, width = 8, height = 8)
    return(p)
}

vaf_vaf_per_sample_interactive_general <- function(sub_data, sample1, sample2, output_dir, points="VAF", col.x, col.y, color="cluster_id"){
    # sub_data = data %>% filter(sample.x == sample1, sample.y==sample2)
    
    sample = sample_group=str_sub(sample1, 1, 5)
    output3 = paste0(output_dir,"Per_sample_interactive/patient_", sample,"/")
    dir.create(output3, showWarnings = F, recursive = T)
    
    p = sub_data %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes_string(x = col.x, y = col.y, tooltip = tooltip, color = color), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x=sample1, y=sample2, color="pyclone clustering") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    # girafe(ggobj = p)
    htmlwidgets::saveWidget(widget = girafe(ggobj = p), 
                            file = paste0(output3, sample1, "_vs_", sample2, "_",points,"-",points,"_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sample1, "_vs_", sample2, "_",points,"-",points,"_",snv_caller,"_interactive_files"), recursive = T)
}

vaf_vaf_combined_sample_general <- function(data, sg, output_dir, points="VAF", col.x, col.y){
    output3 = paste0(output_dir, "Combined_per_patient/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    p = data2 %>%
        ggplot(aes_string(col.x, col.y, color="cluster_id")) +
        geom_point(size=3,alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        # coord_cartesian(xlim=c(0,2), ylim=c(0,2)) +
        labs(x=paste(points,"(bottom name)"), y=paste(points,"(upper name)"), color="pyclone clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~sample.y + sample.x)
    # p
    filename = paste0(output3, sg, "_",points,"-",points,"_",snv_caller,".pdf")
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        ggsave(filename, p, width = 18, height = 6)
    } else if (n_samples<8){
        ggsave(filename, p, width = 18, height = 12)
    } else {
        ggsave(filename, p, width = 24, height = 24)
    }
}

vaf_vaf_combined_sample_interactive_general <- function(data, sg, output_dir, points="VAF", col.x, col.y){
    output3 = paste0(output_dir, "Combined_per_patient_interactive/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    
    p = data2 %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes_string(x = col.x, y = col.y, tooltip = "tooltip", color="cluster_id"), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x=paste(points,"(bottom name)"), y=paste(points,"(upper name)"), color="pyclone clustering") +
        # theme_minimal() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~sample.y + sample.x)
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    # if (n_samples<4){
    if (n_samples<4 || sg %in% c("G74-C")){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=6)
    } else if (n_samples<8){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=12)
    } else {
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=24, height_svg=24)
    }
    htmlwidgets::saveWidget(widget = girafe_plot, 
                            file = paste0(output3, sg, "_",points,"-",points,"_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sg, "_",points,"-",points,"_",snv_caller,"_interactive_files"), recursive = T)
    return(girafe(ggobj = p))
}

vaf_vaf_combined_sample <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "Combined_per_patient/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    p = data2 %>%
        ggplot(aes(t_vaf.x, t_vaf.y, color=cluster_id)) +
        geom_point(size=3,alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="VAF (bottom name)", y="VAF (upper name)", color="pyclone clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~sample.y + sample.x)
    # p
    filename = paste0(output3, sg, "_VAF-VAF_",snv_caller,".png")
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        ggsave(filename, p, width = 18, height = 6)
    } else if (n_samples<8){
        ggsave(filename, p, width = 18, height = 12)
    } else {
        ggsave(filename, p, width = 24, height = 24)
    }
}

vaf_vaf_combined_sample_interactive <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "Combined_per_patient_interactive/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    
    p = data2 %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes(x = t_vaf.x, y = t_vaf.y, tooltip = tooltip, color=cluster_id), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="VAF (bottom name)", y="VAF (upper name)", color="pyclone clustering") +
        # theme_minimal() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), ) +
        facet_wrap(~sample.y + sample.x)
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=6)
    } else if (n_samples<8){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=12)
    } else {
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=24, height_svg=24)
    }
    htmlwidgets::saveWidget(widget = girafe_plot, 
                            file = paste0(output3, sg, "_VAF-VAF_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sg, "_VAF-VAF_",snv_caller,"_interactive_files"), recursive = T)
    return(girafe(ggobj = p))
}

cp_cp_combined_sample <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "CP_vs_CP/Combined_per_patient/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    p = data2 %>%
        ggplot(aes(cellular_prevalence.x, cellular_prevalence.y, color=cluster_id)) +
        geom_point(size=3,alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="CP (bottom name)", y="CP (upper name)", color="pyclone clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~sample.y + sample.x)
    # p
    filename = paste0(output3, sg, "_CP-CP_",snv_caller,".png")
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        ggsave(filename, p, width = 18, height = 6)
    } else if (n_samples<8){
        ggsave(filename, p, width = 18, height = 12)
    } else {
        ggsave(filename, p, width = 24, height = 24)
    }
}

cp_cp_combined_sample_interactive <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "CP_vs_CP/Combined_per_patient_interactive/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    
    p = data2 %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes(x = cellular_prevalence.x, y = cellular_prevalence.y, tooltip = tooltip, color=cluster_id), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="CP (bottom name)", y="CP (upper name)", color="pyclone clustering") +
        # theme_minimal() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), ) +
        facet_wrap(~sample.y + sample.x)
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=6)
    } else if (n_samples<8){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=12)
    } else {
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=24, height_svg=24)
    }
    htmlwidgets::saveWidget(widget = girafe_plot, 
                            file = paste0(output3, sg, "_CP-CP_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sg, "_CP-CP_",snv_caller,"_interactive_files"), recursive = T)
    return(girafe(ggobj = p))
}

ccf_ccf_combined_sample <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "CCF_vs_CCF/Combined_per_patient/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    p = data2 %>%
        ggplot(aes(ccf.x, ccf.y, color=cluster_id)) +
        geom_point(size=3,alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="CFF (bottom name)", y="CFF (upper name)", color="pyclone clustering") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~sample.y + sample.x)
    # p
    filename = paste0(output3, sg, "_CCF-CCF_",snv_caller,".png")
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4){
        ggsave(filename, p, width = 18, height = 6)
    } else if (n_samples<8){
        ggsave(filename, p, width = 18, height = 12)
    } else {
        ggsave(filename, p, width = 24, height = 24)
    }
}

ccf_ccf_combined_sample_interactive <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "CCF_vs_CCF/Combined_per_patient_interactive/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    
    p = data2 %>%
        ggplot() +
        # geom_point_interactive(aes(t_vaf.x, t_vaf.y, tooltip=tooltip, data_id=tooltip))+#, size=3, color="blue",alpha=0.3) +
        geom_point_interactive(aes(x = ccf.x, y = ccf.y, tooltip = tooltip, data_id=tooltip, color=cluster_id), size=3, alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="CFF (bottom name)", y="CFF (upper name)", color="pyclone clustering") +
        # theme_minimal() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), ) +
        facet_wrap(~sample.y + sample.x)
    n_samples = length(unique(c(data2$sample.x, data2$sample.y)))
    if (n_samples<4 || sg %in% c("G74-C")){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=6)
    } else if (n_samples<8){
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=18, height_svg=12)
    } else {
        girafe_plot = ggiraph(ggobj = p, width = 1, zoom_max = 1, tooltip_opacity = 0.75, width_svg=24, height_svg=24)
    }
    htmlwidgets::saveWidget(widget = girafe_plot, 
                            file = paste0(output3, sg, "_CFF-CFF_",snv_caller,"_interactive.html"), 
                            selfcontained = T, libdir = NULL)
    unlink(paste0(output3, sg, "_CFF-CFF_",snv_caller,"_interactive_files"), recursive = T)
    return(girafe(ggobj = p))
}


optimizeParameters <- function(data){
    # Check normal vaf distribution
    ggplot(data, aes(x=n_vaf)) + geom_density()
    
    # Visualise relationships with filtering on VAF
    filtered_variants = map_dfr(
        .x = seq(0.1,1,0.1),
        .f = function(t_minVAF){
            print(t_minVAF)
            filtered_variants = data %>% 
                filter(!is.na(t_vaf)) %>%
                group_by(chr, start, end, sample_group) %>% 
                summarise(t_alt_count=max(t_alt_count), t_vaf=max(t_vaf), n_vaf=max(n_vaf), t_ref_count=max(t_ref_count),
                          n_ref_count=max(n_ref_count), t_depth=max(t_depth), n_depth = max(n_depth), .groups = "drop") %>%
                # filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads, n_vaf <= n_maxVAF, n_ref_count >= n_ref_minReads)
                # filter(t_vaf >= t_minVAF) %>%
                filter(t_vaf <= t_minVAF) %>%
                summarise(mean_t_alt_count = mean(t_alt_count), mean_t_vaf = mean(t_vaf), mean_t_ref_count=mean(t_ref_count),
                          mean_n_vaf = mean(n_vaf), mean_n_ref_count = mean(n_ref_count),
                          mean_t_depth = mean(t_depth), mean_n_depth = mean(n_depth))
            filtered_variants$t_minVAF = t_minVAF
            return(filtered_variants)
    })
    
    ggplot(filtered_variants) +
        # geom_line(aes(t_minVAF, mean_t_alt_count), color = "blue") +
        geom_line(aes(t_minVAF, mean_n_ref_count), color = "red") +
        # geom_line(aes(t_minVAF, mean_t_ref_count), color = "green") +
        # coord_cartesian(ylim=c(0,350))
        coord_cartesian(xlim=c(0,1))
        
    ggplot(filtered_variants) +
        geom_line(aes(t_minVAF, mean_t_depth), color = "blue") +
        geom_line(aes(t_minVAF, mean_n_depth), color = "red")
    
    
    # Visualise relationships with filtering on total tumor depth
    filtered_variants = map_dfr(
        .x = seq(5,100,5),
        .f = function(t_minDP){
            print(t_minDP)
            filtered_variants = data %>% 
                filter(!is.na(t_vaf)) %>%
                group_by(chr, start, end, sample_group) %>% 
                summarise(t_alt_count=max(t_alt_count), t_vaf=max(t_vaf), n_vaf=max(n_vaf), t_ref_count=max(t_ref_count),
                          n_ref_count=max(n_ref_count), t_depth=max(t_depth), n_depth = max(n_depth), .groups = "drop") %>% 
                # dim
                # filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads, n_vaf <= n_maxVAF, n_ref_count >= n_ref_minReads)
                # filter(t_alt_count >= t_minDP) %>%
                # filter(t_depth >= t_minDP) %>%
                # filter(t_alt_count <= t_minDP) %>%
                # filter(n_ref_count <= t_minDP & t_alt_count <= t_minDP) %>%
                filter(n_ref_count <= t_minDP & t_alt_count <= 10) %>%
                # filter(n_ref_count >= t_minDP | t_alt_count >= t_minDP) %>%
                summarise(mean_t_alt_count = mean(t_alt_count), mean_t_vaf = mean(t_vaf), mean_t_ref_count=mean(t_ref_count),
                          mean_n_vaf = mean(n_vaf), mean_n_ref_count = mean(n_ref_count),
                          mean_t_depth = mean(t_depth), mean_n_depth = mean(n_depth), num_variants=n())
            filtered_variants$t_minDP = t_minDP
            return(filtered_variants)
        })
    
    ggplot(filtered_variants) +
        geom_line(aes(t_minDP, mean_t_alt_count), color = "blue") +
        geom_line(aes(t_minDP, mean_n_ref_count), color = "red") +
        # geom_line(aes(t_minDP, mean_t_ref_count), color = "green") +
        ylab("Read Counts")
    
    ggplot(filtered_variants) +
        geom_line(aes(t_minDP, mean_t_depth), color = "blue") +
        geom_line(aes(t_minDP, mean_n_depth), color = "red")
    
    ggplot(filtered_variants) +
        geom_line(aes(t_minDP, mean_t_vaf), color = "blue") +
        geom_line(aes(t_minDP, mean_n_vaf), color = "red") +
        ylab("VAF")
    
    
    
    joined_data = join_vaf_vs_vaf(data)
    sum(is.na(joined_data$t_vaf.x))
    sd(joined_data$t_vaf.x, na.rm = T)
    sd(joined_data$t_vaf.y, na.rm = T)
    
    for (alt_count in seq(2,30,2)){
        d = joined_data %>% filter(t_alt_count.x>alt_count)
        print(sd(d$t_vaf.x, na.rm = T))
    }
    
    for (alt_count in seq(2,30,2)){
        for (n_count in seq(10,50,5)){
            d = filter(data, t_alt_count>alt_count, n_ref_count>n_count)
            print(paste(alt_count,n_count,sd(d$t_vaf, na.rm = T)))
        }
    }
    
    
}


get_combined_data = function(snv_caller, cnv_caller){
        
    snv_path = get_snv_path(snv_caller)
    vcf_path = get_vcf_path(snv_caller)
    cnv_path = get_cnv_path(cnv_caller)
    ploidy_path = get_ploidy_path(cnv_caller)
    
    # get snv data
    if (str_detect(snv_caller,"varscan")){
        # Exclude G74-E samples as they are huge!
        maf = read_mafs_dt(path = snv_path, skip_E = T)
    } else {
        # maf = read_mafs(path = snv_path)
        maf = read_mafs_dt(path = snv_path)
    }
    
    # avoid long names
    snv_caller2 = snv_caller
    snv_caller2 <<- snv_caller
    if (str_detect(snv_caller,"mutect2_force_calling")){
        snv_caller <- "mutect2_force_calling"
        snv_caller <<- "mutect2_force_calling"
    }
    
    dim(maf)
    # snv_data = filter(maf, FILTER == "PASS", !is.na(t_alt_count)) # PASS is no longer valid in force-calling
    # snv_data = filter(maf, !is.na(t_alt_count))
    snv_data = maf %>% dplyr::rename(chr=Chromosome, start=Start_Position, end=End_Position, sample=Tumor_Sample_Barcode)
    dim(snv_data)
    
    # get mutect2 vaf
    mutect2_vaf = read_vaf_from_vcf(vcf_path)
    snv_data = inner_join(snv_data, mutect2_vaf)#, by.x=c("chr","start","sample"), by.y=c("chr","pos","sample"))
    
    # get cnv data
    cnv_data = read_sequenza(path = cnv_path, ploidy_path = ploidy_path)
    purity_data = cnv_data %>% distinct(sample, purity, ploidy)
    dir.create("output/connor_pipeline/purity_ploidy_estimates/", showWarnings = F, recursive = T)
    write.table(x = purity_data, file = "output/connor_pipeline/purity_ploidy_estimates/All_samples_sequenza_purity_ploidy.txt", quote = F, row.names = F, sep = "\t")
    
    # get pyclone data
    # pyclone_data = read_pyclone_data("output/connor_pipeline/pyclone_manual")
    # pyclone_data = read_pyclone_data("output/connor_pipeline/pyclone_noPlasma_noFFPE/")
    pyclone_data = read_pyclone_data(pyclone_info(pyclone_run)[1])
    
    # Merge snv, cnv and pyclone data
    merged = combine_snv_cnv_data(snv_data = snv_data, cnv_data = cnv_data)
    merged = left_join(merged, dplyr::select(pyclone_data, chr, start, end, sample=sample_id, cluster_id, cellular_prevalence)) %>%
        mutate(cluster_id = ifelse(is.na(cluster_id), "Unclustered", cluster_id)) %>%
        dplyr::rename(cluster_id_orig=cluster_id)
    # Fix unclustered variants if defined in another sample
    clustered = merged %>% 
        dplyr::select(chr, start, end, sample_group, cluster_id=cluster_id_orig) %>%
        filter(cluster_id != "Unclustered") %>%
        distinct()
    merged = left_join(merged, clustered) %>% mutate(cluster_id = ifelse(is.na(cluster_id), "Unclustered", cluster_id))
    merged <<- merged
    table(merged$cluster_id)
    table(merged$cluster_id_orig)
    
    # Add sequenza purity
    merged = left_join(dplyr::select(merged,-purity), purity_data)
    
    # Add TitanCNA purity
    titan_df = read_titan()
    titan_purity = titan_df %>% distinct(sample, titan_purity, titan_ploidy)
    write.table(x = titan_purity, file = "output/connor_pipeline/purity_ploidy_estimates/All_samples_titanCNA_purity_ploidy.txt", quote = F, row.names = F, sep = "\t")
    titan_data = dplyr::select(titan_df, chr, titan_start, titan_end, sample, titan_cnv, titan_call, titan_purity, titan_ploidy)
    merged = combine_snv_titan_data(merged, titan_data)
    
    # compute CCF
    merged = manual_compute_ccf(merged, vaf_col = "t_vaf", ccf_col="ccf")
    merged = manual_compute_ccf(merged, vaf_col = "vcf_vaf", ccf_col="vcf_ccf")
    merged = manual_compute_ccf(merged, vaf_col = "t_vaf", ccf_col="ccf_titan", cnv_col = "titan_cnv", purity_col = "titan_purity")
    merged = manual_compute_ccf(merged, vaf_col = "vcf_vaf", ccf_col="vcf_ccf_titan", cnv_col = "titan_cnv", purity_col = "titan_purity")
    
    # maf@data %>% dplyr::count(Tumor_Sample_Barcode)
    merged %>% dplyr::count(chr,start,end,sample_group,Variant_Type) %>% arrange(desc(n)) %>% filter(n>0) %>% dplyr::count(sample_group, Variant_Type)
    
    # merged = filter(merged, sample_group %in% c("G74-D"))
    
    return(merged)
}

get_filtered_variants <- function(){
    merged = get_combined_data(snv_caller,cnv_caller)
    
    filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(7), n_maxVAF = 0.01, n_ref_minReads=c(20))
    filters = filters_config[1,]
    
    filtered = merged
    # filtered = filter(merged, !str_detect(sample, "twist"))
    # filtered = filter(merged, !str_detect(sample, "twist"), !(sample %in% c("G74-C1","G74-E1")))
    # filtered = filter(merged, !(sample %in% c("G74-C1","G74-E1")))
    # filtered = filter(merged, !str_detect(sample, "twist"), cluster_id != "unclustered")
    filtered = filter_variants(filtered, filters)
    # filtered = filter_variants(merged, filters)
    table(filtered$cluster_id)
    
    # output2 = paste0(output, "pyclone-", snv_caller2, "/", snv_caller ,extra,"/")
    output2 = paste0(output, pyclone_info(pyclone_run)[2],"-", snv_caller2, "/", snv_caller ,extra,"/")
    # output2 = paste0(output, pyclone_info(pyclone_run)[2],"-", snv_caller2, "_noPlasma_noFFPE/", snv_caller ,extra,"/")
    
    excel = filtered %>% dplyr::select(sample, gene=Hugo_Symbol, chr, start, end, t_vaf, t_depth, t_ref_count, t_alt_count, 
                               Variant_Classification, Variant_Type, ref=Reference_Allele, alt=Tumor_Seq_Allele2, 
                               HGVSc, HGVSp, Exon_Number, 
                               # start_cnv=start2, end_cnv=end2, 
                               copynumber=tumor.total, cnv.major=tumor.major, cnv.minor=tumor.minor, ploidy, sample_group)
    
    for (s in unique(excel$sample_group)){
    # for (s in c("G74-C")){
        excel_wide = excel %>% 
            filter(sample_group == s) %>% 
            pivot_wider(names_from = sample, values_from = c(t_vaf,t_depth,t_ref_count,t_alt_count,
                                                             copynumber,cnv.major,cnv.minor,ploidy))
        dir.create(paste0(output2,"filtered_variants/"), showWarnings = F, recursive = T)
        write.table(excel_wide, file = paste0(output2,"filtered_variants/", s, "_filtered_variants.txt"), quote = F, row.names = F, sep = "\t")
        # write.table(excel_wide, file = paste0("output/connor_pipeline/Data tables - Excel/manuscript 1 - variant select/", s, "_variant-select.txt"), quote = F, row.names = F, sep = "\t")
        # write_xlsx(excel_wide, path = paste0("output/connor_pipeline/Data tables - Excel/manuscript 1 - variant select/", s, "_variant-select.xlsx"))
    }
    
    # filtered %>% dplyr::select()
    # write.table(filtered, file = "output/connor_pipeline/Data tables - Excel/mutect2_snv_cnv_data.txt", sep = "\t", row.names = F, quote = F)
    
    # Write filtered data to file without tooltip column
    
    # Return
    return(filtered)
}
    
get_joined_data <-function(){
    filtered = get_filtered_variants()
    
    joined_data = join_vaf_vs_vaf(filtered)
    # joined_data = joined_data %>% mutate(cluster_color=ifelse(cluster_id.x=="unclustered",cluster_id.y,cluster_id.x))
    
    for (s in unique(joined_data$sample_group)){
        data2 = filter(joined_data, sample_group == s) %>% dplyr::select(-tooltip)
        dir.create(path = paste0(output2,"joined_data/"), recursive = T, showWarnings = F)
        write.table(data2, file = paste0(output2,"joined_data/",s,"_joined_data.txt"), quote = F, row.names = F, sep = "\t")
    }
    
    return(joined_data)

}

plot_vaf_vaf <- function(){
    joined_data = get_joined_data()
    as.data.table(joined_data)
    
    output2 = paste0(output, pyclone_info(pyclone_run)[2],"-", snv_caller2, "/", snv_caller ,extra,"/")
    
    configs = list(
        c("VAF","t_vaf"),
        c("CCF","ccf"),
        c("VCF_VAF", "vcf_vaf"),
        c("VCF_CCF", "vcf_ccf"),
        c("CP","cellular_prevalence"),
        c("CCF_TitanCNA","ccf_titan"),
        c("VCF_CCF_TitanCNA","vcf_ccf_titan")
    )
    
    for (config in configs){
        print(config)
        for (sg in unique(joined_data$sample_group)){
            print(sg)
            vaf_vaf_combined_sample_general(joined_data, sg, paste0(output2, config[1],"_vs_",config[1],"/"), points = config[1], col.x = paste0(config[2],".x"), col.y = paste0(config[2],".y"))
            # vaf_vaf_combined_sample_interactive_general(joined_data, sg, paste0(output2, config[1],"_vs_",config[1],"/"), points = config[1], col.x = paste0(config[2],".x"), col.y = paste0(config[2],".y"))
        }
    }
    
    pairs = joined_data %>% dplyr::select(sample.x, sample.y) %>% distinct() %>% as.data.frame(., stringsAsFactors=F)
    for (config in configs){
        print(config)
        for (i in 1:nrow(pairs)){
            print(paste(pairs[i,1], pairs[i,2]))
            # sub_data = joined_data %>% filter(sample.x == as.character(pairs[i,1]), sample.y==as.character(pairs[i,2]))
            sub_data = joined_data %>% filter(sample.x == pairs[i,1], sample.y==pairs[i,2])
            vaf_vaf_per_sample_general(sub_data, pairs[i,1], pairs[i,2], output_dir = paste0(output2, config[1],"_vs_",config[1],"/"),
                               points = config[1], col.x = paste0(config[2],".x"), col.y = paste0(config[2],".y"), color="cluster_id")
            # vaf_vaf_per_sample_interactive(sub_data, pairs[i,1], pairs[i,2], output2)
        }
    }
    
    # # plot VAF-VAF all samples per patient
    # for (sg in unique(joined_data$sample_group)){
    #     print(sg)
    #     vaf_vaf_combined_sample_general(joined_data, sg, paste0(output2, "VAF_vs_VAF/"), points = "VAF", col.x = "t_vaf.x", col.y = "t_vaf.y")
    #     vaf_vaf_combined_sample_interactive_general(joined_data, sg, paste0(output2, "VAF_vs_VAF/"), points = "VAF", col.x = "t_vaf.x", col.y = "t_vaf.y")
    # }
    # 
    # # plot CP-CP all samples per patient
    # for (sg in unique(joined_data$sample_group)){
    #     print(sg)
    #     cp_cp_combined_sample(joined_data, sg, output2)
    #     cp_cp_combined_sample_interactive(joined_data, sg, output2)
    # }
    # 
    # # plot CCF-CCF all samples per patient
    # for (sg in unique(joined_data$sample_group)){
    #     print(sg)
    #     ccf_ccf_combined_sample(joined_data, sg, output2)
    #     ccf_ccf_combined_sample_interactive(joined_data, sg, output2)
    # }
    
    
}

plot_vaf_vaf_manuscript <- function(){
    library(patchwork)
    
    filtered = get_filtered_variants()
    joined_data = join_vaf_vs_vaf_keepAll(filtered)
    dim(joined_data)
    joined_data = joined_data %>% mutate(sample.x=str_replace(sample.x,"G74-",""),
                                         sample.y=str_replace(sample.y,"G74-",""))
    
    output2 = paste0(output, pyclone_info(pyclone_run)[2],"-", snv_caller2, "/manuscript_figures/")
    
    config = c("VAF","t_vaf")
    
    draw_combined_plot <- function(pairs, patient=NULL, plots_combined=T, ncols=3, nrows=1, noLegend=F){
        # pairs = data.frame(sample.x = c("C6","C6","C7"),sample.y = c("C7","CP1","CP1"))
        if (is.null(patient)){
            patient = str_sub(pairs$sample.x[1], 1,1)
        }
        plots = list()
        for (i in 1:nrow(pairs)){
            print(paste(pairs[i,1], pairs[i,2]))
            sub_data = joined_data %>% filter(sample.x == pairs[i,1], sample.y==pairs[i,2])
            plots[[i]] = vaf_vaf_per_sample_general(sub_data, pairs[i,1], pairs[i,2], output_dir = paste0(output2, config[1],"_vs_",config[1],"/"),
                                                    points = config[1], col.x = paste0(config[2],".x"), col.y = paste0(config[2],".y"), color="cluster_id")
        }
        if (noLegend){
            p = plots[[1]] + theme(legend.position = "none")
            for (i in 2:length(plots)){
                p = p + (plots[[i]] + theme(legend.position = "none"))
            }
        } else {
            p = plots[[1]]
            for (i in 2:length(plots)){
                p = p + plots[[i]]
            }
        }
        p = p + plot_layout(guides = "collect", ncol = ncols, nrow = nrows)
        # p = plots[[1]] + plots[[2]] + plots[[3]] + plot_layout(guides = "collect", ncol = ncols, nrow = nrows)
        # p
        ggsave(paste0(output2,"combined/Patient_",patient,".pdf"), p, width = 19.75, height = 6)
        if(plots_combined){
            return(p)
        } else {
            return(plots)
        }
    }
    # Patient C
    C = draw_combined_plot(data.frame(sample.x = c("C6","C6","C7"),sample.y = c("C7","CP1","CP1")))
    # Patient E
    E = draw_combined_plot(data.frame(sample.x = c("E4","E4","E12"),sample.y = c("E12","EP1","EP1")))
    # Patient N
    N = draw_combined_plot(data.frame(sample.x = c("N5","N5","N8"),sample.y = c("N8","NP1","NP1")))
    
    # p = C[[1]] + C[[2]] + C[[3]] + E[[1]] + E[[2]] + E[[3]] + N[[1]] + N[[2]] + N[[3]] + plot_layout(ncol=3, nrow = 3)
    p = C / E / N
    p
    ggsave(paste0(output2,"combined/Patients_CEN.pdf"), p, width = 19.75, height = 18, scale=0.8)
    
    # Patient D
    D = draw_combined_plot(data.frame(sample.x = c("D2","D2","D3"),sample.y = c("D3","DP1","DP1")))
    
    
    # Patient F
    F_patient = draw_combined_plot(data.frame(sample.x = c("F2","F2","F6"),sample.y = c("F6","FP1","FP1")))
    # Patient K
    K_patient = draw_combined_plot(data.frame(sample.x = c("K4","K4","K8"),sample.y = c("K8","KP1","KP1")))
    p = F_patient / K_patient
    p
    ggsave(paste0(output2,"combined/Patients_FK.pdf"), p, width = 19.75, height = 12, scale=0.8)

    
    # Patient G
    G = draw_combined_plot(data.frame(sample.x = c("G1","G1","G2"),sample.y = c("G2","GP1","GP1")))
    # Patient M
    M = draw_combined_plot(data.frame(sample.x = c("M2","M2","M6"),sample.y = c("M6","MP1","MP1")))
    p = G / M
    p
    ggsave(paste0(output2,"combined/Patients_GM.pdf"), p, width = 19.75, height = 12, scale=0.8)
    
    # Plasma plasma VAF-VAF
    plasma_vaf = draw_combined_plot(data.frame(sample.x = c("CP1","DP1","EP1","FP1","GP1","KP1","MP1","NP1"),
                                               sample.y = paste0(c("CP1","DP1","EP1","FP1","GP1","KP1","MP1","NP1"),"-twist")),
                                    ncols = 3, nrows = 3, patient="plasma", noLegend = T)
    # plasma_vaf=p
    plasma_vaf = plasma_vaf + plot_layout(guides="auto")
    ggsave(paste0(output2,"combined/plasma-nimblegen_vs_plasma-twist.pdf"), plasma_vaf, width = 19.75, height = 18, scale=0.8)
    ggsave(paste0(output2,"combined/plasma-nimblegen_vs_plasma-twist-noLegends.pdf"), plasma_vaf, width = 18, height = 18, scale=0.8)
    # ggsave(paste0(output2,"combined/plasma-nimblegen_vs_plasma-twist.pdf"), plasma_vaf, width = 19.75, height = 22, scale=0.8)
    
    # Manuscript 1 - Fig 2 - patient E
    E1 = draw_combined_plot(data.frame(sample.x = c("E4","E4","E5"),sample.y = c("E5","E6","E6")),
                            ncols = 3, nrows = 1, patient = "Fig2-D1", noLegend = T)
    ggsave(paste0(output2,"combined/Fig2-D1.pdf"), E1, width = 18, height = 6, scale=0.6)
    E2 = draw_combined_plot(data.frame(sample.x = rep(c("E4","E5","E6"),each=6),
                                       sample.y = rep(c("E8","E9","E10","E11","E12","E13"),3)),
                            ncols = 9, nrows = 2, patient = "Fig2-D2", noLegend = T)
    ggsave(paste0(output2,"combined/Fig2-D2.pdf"), E2, width = 54, height = 12, scale=0.5)
    E3 = draw_combined_plot(data.frame(sample.x = c("E8","E8","E8","E8","E8","E9","E9","E9","E9","E10","E10","E10","E11","E11","E12"),
                                       sample.y = c("E9","E10","E11","E12","E13","E10","E11","E12","E13","E11","E12","E13","E12","E13","E13")),
                            ncols = 8, nrows = 2, patient = "Fig2-D3", noLegend = T)
    ggsave(paste0(output2,"combined/Fig2-D3.pdf"), E3, width = 48, height = 12, scale=0.5)
    # layout = "
    # AAA######
    # BBBBBBBBB
    # BBBBBBBBB
    # CCCCCCCC#
    # CCCCCCC##
    # "
    p = (E1 / E2 / E3) #+ plot_layout(ncol = 9, nrow=5)
    p
    ggsave(paste0(output2,"combined/M1Fig2.pdf"), p, width = 28, height = 22, scale=0.8)
}

plot_heatmap <- function(){
    library(patchwork)
    library(wesanderson) # colors
    library(yarrr) # colors
    
    # filtered = get_filtered_variants()
    merged = get_combined_data(snv_caller,cnv_caller)
    
    filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(7), n_maxVAF = 0.01, n_ref_minReads=c(20))
    filters = filters_config[1,]
    
    create_heatmap <-function(s){
        merged2 = merged %>% 
            # filter(n_ref_count>=20, n_vaf<=0.01, t_alt_count>=7) %>%
            filter(n_ref_count>=20, n_vaf<=0.01, t_vaf>=0.01, t_alt_count>=2) %>%
            # filter(n_ref_count>=20, n_vaf<=0.01, ccf>=0.03, t_alt_count>=3) %>%
            filter(sample_group == s, cluster_id != "Unclustered") %>%
            filter(!(sample %in% c("G74-C1","G74-CP1","G74-CP1-twist","G74-DP1","G74-DP1-twist","G74-E1","G74-EP1","G74-EP1-twist"))) %>%
            mutate(cluster_order = -1,
                   sample=str_replace(sample,"G74-",""))
        print(dim(merged2))
        
        if (s == "G74-E"){
            merged2$sample = factor(merged2$sample, levels=paste0("E",c(4,5,6,8,9,10,11,12,13)))
            merged2$cluster_order[merged2$cluster_id == 4] = 0
            merged2$cluster_order[merged2$cluster_id == 10] = 1
            merged2$cluster_order[merged2$cluster_id == 2] = 2
            merged2$cluster_order[merged2$cluster_id == 11] = 3
            merged2$cluster_order[merged2$cluster_id == 12] = 4
            merged2$cluster_order[merged2$cluster_id == 0] = 5
            merged2$cluster_order[merged2$cluster_id == 7] = 6
            merged2$cluster_order[merged2$cluster_id == 8] = 7
            merged2$cluster_order[merged2$cluster_id == 9] = 8
            merged2$cluster_order[merged2$cluster_id == 5] = 9
            merged2$cluster_order[merged2$cluster_id == 6] = 10
            merged2$cluster_order[merged2$cluster_id == 1] = 11
        } else if (s == "G74-C"){
            merged2$cluster_order[merged2$cluster_id == 3] = 0
            merged2$cluster_order[merged2$cluster_id == 0] = 1
            merged2$cluster_order[merged2$cluster_id == 1] = 2
            merged2$cluster_order[merged2$cluster_id == 2] = 3
            merged2$cluster_order[merged2$cluster_id == 4] = 4
        } else if (s == "G74-D"){
            merged2$cluster_order[merged2$cluster_id == 0] = 0
            merged2$cluster_order[merged2$cluster_id == 4] = 1
            merged2$cluster_order[merged2$cluster_id == 2] = 2
            merged2$cluster_order[merged2$cluster_id == 3] = 3
            merged2$cluster_order[merged2$cluster_id == 1] = 4
            merged2$cluster_order[merged2$cluster_id == 5] = 5
        }
        
        # merged2 = filter_variants(merged2, filters)
        
        merged2 = merged2 %>% mutate(cluster_id=as.numeric(cluster_id),
               # Hugo_Symbol=fct_reorder(Hugo_Symbol, cluster_id, .desc = T)) %>%
               # pos=fct_reorder(paste0(chr,"-",start), cluster_id, .desc = T))
               pos=fct_reorder(paste0(chr,"-",start), cluster_order, .desc = T))
        # pos=paste0(chr,"-",start))
        
        p = merged2 %>%
            ggplot(aes(x=sample, y=pos, fill=as.factor(cluster_id), alpha=t_vaf)) +
            # ggplot(aes(x=sample, y=pos, fill=as.factor(cluster_id))) +
            # ggplot(aes(x=sample, y=pos, fill=as.factor(cluster_id), alpha=vcf_vaf)) +
            # ggplot(aes(x=sample, y=Hugo_Symbol, fill=as.factor(cluster_id))) +
            # ggplot(aes(x=sample, y=Hugo_Symbol, fill=vcf_vaf)) +
            geom_tile(na.rm = F) +
            # scale_fill_manual(wes_palette("Zissou1", n = 5)) +
            # scale_fill_manual(values = unname(piratepal("basel", trans = 0.2))) +
            # scale_fill_brewer(palette = "Pastel1") +
            # scale_fill_brewer(palette = "Pastel2") +
            # scale_fill_brewer(palette = "Accent") +
            # scale_fill_brewer(palette = "Set1") +
            # scale_fill_brewer(palette = "Set2") +
            # scale_fill_brewer(palette = "Set3") +
            scale_y_discrete(breaks = merged2$pos, labels = merged2$Hugo_Symbol) +
            # scale_y_continuous(name = merged$Hugo_Symbol) +
            labs(x = "Sample", y = "Gene", fill = "PyClone Clusters", alpha = "VAF") +
            ggtitle(paste("Patient",str_replace(s,"G74-",""))) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold", hjust = 0.5))
        p
        return(p)
    }
    # create_heatmap("G74-D")
    C = create_heatmap("G74-C") + theme(axis.text.y = element_text(size = 8),text = element_text(size = 14))
    D = create_heatmap("G74-D") + theme(text = element_text(size = 14))
    E = create_heatmap("G74-E") + theme(axis.text.y = element_blank(), text = element_text(size = 14))
    E
    (C | D)
    
    CD = (C | D)
    ggsave(filename = "output/connor_pipeline/heatmap/G74-C_and_G74-D_heatmap_VAF_intensities.pdf", plot = CD, height = 18, width = 14)
    
    ggsave(filename = "output/connor_pipeline/heatmap/G74-E_heatmap_VAF_intensities.pdf", 
           plot = E, height = 18, width = 14)
    
    CDE = (C | D | E) + plot_layout(widths = c(1,1,2))
    ggsave(filename = "output/connor_pipeline/heatmap/G74-CDE_heatmap_VAF_intensities.pdf", 
           plot = CDE, height = 18, width = 24, scale = 0.85)
    
    CDE_newline = (C | D) / E + plot_layout(heights = c(7.5,4))
    ggsave(filename = "output/connor_pipeline/heatmap/G74-CDE_on-top_heatmap_VAF_intensities.pdf", 
           plot = CDE_newline, height = 30, width = 24, scale=0.7)
}

