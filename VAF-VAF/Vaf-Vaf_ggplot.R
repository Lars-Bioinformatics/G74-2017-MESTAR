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
snv_caller = "mutect2"
snv_caller = "varscan2"
snv_caller = "varscan2_1000g-pon"
snv_caller = "varscan2_1000g-pon_q30"
snv_caller = "varscan2_1000g-pon_default-settings"
snv_caller = "varscan2_1000g-pon_gnomad"

cnv_caller = "sequenza"


# Choose output folder
# output = "output/connor_pipeline/VAF_VAF_plots2/"
output = "VAF_VAF_plots2/"

run_vaf_vaf(snv_caller, cnv_caller)


# Run several snv_callers and cnv_callers
for (snv_caller in c("mutect2_joint_calling_pass")){
    for (cnv_caller in c("sequenza")){
        # snv_path = get_snv_path(snv_caller)
        # cnv_path = get_cnv_path(cnv_caller)
        run_vaf_vaf(snv_caller, cnv_caller)
    }
}


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
    
    if (snv_caller=="varscan2") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf/"
    if (snv_caller=="varscan2_1000g-pon") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered/"
    if (snv_caller=="varscan2_1000g-pon_q30") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered_q30/"
    if (snv_caller=="varscan2_1000g-pon_default-settings") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g_pon_filtered_default-settings/"
    if (snv_caller=="varscan2_1000g-pon_gnomad") snv_path = "data/connor_pipeline/SNV/maf_files/varscan2_maf_somatic_1000g-pon_gnomad-filtered/"
    
    return(snv_path)
}

get_cnv_path <- function(cnv_caller){
    if (cnv_caller == "sequenza") cnv_path = "data/connor_pipeline/CNV/Sequenza/sequenza_segments/"
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
                if (is.na(ploidy)){
                    data = manually_calculate_ploidy(data)
                } else {
                    ploidy = as.numeric(ploidy)
                    data$ploidy = ploidy
                    data$ploidy_rounded = round(ploidy)
                }
                
            }
            data
        }
    )
    
    return(data)
}

###############################################################################
### Combined functions for snv-cnv data object
###############################################################################
combine_snv_cnv_data <- function(snv_data, cnv_data){
    
    snv_data = snv_data %>% dplyr::rename(sample=Tumor_Sample_Barcode, chr=Chromosome, start=Start_Position, end=End_Position)
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
        dplyr::select(sample, chr, start, end, t_vaf, t_alt_count, t_ref_count, n_ref_count, Hugo_Symbol, tumor.total, tumor.major, tumor.minor, sample_group, HGVSc, HGVSp) %>%
        # mutate(sample_group=str_sub(sample, 1, 5)) %>%
        inner_join(., ., by=c("chr", "start", "end", "sample_group", "Hugo_Symbol", "HGVSc", "HGVSp")) %>%
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

vaf_vaf_combined_sample <- function(data, sg, output_dir){
    output3 = paste0(output_dir, "Combined_per_patient/")
    dir.create(output3, showWarnings = F, recursive = T)
    data2 = data %>% filter(sample_group==sg)
    p = data2 %>%
        ggplot(aes(t_vaf.x, t_vaf.y)) +
        geom_point(size=3, color="blue",alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="VAF (bottom name)", y="VAF (upper name)") +
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
        geom_point_interactive(aes(x = t_vaf.x, y = t_vaf.y, tooltip = tooltip), size=3, color="blue",alpha=0.3) +
        geom_abline(intercept=0, slope = 1, linetype="dotted") +
        # scale_color_manual(values = brewer.pal(3, "Set1")[1:2], breaks = c("Yes","No")) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(x="VAF (bottom name)", y="VAF (upper name)") +
        # theme_minimal() +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
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


# Mutect2
run_vaf_vaf = function(snv_caller, cnv_caller){
        
    snv_path = get_snv_path(snv_caller)
    cnv_path = get_cnv_path(cnv_caller)
    
    # get snv data
    if (str_detect(snv_caller,"varscan")){
        # Exclude G74-E samples as they are huge!
        # files = setdiff(list.files(snv_path, full.names = T), list.files(snv_path, pattern = "G74-E", full.names = T))
        # files = list.files(snv_path, pattern = "G74-D", full.names = T)
        # maf = merge_mafs(mafs = files)
        # maf = fix_maf_sampleNames_and_add_vaf(maf)
        maf = read_mafs_dt(path = snv_path, skip_E = T)
    } else {
        # maf = read_mafs(path = snv_path)
        maf = read_mafs_dt(path = snv_path)
    }
    
    # avoid long names
    snv_caller2 = snv_caller
    if (str_detect(snv_caller,"mutect2_force_calling")){
        snv_caller <- "mutect2_force_calling"
        snv_caller <<- "mutect2_force_calling"
    }
    
    dim(maf)
    # snv_data = filter(maf, FILTER == "PASS", !is.na(t_alt_count)) # PASS is no longer valid in force-calling
    # snv_data = filter(maf, !is.na(t_alt_count))
    snv_data = maf
    dim(snv_data)
    # dim(maf@data)
    # snv_data = filter(maf@data, FILTER == "PASS", !is.na(t_alt_count))
    # table(filter(maf@data, FILTER != "PASS")$FILTER)
    # snv_data = filter(maf@data, FILTER == "PASS")
    # snv_data = maf@data
    # snv_data = filter(maf@data, !is.na(t_alt_count))
    
    
    # get cnv data
    cnv_data = read_sequenza(path = cnv_path, manual_ploidy = T)
    
    # Merge snv and cnv data
    merged = combine_snv_cnv_data(snv_data = snv_data, cnv_data = cnv_data)
    merged <<- merged
    
    # maf@data %>% dplyr::count(Tumor_Sample_Barcode)
    merged %>% dplyr::count(chr,start,end,sample_group,Variant_Type) %>% arrange(desc(n)) %>% filter(n>0) %>% dplyr::count(sample_group, Variant_Type)
    data = merged
    
    # merged = filter(merged, sample_group %in% c("G74-D"))
    
    # filters_config = expand_grid(t_minVAF=c(0.1,0.2,0.3), t_alt_minReads=c(1,3,5,10), n_maxVAF = 0.01, n_ref_minReads=c(10,20))
    # filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(5,7,10,12), n_maxVAF = 0.01, n_ref_minReads=c(20,50))
    # filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(5,7,10,12), n_maxVAF = 0.01, n_ref_minReads=c(20))
    filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(7), n_maxVAF = 0.01, n_ref_minReads=c(20))
    # filters_config = expand_grid(t_minVAF=c(0), t_alt_minReads=c(10), n_maxVAF = 0.01, n_ref_minReads=c(10,20,30,40,50))
    # filters_config = expand_grid(t_minVAF=c(0.1), t_alt_minReads=c(10), n_maxVAF = 0.01, n_ref_minReads=10)
    # filters_config = rbind(filters_config, c(0,0,1,0)) # no filters
    
    print("hej1")
    cl <- makeCluster(detectCores()-1)
    # cl <- makeCluster(2)
    clusterExport(cl, varlist = c("merged", "filter_variants", "join_vaf_vs_vaf", "vaf_vaf_combined_sample", "output", "snv_caller", "snv_caller2",
                                  "vaf_vaf_combined_sample_interactive","vaf_vaf_per_sample","vaf_vaf_per_sample_interactive"))
    print("hej2")
    # apply(filters_config, 1, function(filters){
    parApply(cl, filters_config, 1, function(filters){
        library(tidyverse)
        library(ggiraph)
        # filters = filters_config[1,]
        # print(filters)
        # rowl <<- filters
        # merged=filter(merged,n_ref_count<500, t_ref_count<1000)
        filtered = filter_variants(merged, filters)
        
        output2 = paste0(output, snv_caller2, "/", snv_caller ,extra,"/")
        
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
            write.table(excel_wide, file = paste0(output2,"filtered_variants/", s, "_filtered_variants.txt"), quote = F, row.names = F, sep = "\t")
            # write.table(excel_wide, file = paste0("output/connor_pipeline/Data tables - Excel/manuscript 1 - variant select/", s, "_variant-select.txt"), quote = F, row.names = F, sep = "\t")
            # write_xlsx(excel_wide, path = paste0("output/connor_pipeline/Data tables - Excel/manuscript 1 - variant select/", s, "_variant-select.xlsx"))
        }
        
        
        # filtered %>% dplyr::select()
        # write.table(filtered, file = "output/connor_pipeline/Data tables - Excel/mutect2_snv_cnv_data.txt", sep = "\t", row.names = F, quote = F)
        
        joined_data = join_vaf_vs_vaf(filtered)
        for (s in unique(joined_data$sample_group)){
            data2 = filter(joined_data, sample_group == s) %>% dplyr::select(-tooltip)
            dir.create(path = paste0(output2,"joined_data/"), recursive = T, showWarnings = F)
            write.table(data2, file = paste0(output2,"joined_data/",s,"_joined_data.txt"), quote = F, row.names = F, sep = "\t")
        }
        
        for (sg in unique(joined_data$sample_group)){
            print(sg)
            vaf_vaf_combined_sample(joined_data, sg, output2)
            if (sg != "G74-E"){
                vaf_vaf_combined_sample_interactive(joined_data, sg, output2)
            }
        }
        
        # pairs = joined_data %>% dplyr::select(sample.x, sample.y) %>% distinct() %>% as.data.frame(., stringsAsFactors=F)
        # for (i in 1:nrow(pairs)){
        #     print(paste(pairs[i,1], pairs[i,2]))
        #     # sub_data = joined_data %>% filter(sample.x == as.character(pairs[i,1]), sample.y==as.character(pairs[i,2]))
        #     sub_data = joined_data %>% filter(sample.x == pairs[i,1], sample.y==pairs[i,2])
        #     vaf_vaf_per_sample(sub_data, pairs[i,1], pairs[i,2], output2)
        #     vaf_vaf_per_sample_interactive(sub_data, pairs[i,1], pairs[i,2], output2)
        # }
        
    })
    stopCluster(cl)
}




# PCA analysis
PCA <- function(){
    # PCA
    standardize <- function(z) {
        rowmed <- apply(z, 1, median)
        rowmad <- apply(z, 1, mad)  # median absolute deviation
        rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
        rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
        return(rv)
    }
    
    # dat2<-joined_data %>% filter(sample.x=="G74-D2", sample.y=="G74-D3")
    dat = data %>% 
        filter(sample_group=="G74-D", !is.na(t_vaf)) %>%
        mutate(var=paste0(chr,start)) %>%
        select(var, sample, t_vaf) %>%
        # reshape(idvar = var, timevar = sample, direction = "wide")
        spread(key = sample, value=t_vaf)
    dat2=dat[,c(2,3,4)]
    dat2=dat[,3:12]
    dim(dat2)
    dat2=na.omit(dat2)
    dat2<-scale(dat2)
    
    fviz_eig(dat.pca)
    
    set.seed(123)
    # dat.pca <- prcomp(t(dat2), center = TRUE,scale. = TRUE)
    dat.pca <- prcomp(t(dat2), center = F,scale = F)
    print(dat.pca)
    x<-dat.pca$rotation[,1]
    y<-dat.pca$rotation[,2]
    z<-dat.pca$rotation[,3]
    # z<-dat.pca$x[,3]
    plot(x,y)
    library(rgl)
    
    smp<-info[,2]
    
    x=dat2$`G74-D2`
    y=dat2$`G74-D3`
    z=dat2$`G74-DP1`
    plot3d(x,y,z, size=7, main="BASELINE SAMPLES")
    
    # with(t(dat.pca),plot3d(x,y,z, col=info[,5], size=7, main="BASELINE SAMPLES"))
    # with(t(dat.pca),text3d(x,y,z,smp,pos=3,cex=0.75))
    with(t(dat.pca),plot3d(x,y,z, size=7, main="BASELINE SAMPLES"))
    with(t(dat2),plot3d(x,y,z, size=7, main="BASELINE SAMPLES"))
    with(t(dat.pca),text3d(x,y,smp,pos=3,cex=0.75))
}





# sequenza = readRDS("Sequenza_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_vs_G65-B10-BA34_nimblegen-medexome_HYLKFDSXX_bin50/Sequenza_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_vs_G65-B10-BA34_nimblegen-medexome_HYLKFDSXX_sequenza_cp_table.RData")
# load("Sequenza_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_vs_G65-B10-BA34_nimblegen-medexome_HYLKFDSXX_bin50/Sequenza_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_vs_G65-B10-BA34_nimblegen-medexome_HYLKFDSXX_sequenza_extract.RData")
# `Sequenza_G65-T10-BA34_nimblegen-medexome_HYLKFDSXX_vs_G65-B10-BA34_nimblegen-medexome_HYLKFDSXX_sequenza_extract`$ratio$chr1






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
