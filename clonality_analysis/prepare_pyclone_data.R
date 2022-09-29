save.image(snakemake@output[[2]])

library(maftools)
library(data.table)
library(tidyverse)

# Input
# input_snv = list.files(path = "SNV/maf_files/mutect2_joint_calling_maf/", pattern = "D2", full.names = T)
# input_cnv = list.files(path = "CNV/Sequenza/sequenza_segments/", pattern = "D2", full.names = T)
print(snakemake@input)
input_snv = snakemake@input$snv_data
input_cnv = snakemake@input$cnv_data
print(paste("input snv:", input_snv))
print(paste("input cnv:", input_cnv))

# Helper functions
manually_calculate_ploidy <- function(data){
    # Tumor ploidy
    data2 = data %>% mutate(segment.size=stop-start)
    genome_lengths = data2 %>% group_by(sample) %>% summarise(genome_length = sum(segment.size))
    data2 = merge(data2, genome_lengths, by = "sample")
    # ploidy = length_of_segment * total_cnv_in_that segment / length_of_genome
    ploidy = data2 %>% group_by(sample) %>% summarise(ploidy = sum(segment.size*total_cn/genome_length),
                                                      ploidy_rounded = round(sum(segment.size*total_cn/genome_length))
    )
    data = merge(data, ploidy, by = "sample")
    return(data)
}

read_sequenza <- function(cnv_files, ploidy_path="", manual_ploidy=F){
    # cnv_files <- list.files(pattern = "segments", path = path)
    data <- map_dfr(
        .x = cnv_files,
        .f = function(file){
            print(file)
            sample = str_split(basename(file), "_", simplify = T)[2]
            print(sample)
            data = read.table(file, header = T) %>% 
                dplyr::select(chr=chromosome, start=start.pos, stop=end.pos, total_cn=CNt, major_cn=A, minor_cn=B) %>%
                # dplyr::select(chr=chromosome, start=start.pos, stop=end.pos, tumor.total=CNt, tumor.major=A, tumor.minor=B)
                mutate(normal_cn=2)
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

combine_snv_cnv_data <- function(snv_data, cnv_data){
    
    snv_data = snv_data %>% dplyr::rename(sample=Tumor_Sample_Barcode, chr=Chromosome, start=Start_Position, end=End_Position)
    cnv_data = cnv_data %>% dplyr::rename(start2=start, end2=stop) %>% as.data.table()
    
    # Set keys used for finding overlaps
    setkey(snv_data, sample, chr, start, end)
    setkey(cnv_data, sample, chr, start2, end2)
    # setkey(snv_data, chr, start, end)
    # setkey(cnv_data, chr, start2, end2)
    
    combined = foverlaps(snv_data,cnv_data) %>%
        mutate(sample_group=str_sub(sample, 1, 5))
    
    return(combined)
}

filter_variants <- function(merged, filters){
    print(filters)
    t_minVAF=as.numeric(filters[1])
    t_alt_minReads=as.numeric(filters[2])
    n_maxVAF=as.numeric(filters[3])
    n_ref_minReads=as.numeric(filters[4])

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

# Read data

snv_data = rbindlist(lapply(input_snv, fread, sep="\t")) %>%
    mutate(t_vaf = t_alt_count/t_depth, n_vaf = n_alt_count/n_depth) %>%
    mutate(material = ifelse(str_detect(Tumor_Sample_Barcode,"P"), "plasma", "tumor"))
print(dim(snv_data))

# gender = ifelse("chrY" %in% snv_data$Chromosome, "male", "female")
# cnv_data = rbindlist(lapply(input_cnv, fread, sep="\t")) %>%
#     dplyr::select(chr=chromosome, start=start.pos, stop=end.pos, total_cn=CNt, major_cn=A, minor_cn=B) %>%
#     mutate(normal_cn=ifelse(chr %in% c("chrX", "chrY") & gender=="male", 1, 2))
# print(dim(cnv_data))

cnv_data = read_sequenza(input_cnv)
head(cnv_data)

filters=c(t_minVAF=0, t_alt_minReads=7, n_maxVAF = 0.01, n_ref_minReads=20)
pyclone_ready = combine_snv_cnv_data(snv_data = snv_data, cnv_data = cnv_data) %>%
    filter_variants(filters) %>%
    mutate(mutation_id = paste0(chr,"_",start,"_",end)) %>%
    dplyr::select(mutation_id, ref_counts=t_ref_count, var_counts=t_alt_count, normal_cn, minor_cn, major_cn, sample) %>%
    drop_na() %>%
    filter(major_cn>0)

# sample = unique(snv_data$Tumor_Sample_Barcode)
write.table(pyclone_ready, file = snakemake@output[[1]], quote = F, row.names = F, sep = "\t")