library(vcfR)
library(tidyverse)
library(data.table)
library(patchwork)

setwd("~/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/")

output = "output/connor_pipeline/compare_VAF_calculations/"

read_mafs_dt <- function(path, skip_E=F){
    if (skip_E){ # Varscan
        maf.files = setdiff(list.files(path, full.names = T), list.files(path, pattern = "G74-E", full.names = T))
    } else {
        maf.files = list.files(path, pattern = ".maf", full.names = T)
    }
    data = rbindlist(lapply(maf.files, fread, sep="\t"))
 
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
            setNames(data.frame(chr_pos[,1],as.numeric(chr_pos[,2]), unname(mutect2_vaf[,vcf_sample]), sample), c("chr","pos","vcf_vaf","sample"))
        }
    )
}

filter_variants <- function(merged, filters, tumorOnly=F){
    if (tumorOnly){
        print(filters)
        t_minVAF=as.numeric(filters[1])
        t_alt_minReads=as.numeric(filters[2])
        
        # Filter
        filtered_variants = merged %>%
            filter(material=="tumor") %>%
            group_by(chr, start, end, sample_group) %>% 
            summarise(t_alt_count=max(t_alt_count), t_vaf=max(t_vaf)) %>%
            filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads)
        # filter(t_vaf >= t_minVAF, t_alt_count >= t_alt_minReads)
        
        data = semi_join(merged, filtered_variants, by = c("chr", "start", "end", "sample_group"))
        return(data)
        
    } else {
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
}

snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_ctDNA_pipeline_maf/"
maf_matched = read_mafs_dt(path = snv_path)
dim(maf_matched)
sum(is.na(maf_matched$t_vaf))
# maf_matched$t_vaf[is.na(maf_matched$t_vaf)] = 0

maf %>% filter(Tumor_Sample_Barcode == "G74-N4") %>% dplyr::select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, t_vaf) %>% head
maf %>% dplyr::select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, t_vaf) %>% filter(Start_Position==1242596)

filters=c(t_minVAF=0, t_alt_minReads=7, n_maxVAF = 0.01, n_ref_minReads=20)
maf_matched_filtered = maf_matched %>% 
    rename(chr=Chromosome,start=Start_Position,end=End_Position,sample=Tumor_Sample_Barcode) %>% 
    mutate(sample_group=str_sub(sample, 1, 5)) %>%
    filter_variants(filters)
dim(maf_matched_filtered)

snv_path = "data/connor_pipeline/SNV/maf_files/mutect2_ctDNA_pipeline_maf-tumorOnly/"
maf_tumorOnly = read_mafs_dt(path = snv_path)
dim(maf_tumorOnly)
sum(is.na(maf_tumorOnly$t_vaf))
# maf_tumorOnly$t_vaf[is.na(maf_tumorOnly$t_vaf)] = 0

filters=c(t_minVAF=0, t_alt_minReads=7, n_maxVAF = 0.01, n_ref_minReads=20)
maf_tumorOnly_filtered = maf_tumorOnly %>% 
    rename(chr=Chromosome,start=Start_Position,end=End_Position,sample=Tumor_Sample_Barcode) %>% 
    mutate(sample_group=str_sub(sample, 1, 5)) %>%
    filter_variants(filters, tumorOnly = T)
dim(maf_tumorOnly_filtered)

vcf_path = "data/connor_pipeline/SNV/vcf_files/mutect2_ctDNA_pipeline/vcf_matched/"
vcf_matched = read_vaf_from_vaf(vcf_path)
sum(is.na(vcf_matched$vaf))
head(vcf_matched)

vcf_path = "data/connor_pipeline/SNV/vcf_files/mutect2_ctDNA_pipeline/vcf_tumorOnly/"
vcf_tumorOnly = read_vaf_from_vaf(vcf_path)
sum(is.na(vcf_tumorOnly$vaf))
head(vcf_tumorOnly)

dim(maf_matched)
dim(maf_tumorOnly)
merged = inner_join(dplyr::select(maf_tumorOnly_filtered,chr,pos=start,sample,vaf_reads_matched=t_vaf), dplyr::select(maf_matched,chr=Chromosome,pos=Start_Position,vaf_reads_tumorOnly=t_vaf,sample=Tumor_Sample_Barcode))
merged = inner_join(merged, dplyr::rename(as_tibble(vcf_matched),vaf_mutect2_matched=vcf_vaf))
merged = inner_join(merged, dplyr::rename(as_tibble(vcf_tumorOnly),vaf_mutect2_tumorOnly=vcf_vaf))
head(merged)
dim(merged)

mean(na.omit(maf_matched$t_vaf))
mean(na.omit(maf_tumorOnly$t_vaf))
mean(na.omit(vcf_matched$vcf_vaf))
mean(na.omit(vcf_tumorOnly$vcf_vaf))

sd(na.omit(maf_matched$t_vaf))
sd(na.omit(maf_tumorOnly$t_vaf))
sd(na.omit(vcf_matched$vcf_vaf))
sd(na.omit(vcf_tumorOnly$vcf_vaf))

p = ggplot() + 
    geom_qq(aes(sample=maf_matched$t_vaf), color="blue") +
    geom_qq(aes(sample=maf_tumorOnly$t_vaf), color="green") +
    geom_qq(aes(sample=vcf_matched$vcf_vaf), color="red") +
    geom_qq(aes(sample=vcf_tumorOnly$vcf_vaf), color="orange")
p
ggsave(paste0(output,"VAF_qq.png"), plot = p, width = 12, height = 8)

p1 = ggplot(merged) + 
    geom_point(aes(x=vaf_reads_matched, y=vaf_mutect2_matched), color="blue")
p2 = ggplot(merged) + 
    geom_point(aes(x=vaf_reads_matched, y=vaf_reads_tumorOnly), color="green")
p5 = ggplot(merged) + 
    geom_point(aes(x=vaf_reads_matched, y=vaf_mutect2_tumorOnly), color="purple")
p3 = ggplot(merged) + 
    geom_point(aes(x=vaf_reads_tumorOnly, y=vaf_mutect2_tumorOnly), color="red")
p4 = ggplot(merged) + 
    geom_point(aes(x=vaf_mutect2_matched, y=vaf_mutect2_tumorOnly), color="orange")
p6 = ggplot(merged) + 
    geom_point(aes(x=vaf_reads_tumorOnly, y=vaf_mutect2_matched), color="brown")

p = (p1 + p2 + p5) /
    (p3 + p4 + p6)
p
ggsave(paste0(output,"VAF-VAF_comparison.png"), plot = p, width = 16, height = 9)


merged_sorted = merged %>% 
    arrange(vaf_reads_matched) %>%
    mutate(id1=row_number()) %>%
    arrange(vaf_reads_tumorOnly) %>%
    mutate(id2=row_number()) %>%
    arrange(vaf_mutect2_matched) %>%
    mutate(id3=row_number()) %>%
    arrange(vaf_mutect2_tumorOnly) %>%
    mutate(id4=row_number())

p = ggplot(merged_sorted) +
    geom_point(aes(x=id1, y=vaf_reads_matched), color="blue") +
    geom_point(aes(x=id2, y=vaf_reads_tumorOnly), color="green") +
    geom_point(aes(x=id3, y=vaf_mutect2_matched), color="red") +
    geom_point(aes(x=id4, y=vaf_mutect2_tumorOnly), color="orange")
p
ggsave(paste0(output,"VAF_increasing_order.png"), plot = p, width = 12, height = 8)

p = vaf_vcf_matched %>%
    ggplot() +
    # geom_histogram(aes(x=vaf))
    geom_jitter(aes(x=chr, y=vaf))
print(p)

