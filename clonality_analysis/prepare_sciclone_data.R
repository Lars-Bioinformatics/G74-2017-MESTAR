library(tidyverse)
library(maftools)

# Input 
setwd("/Users/lars/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/")

# input_snvs = list.files(path = "data/connor_pipeline/SNV/maf_files/mutect2_joint_calling_maf/", pattern = "G74-D", full.names = T)
input_snvs = list.files(path = "data/connor_pipeline/SNV/maf_files/mutect2_joint_calling_maf/", full.names = T)
# input_cnvs = list.files(path = "data/connor_pipeline/CNV/Sequenza/sequenza_segments/", pattern = "G74-D", full.names = T)
input_cnvs = list.files(path = "data/connor_pipeline/CNV/Sequenza/sequenza_segments/", full.names = T)

# Read data
for (i in 1:length(input_snvs)){
    input_snv = input_snvs[i]
    input_cnv = input_cnvs[i]
    
    maf = merge_mafs(mafs = input_snv)
    snv_data = filter(maf@data, FILTER == "PASS", !is.na(t_alt_count)) %>%
        mutate(vaf = t_alt_count/t_depth*100,chr=gsub("chr","",Chromosome)) %>%
        dplyr::select(chr, pos=Start_Position, ref_reads=t_ref_count, 
                      var_reads=t_alt_count, vaf)
    
    sequenza = read.table(input_cnv, header = T)
    cnv_data = sequenza %>%
        mutate(logR=log2(2 * depth.ratio) - 1, chr=gsub("chr","",chromosome)) %>%
        dplyr::select(chr, start=start.pos, stop=end.pos, segment_mean=logR) %>%
        filter(chr!="chrY")
    
    regions_to_exclude = sequenza %>% 
        filter(CNt==2, B==0, chromosome!="chrY") %>%
        select(chr=chromosome, start=start.pos, end=end.pos)
    
    # Write new input files for sciclone
    sample = unique(maf@data$Tumor_Sample_Barcode)
    write.table(snv_data, file = paste0("output/connor_pipeline/sciclone/input_tsv/",sample,"_snv_sciclone_input.tsv"), quote = F, row.names = F)
    write.table(cnv_data, file = paste0("output/connor_pipeline/sciclone/input_tsv/",sample,"_cnv_sciclone_input.tsv"), quote = F, row.names = F)
    write.table(regions_to_exclude, file = paste0("output/connor_pipeline/sciclone/input_tsv/",sample,"_regions-to-exclude_sciclone_input.tsv"), quote = F, row.names = F)
}
