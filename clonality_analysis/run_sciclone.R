# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("IRanges")
# BiocManager::install("limma")
# #install devtools if you don't have it already
# install.packages("devtools")
# library(devtools)
# install_github("genome/bmm")
# install_github("genome/sciClone")

library(sciClone)
library(tidyverse)

setwd("/Users/lars/OneDrive - Syddansk Universitet/Projekter/G74-2017-MESTAR/output/connor_pipeline/sciclone/")

patients = paste0("G74-",c("C","D","E","F","G","K","M","N"))
patients = paste0("G74-",c("F","G","K","M","N"))
patients = paste0("G74-",c("M","N"))

for (patient in patients){
    # patient = patients[1]
    print(patient)
    
    snv_files = list.files("input_tsv", pattern = paste0(patient,".+snv"), full.names = T)
    snv_data = lapply(snv_files, read.table, header=T)
    cnv_files = list.files("input_tsv", pattern = paste0(patient,".+cnv"), full.names = T)
    cnv_data = lapply(cnv_files, read.table, header=T)
    regions_exclude_files = list.files("input_tsv", pattern = paste0(patient,".+regions-to-exclude"), full.names = T)
    regions_exclude_data = do.call(rbind,lapply(regions_exclude_files, read.table, header=T)) %>% arrange(chr,start)
    sample_names = word(snv_files, 3, sep = "_|/")
    
    #3d clustering using three samples:
    sc = sciClone(vafs=snv_data,
                  copyNumberCalls=cnv_data,
                  sampleNames=sample_names, 
                  minimumDepth = 8,
                  cnCallsAreLog2 = T,
                  regionsToExclude=regions_exclude_data)
    
    #create output
    dir.create(paste0(patient,"_results"), showWarnings = F)
    writeClusterTable(sc, paste0(patient,"_results/clusters2"))
    sc.plot1d(sc, paste0(patient,"_results/clusters2.1d.pdf"))
    sc.plot2d(sc, paste0(patient,"_results/clusters2.2d.pdf"))
    # sc.plot3d(sc, sc@sampleNames, size=700, outputFile=paste0(patient,"_results/clusters3.3d.gif"))
}
