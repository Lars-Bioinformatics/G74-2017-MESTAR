WORK = "/work/Data/"
SNV_INPUT = WORK+"Connor/mutect2_joint_calling_try3/vcf_filterFlag_PASS/maf_files/"
CNV_INPUT = WORK+"Connor/sequenza_segments/"
PLOIDY_INPUT = WORK+"Connor/sequenza_ploidy/"
OUTPUT = WORK+"Connor/pyclone/"

configfile: WORK+"samples-matched.yaml"

onstart:
    shell("mkdir -p "+OUTPUT)
    shell("mkdir -p "+OUTPUT+"input_tsv")


PATIENTS = [patient for patient in config]
PATIENTS.remove("G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY")
print(PATIENTS)


# Install Sciclone
#install IRanges from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IRanges")
#install devtools if you don't have it already
install.packages("devtools")
library(devtools)
install_github("genome/bmm")
install_github("genome/sciClone")








def getTumorContent(wildcards):
    tumor_content = []
    normal=config[wildcards.patient]["normal"]
    for tumor in config[wildcards.patient]["tumor"]:
        with open(PLOIDY_INPUT+"Sequenza_"+tumor+"_vs_"+normal+"_confints_CP.txt", "r") as f:
            # print(f)
            for i, line in enumerate(f):
                if i == 2: # third line
                    tumor_content.append(line.split("\t")[0])
    # print(checkpoints.sequenza_ploidy.get(**wildcards).output)
    # print(checkpoints.sequenza_ploidy.get(**wildcards).output[1])
    # with open(checkpoints.sequenza_ploidy.get(**wildcards).output[1], 'r') as f:
    return(' '.join(tumor_content))


rule all:
    input:
        expand(OUTPUT+"{patient}_pyclone_manual/tables/cluster.tsv", patient="G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"),
        expand(OUTPUT+"{patient}_pyclone_manual/tables/loci.tsv", patient="G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"),
        expand(OUTPUT+"{patient}_pyclone", patient=PATIENTS),
        expand(OUTPUT+"{patient}_pyclone_connected", patient=PATIENTS),


# checkpoint sequenza_ploidy:
#     output: 
#         PLOIDY_INPUT+"Sequenza_{tumor}_vs_{normal}_confints_CP.txt"


rule prepare_pyclone_data:
    input:
        snv_data = SNV_INPUT+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.maf",
        cnv_data = CNV_INPUT+"Sequenza_{tumor}_vs_{normal}_segments.txt"
    output:
        pyclone_tsv = OUTPUT+"input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv"
    script:
        WORK+"scripts/prepare_pyclone_data.R"


rule run_pyclone:
    input:
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
        sequenza_ploidy = lambda wildcards: expand(PLOIDY_INPUT+"Sequenza_{tumor}_vs_{normal}_confints_CP.txt", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
    output:
        out_dir = directory(OUTPUT+"{patient}_pyclone/")
    params:
        samples = lambda wildcards: [tumor.split("_")[0] for tumor in config[wildcards.patient]["tumor"]],
        tumor_contents = getTumorContent
    shell:
        """
        mkdir -p ~/.config/matplotlib
        echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

        conda run --name pyclone \
        PyClone run_analysis_pipeline \
            --in_files {input.pyclone_tsv} \
            --working_dir {output.out_dir} \
            --tumour_contents {params.tumor_contents} \
            --samples {params.samples} \
            --plot_file_format pdf
        """
# eval "$(conda shell.bash hook)"
# conda activate pyclone

rule run_pyclone_connected:
    input:
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
        sequenza_ploidy = lambda wildcards: expand(PLOIDY_INPUT+"Sequenza_{tumor}_vs_{normal}_confints_CP.txt", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
    output:
        out_dir = directory(OUTPUT+"{patient}_pyclone_connected/")
    params:
        samples = lambda wildcards: [tumor.split("_")[0] for tumor in config[wildcards.patient]["tumor"]],
        tumor_contents = getTumorContent
    shell:
        """
        mkdir -p ~/.config/matplotlib
        echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

        conda run --name pyclone \
        PyClone run_analysis_pipeline \
            --in_files {input.pyclone_tsv} \
            --working_dir {output.out_dir} \
            --tumour_contents {params.tumor_contents} \
            --samples {params.samples} \
            --plot_file_format pdf \
            --init_method connected
        """


# Semi manual run of pyclone for large mutations/many samples
# Causes the plotting code to break

rule run_pyclone_setup:
    input:
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
        sequenza_ploidy = lambda wildcards: expand(PLOIDY_INPUT+"Sequenza_{tumor}_vs_{normal}_confints_CP.txt", 
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
    output:
        out_dir = directory(OUTPUT+"{patient}_pyclone_manual/"),
        config = OUTPUT+"{patient}_pyclone_manual/config.yaml"
    params:
        samples = lambda wildcards: [tumor.split("_")[0] for tumor in config[wildcards.patient]["tumor"]],
        tumor_contents = getTumorContent
    shell:
        """
        conda run --name pyclone \
        PyClone setup_analysis \
            --in_files {input.pyclone_tsv} \
            --working_dir {output.out_dir} \
            --tumour_contents {params.tumor_contents} \
            --samples {params.samples}
        """

rule run_pyclone_run:
    input:
        config = OUTPUT+"{patient}_pyclone_manual/config.yaml"
    output:
        out_dir = directory(OUTPUT+"{patient}_pyclone_manual/trace"),
    shell:
        """
        conda run --name pyclone \
        PyClone run_analysis \
            --config_file {input.config} \
            --seed 42
        """

rule run_pyclone_build_table_cluster:
    input:
        config = OUTPUT+"{patient}_pyclone_manual/config.yaml",
        out_dir = OUTPUT+"{patient}_pyclone_manual/trace/",
    output:
        cluster_tsv = OUTPUT+"{patient}_pyclone_manual/tables/cluster.tsv",
    shell:
        """
        conda run --name pyclone \
        PyClone build_table \
            --config_file {input.config} \
            --out_file {output.cluster_tsv} \
            --table_type cluster
        """

rule run_pyclone_build_table_loci:
    input:
        config = OUTPUT+"{patient}_pyclone_manual/config.yaml",
        out_dir = OUTPUT+"{patient}_pyclone_manual/trace/",
    output:
        loci_tsv = OUTPUT+"{patient}_pyclone_manual/tables/loci.tsv",
    shell:
        """
        conda run --name pyclone \
        PyClone build_table \
            --config_file {input.config} \
            --out_file {output.loci_tsv} \
            --table_type loci
        """