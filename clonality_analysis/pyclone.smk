from snakemake.utils import R

WORK = "/work/Data/"
# SNV_INPUT = WORK+"Connor/mutect2_joint_calling_try3/vcf_filterFlag_PASS/maf_files/"
# SNV_INPUT = WORK+"Connor/mutect2_ctDNA_pipeline/mutect2_matched_force-calling/vcf_final/maf_files/"
SNV_INPUT = WORK+"Connor/mutect2_ctDNA_pipeline_noFFPE/mutect2_matched_force-calling/vcf_final/maf_files/"
CNV_INPUT = WORK+"Connor/sequenza_grouped_files/sequenza_segments/"
PLOIDY_INPUT = WORK+"Connor/sequenza_grouped_files/sequenza_ploidy/"
OUTPUT = WORK+"Connor/pyclone_noPlasma_noFFPE_alsoNotInJointCalling/"
# OUTPUT = WORK+"Connor/pyclone_noPlasma_noFFPE/"
# OUTPUT = WORK+"Connor/pyclone/"

# configfile: WORK+"samples-matched_noTwist.yaml"
configfile: WORK+"samples-pyclone_matched_noTwist.yaml"

onstart:
    shell("mkdir -p "+OUTPUT)
    shell("mkdir -p "+OUTPUT+"input_tsv")


ALL_PATIENTS = [patient for patient in config]
PATIENTS = [patient for patient in config]
PATIENTS.remove("G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY")
print(ALL_PATIENTS)
print(PATIENTS)


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
        # expand(OUTPUT+"{patient}_pyclone/input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
        #     patient="G74-N_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY",
        #     tumor="G74-N4_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY",
        #     normal="G74-NN1_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY")
        # expand(OUTPUT+"{patient}_pyclone_manual/tables/cluster.tsv", patient="G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"),
        # expand(OUTPUT+"{patient}_pyclone_manual/tables/loci.tsv", patient="G74-E_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"),
        expand(OUTPUT+"{patient}_pyclone_manual/tables/cluster.tsv", patient=ALL_PATIENTS),
        expand(OUTPUT+"{patient}_pyclone_manual/tables/loci.tsv", patient=ALL_PATIENTS),
        # expand(OUTPUT+"{patient}_pyclone", patient=PATIENTS),
        # expand(OUTPUT+"{patient}_pyclone_connected", patient=PATIENTS),


# checkpoint sequenza_ploidy:
#     output: 
#         PLOIDY_INPUT+"Sequenza_{tumor}_vs_{normal}_confints_CP.txt"


# rule prepare_pyclone_data:
#     input:
#         # snv_data = SNV_INPUT+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.maf",
#         snv_data = SNV_INPUT+"{tumor}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
#         cnv_data = CNV_INPUT+"Sequenza_{tumor}_vs_{normal}_segments.txt"
#     output:
#         pyclone_tsv = OUTPUT+"input_tsv/{tumor}_vs_{normal}_input_pyclone.tsv"
#     script:
#         WORK+"scripts/prepare_pyclone_data.R"

rule prepare_pyclone_data:
    input:
        snv_data = lambda wildcards: expand(SNV_INPUT+"{tumor}_vs_{normal}__somatic_mutect2-forced_filterFlag_norm_final.maf",
                                               tumor=config[wildcards.patient]["tumor"],
                                               normal=config[wildcards.patient]["normal"]),
        cnv_data = lambda wildcards: expand(CNV_INPUT+"Sequenza_{tumor}_vs_{normal}_segments.txt",
                                            tumor=config[wildcards.patient]["tumor"],
                                            normal=config[wildcards.patient]["normal"]),
    output:
        patient_tsv = OUTPUT+"input_tsv/{patient}_input_pyclone.tsv",
        debug_image = OUTPUT+"r_debug_images/{patient}_r_debug_image.RData"
    script:
        WORK+"scripts/prepare_pyclone_data.R"

rule extract_sample:
    input:
        patient_tsv = OUTPUT+"input_tsv/{patient}_input_pyclone.tsv"
    output:
        pyclone_tsv = OUTPUT+"input_tsv/{patient}_tsv/{tumor}_vs_{normal}_input_pyclone.tsv"
    params:
        tumor = lambda wildcards: wildcards.tumor.split("_")[0]
    shell:
        """
        cat {input} | egrep '^mut|{params.tumor}' > {output}
        """
        # """
        # with open(input[0], 'r') as f:
        #     #header = read.line()
        #     with open(output[0],'w') as out:
        #         for line in f:
        #             line = line.split("\t").trim()
        #             print(line)
        #             if line[7] == "wildcards.tumor":
        #                 out.write('\t'.join(line[:7]))
        # """
        # R("""
        # data = read.table({input}, header=T)
        # data = data[data$sample=={wildcards.tumor},]
        # write.table(data, {output}, quotes = F, row.names = F, sep = "\t")
        # """)
    # shell:
    #     """
    #     awk -F'\t' '{{if ($7=={wildcards.tumor}) printf("$1\t$2\t$3\t$4\t$5\t$6")}}' {input} > {output}
    #     """

rule run_pyclone:
    input:
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{{patient}}_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
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
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{{patient}}_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
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
        pyclone_tsv = lambda wildcards: expand(OUTPUT+"input_tsv/{{patient}}_tsv/{tumor}_vs_{normal}_input_pyclone.tsv", 
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