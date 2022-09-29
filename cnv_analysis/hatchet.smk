from snakemake.utils import R

WORK = "/work/Data/"
INPUT_BAM = WORK+"Connor/bam/"
OUTPUT = WORK+"Connor/hatchet/"

configfile: WORK+"samples-matched_noTwist.yaml"

ref_build = "hg38"

# Resources - paths inside docker
if ref_build == "hg19":
    # GRCh37 paths
    ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
    chr_prefix = "False"
else:
    # GRCh38 paths
    ref = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta"
    chr_prefix = "True"


# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G56-sampleA1_nimblegen-medexome_HGL2LDSXX"]
# SAMPLES = "G74-C_thruplex-tag-seq-hv-nimblegen-medexome_HVFFJDSXY"
print(SAMPLES)

processes = 21


rule all:
    input:
        # expand(OUTPUT+"{sample}_hatchet_cnv/{sample}_hatchet.ini", sample=SAMPLES)
        expand(OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/results", sample=SAMPLES, maxcov=[300,500])
        # expand(OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/results", sample=SAMPLES, maxcov=[300])



###############################################################################
#### Sequenza
###############################################################################
rule create_hatchet_config:
    input:
        normal=lambda wildcards: expand(INPUT_BAM+"{normal}.connor.recalibrated.bam", normal=config[wildcards.sample]["normal"]),
        tumors=lambda wildcards: expand(INPUT_BAM+"{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
    output:
        hatchet_config=OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/{sample}_hatchet.ini"
    params:
        tumor_names=lambda wildcards: [tumor.split("_")[0] for tumor in config[wildcards.sample]["tumor"]],
        # tumor_names="hej",
        output=OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/"
    shell:
        """
        cat << EOT >> {output}
[run]
count_reads=True
genotype_snps=True
count_alleles=True
combine_counts=True
cluster_bins=True
plot_bins=True
compute_cn=True
plot_cn=True
reference={ref}
output={params.output}
processes={processes}
normal={input.normal}
bams={input.tumors}
samples={params.tumor_names}

[count_reads]
size=50kb

[genotype_snps]
mincov=8
maxcov={wildcards.maxcov}
reference_version={ref_build}
chr_notation={chr_prefix}

[count_alleles]
mincov=8
maxcov={wildcards.maxcov}
EOT
        """


rule run_hatchet:
    input:
        normal=lambda wildcards: expand(INPUT_BAM+"{normal}.connor.recalibrated.bam", normal=config[wildcards.sample]["normal"]),
        tumors=lambda wildcards: expand(INPUT_BAM+"{tumor}.connor.recalibrated.bam", tumor=config[wildcards.sample]["tumor"]),
        hatchet_config=OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/{sample}_hatchet.ini"
    output:
        hatchet_cnv=directory(OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}/results")
    params:
        workdir=OUTPUT+"{sample}_hatchet_cnv_maxcov{maxcov}"
    resources: cpus=21, mem_mb=120000
    shell:
        """
        echo "Running hachet CNV"
        export GRB_LICENSE_FILE="/work/Data/Connor/hatchet/gurobi.lic"
        conda run --name hatchet \
        python -m hatchet run {input.hatchet_config}
        echo "Finished"
        """
