__title__ = "Pipeline for Somatic Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "10/09/2019"
__version__ = "1.0"

import time

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/Data/samples-pairwise.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
resource_path = "/work/sduvarcall/resources/hg38/"
ref = resource_path+"Homo_sapiens_assembly38.fasta"
gnomad = resource_path+"af-only-gnomad.hg38.vcf.gz"
interval_list = resource_path+"MedExome_hg38_capture_targets.interval_list"
target_regions = resource_path+"MedExome_target_regions/target_regions/"
common_variants = resource_path+"small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
# pon_location = "somatic_panel_of_normals/"
pon_1000g = resource_path+"1000g_pon.hg38.vcf.gz"

# Sample information
# PAIR, = glob_wildcards("{sample}.connor.recalibrated.bam")
# PAIR = ["ECV2-4-plasma180205"]
PAIR = [pair for pair in config]
# print(PAIR)


NORMALS = list(set([config[sample]["normal"] for sample in config]))
print(NORMALS)


INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)
print(INTERVALS)


# SAMPLES = "ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"

#########################################################
####                      Output                     ####
#########################################################
log_file = "log_file_somatic.txt"
output_somatic = "mutect2_matched-default_settings/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
# mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    shell("mkdir -p "+output_somatic+"split/{{vcf,bam,f1r2}}_split")
    shell("mkdir -p "+output_somatic+"vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm}}")
    shell("mkdir -p "+output_somatic+"realigned_bam/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# print([expand(output_somatic+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
#             normal=config[pair]["normal"],
#             tumor=config[pair]["tumor"]) for pair in PAIR])


rule all_pairs:
    input:
        [expand(output_somatic+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR]

# #########################################################
# ####       Create Somatic Panel of Normals           ####
# #########################################################
# '''
# Run Mutect2 in tumor-only mode for each normal sample
# '''
# rule Mutect2_tumor_only_pon:
#     input:
#         bam="bam/{sample}.connor.recalibrated.bam",
#         intervals=target_regions+"{interval}.interval_list"
#     output:
#         vcf=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.gz",
#         vcf_stats=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.gz.stats",
#     resources: cpus=3, mem=18000
#     shell:
#         """
#         gatk --java-options {mem} Mutect2 \
#         -R {ref} \
#         -I {input.bam} \
#         -max-mnp-distance 0 \
#         --native-pair-hmm-threads {threads} \
#         -L {input.intervals} \
#         -O {output.vcf}
#         """

# rule merge_normal_vcf:
#     input:
#         vcf_subfile=expand(pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz", interval=INTERVALS)
#     output:
#         vcf=pon_location+"{sample}_for_pon.vcf.gz",
#         tbi=pon_location+"{sample}_for_pon.vcf.gz.tbi"
#     params:
#         vcf_subfile=expand("-I "+pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf.gz", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}

#         tabix -p vcf {output.vcf}
#         """

# rule GenomicsDB:
#     input:
#         vcf=expand(pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS)
#     output:
#         db=directory(pon_location+"pon_db")
#     params:
#         vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS)
#     shell:
#         """
#         gatk --java-options {mem} GenomicsDBImport \
#             -R {ref} \
#             -L {interval_list} \
#             --genomicsdb-workspace-path {output.db} \
#             --max-num-intervals-to-import-in-parallel {threads} \
#             --merge-input-intervals true \
#             {params.vcf}
#         """

# '''
# Combine the normal calls using CreateSomaticPanelOfNormals.
# '''
# rule CreateSomaticPanelOfNormals:
#     input:
#         db=pon_location+"pon_db"
#         # vcf=expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
#     output:
#         pon=pon_location+"pon.vcf.gz"
#     # params:
#     #     vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
#     shell:
#         """
#         gatk --java-options {mem} CreateSomaticPanelOfNormals \
#         -R {ref} \
#         -V gendb://{input.db} \
#         -O {output}
#         """
#         # {params.vcf} \


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
    input:
        normal="bam/{normal}.connor.recalibrated.bam",
        tumor="bam/{tumor}.connor.recalibrated.bam",
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz",
        idx=output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.tbi",
        vcf_stats=output_somatic+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{interval}__split.vcf.gz.stats",
        bam_subfile=output_somatic+"split/bam_split/{tumor}_vs_{normal}__{interval}_mutect2.bam",
        f1r2=output_somatic+"split/f1r2_split/{tumor}_vs_{normal}_f1r2__{interval}__split.tar.gz"
    resources: cpus=3, mem_mb=18000
    shell:
        """
        gatk --java-options -Xmx{resources.mem_mb}g Mutect2 \
        -R {ref} \
        -I {input.tumor} \
        -I {input.normal} \
        -normal {wildcards.normal} \
        -pon {pon_1000g} \
        --germline-resource {gnomad} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -bamout {output.bam_subfile} \
        -L {input.intervals} \
        -O {output.vcf}
        """

rule merge_somatic_vcf:
    input:
        vcf_subfile=expand(output_somatic+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    output:
        vcf=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
    params:
        vcf_subfile=expand("-I "+output_somatic+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}

        tabix -p vcf {output.vcf}
        """

rule merge_mutect2_bam:
  input:
        bam_subfile=expand(output_somatic+"split/bam_split/{{tumor}}_vs_{{normal}}__{interval}_mutect2.bam", interval=INTERVALS)
  output:
        bam=output_somatic+"realigned_bam/{tumor}_vs_{normal}_somatic_mutect2.bam"
  params:
        bam_subfile=expand("-I " + output_somatic + "split/bam_split/{{tumor}}_vs_{{normal}}__{interval}_mutect2.bam", interval=INTERVALS)
  shell:
        """
        gatk --java-options {mem} MergeSamFiles \
        {params.bam_subfile} \
        --VALIDATION_STRINGENCY LENIENT \
        --USE_THREADING true \
        -O {output.bam} \
        --CREATE_INDEX true
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand(output_somatic+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    output:
        vcf_stats=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
    params:
        vcf_subfile=expand("-stats "+output_somatic+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{interval}__split.vcf.gz.stats", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

#########################################################
####          Learn Read Orientation Bias            ####
####                                                 ####
#### Note: Used to fix orientation bias artifacts    ####
####       from Formalin-Fixed Paraffin-Embedded     ####
####       (FFPE) samples - i.e. not needed for      ####
####       frozen tissue                             ####
#########################################################
rule learnReadOrientationModel:
    input:
        f1r2=expand(output_somatic+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    output:
        f1r2_model=output_somatic+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+output_somatic+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    shell:
        """
        gatk LearnReadOrientationModel \
        {params.f1r2} \
        -O {output}
        """

#########################################################
####           Create Contamination table            ####
#########################################################
rule GetPileupSummaries:
    input:
        bam="bam/{sample}.connor.recalibrated.bam"
    output:
        pileup=output_somatic+"pileups/{sample}_pileup.table"
    shell:
        """
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule CalculateContamination:
    input:
        normal=output_somatic+"pileups/{normal}_pileup.table",
        tumor=output_somatic+"pileups/{tumor}_pileup.table"
    output:
        contamination=output_somatic+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=output_somatic+"contaminations_tables/{tumor}_vs_{normal}_contamination.table"
    output:
        joined_table=output_somatic+"contaminations_tables/{tumor}_vs_{normal}_merged_contamination.table"
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
rule FilterMutectCalls:
    input:
        vcf=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
        idx=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
        bam=output_somatic+"realigned_bam/{tumor}_vs_{normal}_somatic_mutect2.bam",
        stats=output_somatic+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
        contamination=output_somatic+"contaminations_tables/{tumor}_vs_{normal}_merged_contamination.table",
        read_orientation=output_somatic+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    output:
        vcf=output_somatic+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
        tbi=output_somatic+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz.tbi"
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {interval_list} \
        -O {output.vcf}
        """

rule vcf_pass_only:
    input:
        vcf=output_somatic+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
    output:
        vcf_pass=output_somatic+"vcf_filterFlag_PASS/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
    shell:
        """
        bcftools view -i "%FILTER='PASS' | %FILTER='.'" {input.vcf} | bgzip -c > {output.vcf_pass}
        tabix -p vcf {output.vcf_pass}
        """

rule vcf_norm_decompose:
    input:
        vcf=output_somatic+"vcf_filterFlag/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        norm_vcf=output_somatic+"vcf_filterFlag_norm/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf.gz"
    shell:
        """
        zless {input.vcf} \
        | vt decompose -s - \
        | vt normalize -r {ref} - \
        | bgzip -c > {output.norm_vcf}
        """
        # vcf index: tabix -p vcf {output.norm_vcf}
        # Alternative to vt
        # bcftools norm -m-any {input.vcf} > {output.vcf}