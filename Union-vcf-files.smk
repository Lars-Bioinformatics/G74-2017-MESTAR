__title__ = "Pipeline for Somatic JOINT Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "16/12/2019"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/Data/samples-matched.yaml"

# Explicit paths for external input files
resource_path = "/work/sduvarcall/resources/hg38/"
ref = resource_path+"Homo_sapiens_assembly38.fasta"
gnomad = resource_path+"af-only-gnomad.hg38.vcf.gz"
interval_list = resource_path+"MedExome_hg38_capture_targets.interval_list"
target_regions = resource_path+"MedExome_target_regions/target_regions/"
common_variants = resource_path+"small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
pon_1000g = resource_path+"1000g_pon.hg38.vcf.gz"


SAMPLES = [sample for sample in config]
TUMORS = [config[sample]["tumor"] for sample in config]
NORMALS = [config[sample]["normal"] for sample in config]

print(SAMPLES)
print(TUMORS)
print(NORMALS)

#########################################################
####                      Output                     ####
#########################################################
# output_somatic = "mutect2_force_calling_mergedMatchedAndJointCalling/"
output_somatic = "mutect2_force_calling_mergedMatched/"
# output_somatic = "mutect2_force_calling_mergedMatchedAndJointCalling/"


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

onstart:
    shell("mkdir -p "+output_somatic+"union_vcf/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# print(expand(output_somatic+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=SAMPLES))

rule all_pairs:
    input:
        expand(output_somatic+"union_vcf/{sample}_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)
        # expand(output_somatic+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz", sample=SAMPLES)


#1st pipe, splits multi-allelic calls into separate variant calls
#2nd pipe, left-aligns indels and issues warnings when the REF base in your VCF does not match the base in the supplied FASTA reference genome
# rule vcf_norm:
#     input:
#         vcf=output_somatic+"first_variant_calling_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
#     output:
#         norm_vcf=output_somatic+"union_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
#     shell:
#         """
#         bcftools norm -m-any {input.vcf} | \
#             bcftools norm -Ou --check-ref w -f {ref} | \
#             bcftools filter -Ou -i 'FILTER="PASS"' | \
#             bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
#         bcftools index -t {output.norm_vcf}
#         """
#         # vcf index: tabix -p vcf {output.norm_vcf}
#         # Alternative to vt
#         # bcftools norm -m-any {input.vcf} > {output.vcf}

rule vcf_joint_calling_norm:
    input:
        vcf=output_somatic+"first_variant_calling_vcf/{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
    output:
        norm_vcf=output_somatic+"union_vcf/norm_vcf/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
    shell:
        """
        bcftools norm -m-any {input.vcf} | \
            bcftools norm -Ou --check-ref w -f {ref} | \
            bcftools filter -Ou -i 'FILTER="PASS"' | \
            bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' > {output.norm_vcf}
        bcftools index -t {output.norm_vcf}
        """

rule vcf_merge:
    input:
        norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/norm_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"])
    output:
        merged_vcf=output_somatic+"union_vcf/{sample}_mutect2_matched_filterFlag_norm_PASSonly_merged.vcf.gz"
    shell:
        """
        bcftools merge -Oz -m none --force-samples {input.norm_vcfs} > {output.merged_vcf}
        bcftools index -t {output.merged_vcf}
        """


rule vcf_merge_matchedAndJoint:
    input:
        norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/norm_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"]),
        norm_vcf_joint=output_somatic+"union_vcf/norm_vcf/{sample}_somatic_mutect2_filterFlag_norm_PASSonly.vcf.gz"
    output:
        merged_vcf=output_somatic+"union_vcf/{sample}_mutect2_matchedAndJoint_filterFlag_norm_PASSonly_merged.vcf.gz"
    shell:
        """
        bcftools merge -Oz -m none --force-samples {input.norm_vcfs} {input.norm_vcf_joint} > {output.merged_vcf}
        bcftools index -t {output.merged_vcf}
        """

rule vcf_union:
    input:
        norm_vcfs=lambda wildcards: expand(output_somatic+"union_vcf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf.gz", tumor=config[wildcards.sample]["tumor"], normal=config[wildcards.sample]["normal"])
    output:
        merged_vcf=output_somatic+"sites/sites.txt"
    shell:
        """
        bcftools isec -Oz -n+1 -p sites {input.norm_vcfs}
        """

# bcftools norm -m-any G74-N*gz | bcftools norm -Ou --check-ref w -f /work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta | bcftools filter -Ou -i 'FILTER="PASS"' | grep -v ^# | wc -l