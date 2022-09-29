# Run script from folder with fastq files.

from os import getcwd
import time

# SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")
#SAMPLES = "PC1-10-blod_normal_tagseq-medexome"
#SAMPLES = "PC1-14-recidiv_tumor_tagseq-medexome"
FAMNAME = getcwd().rsplit("/",1)[1]
print(FAMNAME)
OUTPUT = "takara_umi/"

resource_path = "/work/sduvarcall/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bed = resource_path + "MedExome_target_regions/MedExome_hg38_capture_targets.bed"
interval_list = resource_path + "MedExome_target_regions/MedExome_hg38_capture_targets.interval_list"
rcUMI = "/work/sduvarcall/resources/preprocessing/TruSeq3-PE-2with_rcUMI.fa"
expectedUMIs = "/work/sduvarcall/resources/preprocessing/expectedUMI.txt"

### Can be run without bed file

mem = 12
timeFormat = time.strftime("%Y_%m_%d-%X")

onstart:
    shell("mkdir -p " + OUTPUT + " logs_slurm")


rule all:
  input:
    # expand("MarkDuplicates/bam/{sample}.recalibrated.bam", sample=SAMPLES),
    # expand("bam/{sample}.connor.recalibrated.bam", sample=SAMPLES),
    OUTPUT+"quality_control/takara_umi_multiqc_report.html"



###############################################################################
####              Takara ThruPLEX Tag-Seq HV with UMI pipeline             ####
###############################################################################

# Step A. Trim adapters and reverse complement UMIs
rule adapter_trimming:
    input:
        f1="fastq/{sample}_R1.fastq.gz",
        f2="fastq/{sample}_R2.fastq.gz",
    output:
        trimmed_f1 = OUTPUT+"fastq_trimmed/{sample}_trimmed_R1.fastq.gz",
        trimmed_f2 = OUTPUT+"fastq_trimmed/{sample}_trimmed_R2.fastq.gz",
        unpaired_f1 = OUTPUT+"fastq_trimmed/{sample}_unpaired_R1.fastq.gz",
        unpaired_f2 = OUTPUT+"fastq_trimmed/{sample}_unpaired_R2.fastq.gz",
    shell:
        """
        trimmomatic -Xmx{mem}g PE \
        {input.f1} {input.f2} \
        {output.trimmed_f1} {output.unpaired_f1} \
        {output.trimmed_f2} {output.unpaired_f2} \
        ILLUMINACLIP:{rcUMI}:1:10:5:9:true MINLEN:20
        """


# Step B. Generate unmapped BAM with UMI in RX tag
rule unmapped_bam:
    input:
        trimmed_f1 = OUTPUT+"fastq_trimmed/{sampleid}_{protocol}_{flowcell}_trimmed_R1.fastq.gz",
        trimmed_f2 = OUTPUT+"fastq_trimmed/{sampleid}_{protocol}_{flowcell}_trimmed_R2.fastq.gz",
    output:
        unmapped_bam = OUTPUT+"unmapped_bam/{sampleid}_{protocol}_{flowcell}_unmapped.bam"
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{protocol}_{flowcell}",
        rgpl = "illumina",
        rgpu = "{flowcell}"
    shell:
        """
        fgbio -Xmx{mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap FastqToBam \
        -i {input.trimmed_f1} {input.trimmed_f2} \
        -r 7M1S+T 7M1S+T \
        --sample {params.rgsm} \
        --library {params.rglb} \
        --read-group-id {params.rgid} \
        --platform {params.rgpl} \
        --platform-unit {params.rgpu} \
        -o {output.unmapped_bam} \
        -s true
        """

rule samToFastq:
    input:
        unmapped_bam = OUTPUT+"unmapped_bam/{sample}_unmapped.bam"
    output:
        unmapped_f1 = temp(OUTPUT+"fastq_minusUMI/{sample}_minusUMI_R1.fastq.gz"),
        unmapped_f2 = temp(OUTPUT+"fastq_minusUMI/{sample}_minusUMI_R2.fastq.gz"),
    shell:
        """
        gatk --java-options -Xmx{mem}G SamToFastq \
        --INPUT {input.unmapped_bam} \
        --FASTQ {output.unmapped_f1} \
        --SECOND_END_FASTQ {output.unmapped_f2}
        """

# Step C. Align the new FASTQ files (with removed UMIs)
rule bwa_mem:
    input:
        unmapped_f1 = OUTPUT+"fastq_minusUMI/{sampleid}_{protocol}_{flowcell}_minusUMI_R1.fastq.gz",
        unmapped_f2 = OUTPUT+"fastq_minusUMI/{sampleid}_{protocol}_{flowcell}_minusUMI_R2.fastq.gz",
    output:
        bam = OUTPUT+"sorted_bam/{sampleid}_{protocol}_{flowcell}.sorted.bam",
        # bai = temp(OUTPUT+"sorted_bam/{sampleid}_{protocol}_{flowcell}.sorted.bai")
    resources: cpus=9, mem=50000
    log:
        "logs/bwa_logs/{sampleid}_{protocol}_{flowcell}.bwa_mem.log"
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{protocol}_{flowcell}",
        rgpl = "illumina",
        rgpu = "{flowcell}"
    shell:
        """
        bwa mem {ref} {input} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        -M \
        -t {resources.cpus} | \
        gatk --java-options -Xmx{mem}G SortSam \
        --INPUT /dev/stdin \
        --OUTPUT {output.bam} \
        --VALIDATION_STRINGENCY LENIENT \
        --SORT_ORDER queryname \
        --TMP_DIR tmp \
        --CREATE_INDEX TRUE
        """

rule mergeBamAlignment:
    input:
        bam = OUTPUT+"sorted_bam/{sample}.sorted.bam",
        unmapped_bam = OUTPUT+"unmapped_bam/{sample}_unmapped.bam",
    output:
        aligned_bam = temp(OUTPUT+"umi_bam/{sample}.sorted.aligned.umi.bam"),
        aligned_bai = temp(OUTPUT+"umi_bam/{sample}.sorted.aligned.umi.bai")
    shell:
        """
        gatk --java-options -Xmx{mem}G MergeBamAlignment \
        --ALIGNED {input.bam} \
        --UNMAPPED {input.unmapped_bam} \
        --OUTPUT {output.aligned_bam} \
        --REFERENCE_SEQUENCE {ref} \
        --SORT_ORDER coordinate \
        --ALIGNER_PROPER_PAIR_FLAGS TRUE \
        --ALIGNED_READS_ONLY TRUE \
        --CREATE_INDEX TRUE \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --MAX_INSERTIONS_OR_DELETIONS -1
        """

# Step D. Group reads per UMI and filter
rule correct_UMIs:
    input:
        aligned_bam = OUTPUT+"umi_bam/{sample}.sorted.aligned.umi.bam",
        aligned_bai = OUTPUT+"umi_bam/{sample}.sorted.aligned.umi.bai"
    output:
        corrected_bam = OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.bam",
        corrected_bai = OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.bai",
        metrics = OUTPUT+"quality_control/{sample}_correctUMIs_metrics.txt"
    shell:
        """
        fgbio -Xmx{mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap CorrectUmis \
        -i {input.aligned_bam} \
        -o {output.corrected_bam} \
        -M {output.metrics} \
        -m 2 -d 2 \
        -U {expectedUMIs}
        """

rule groupReadsByUmi:
    input:
        corrected_bam = OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.bam",
        corrected_bai = OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.bai"
    output:
        grouped_bam = temp(OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.grouped.bam")
    shell:
        """
        fgbio -Xmx{mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap GroupReadsByUmi \
        -i {input.corrected_bam} \
        -o {output.grouped_bam} \
        -s paired -m 20
        """


rule callConsensusReads:
    input:
        grouped_bam = OUTPUT+"umi_bam/{sample}.sorted.aligned.corrected_umi.grouped.bam"
    output:
        consensus_unmapped_bam = temp(OUTPUT+"umi_bam/{sample}.umi.grouped.consensus.bam")
    shell:
        """
        fgbio -Xmx{mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap CallMolecularConsensusReads \
        -i {input.grouped_bam} \
        -o {output.consensus_unmapped_bam} \
        --error-rate-post-umi=25 --min-read=1
        """


rule filterConsensusReads:
    input:
        consensus_unmapped_bam = OUTPUT+"umi_bam/{sample}.umi.grouped.consensus.bam"
    output:
        filtered_unmapped_bam = temp(OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.bam"),
        filtered_unmapped_bai = temp(OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.bai")
    shell:
        """
        fgbio -Xmx{mem}g -XX:+AggressiveOpts -XX:+AggressiveHeap FilterConsensusReads \
        -i {input.consensus_unmapped_bam} \
        -o {output.filtered_unmapped_bam} \
        -r {ref} \
        -M 1 \
        -E 0.05 -e 0.1 \
        -N 30 -n 0.1
        """


# Step E. Align grouped and filtered reads
rule sort_filtered_bam:
    input:
        filtered_unmapped_bam = OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.bam",
        filtered_unmapped_bai = OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.bai"
    output:
        filtered_sorted_unmapped_bam = OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.sorted.bam"
    shell:
        """
        gatk --java-options -Xmx{mem}G SortSam \
        --INPUT {input.filtered_unmapped_bam} \
        --OUTPUT {output.filtered_sorted_unmapped_bam} \
        --VALIDATION_STRINGENCY LENIENT \
        --SORT_ORDER queryname \
        --TMP_DIR tmp \
        --CREATE_INDEX true
        """

rule samToFastq_final:
    input:
        filtered_sorted_unmapped_bam = OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.sorted.bam"
    output:
        consensus_filtered_f1 = OUTPUT+"fastq_minusUMI/{sample}_minusUMI_consensus_filtered_R1.fastq.gz",
        consensus_filtered_f2 = OUTPUT+"fastq_minusUMI/{sample}_minusUMI_consensus_filtered_R2.fastq.gz",
    shell:
        """
        gatk --java-options -Xmx{mem}G SamToFastq \
        --INPUT {input.filtered_sorted_unmapped_bam} \
        --FASTQ {output.consensus_filtered_f1} \
        --SECOND_END_FASTQ {output.consensus_filtered_f2}
        """

# Step C. Align the new FASTQ files (with removed UMIs) with Bowtie 2
rule bwa_mem_final:
    input:
        consensus_filtered_f1 = OUTPUT+"fastq_minusUMI/{sampleid}_{protocol}_{flowcell}_minusUMI_consensus_filtered_R1.fastq.gz",
        consensus_filtered_f2 = OUTPUT+"fastq_minusUMI/{sampleid}_{protocol}_{flowcell}_minusUMI_consensus_filtered_R2.fastq.gz",
    output:
        bam = OUTPUT+"sorted_bam/{sampleid}_{protocol}_{flowcell}.sorted.consensus.filtered.bam",
        # bai = temp(OUTPUT+"sorted_bam/{sampleid}_{protocol}_{flowcell}.sorted.consensus.filtered.bai")
    resources: cpus=9, mem=50000
    log:
        "logs/bwa_logs/{sampleid}_{protocol}_{flowcell}.bwa_mem.log"
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{protocol}_{flowcell}",
        rgpl = "illumina",
        rgpu = "{flowcell}"
    shell:
        """
        bwa mem {ref} {input} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        -M \
        -t {resources.cpus} | \
        gatk --java-options -Xmx{mem}G SortSam \
        --INPUT /dev/stdin \
        --OUTPUT {output.bam} \
        --VALIDATION_STRINGENCY LENIENT \
        --SORT_ORDER queryname \
        --TMP_DIR tmp \
        --CREATE_INDEX TRUE
        """

rule mergeBamAlignment_final:
    input:
        bam = OUTPUT+"sorted_bam/{sample}.sorted.consensus.filtered.bam",
        filtered_sorted_unmapped_bam = OUTPUT+"umi_bam/{sample}.umi.consensus.filtered.sorted.bam"
    output:
        bam = OUTPUT+"bam/{sample}.umi.consensus.aligned.bam",
        bai = OUTPUT+"bam/{sample}.umi.consensus.aligned.bai"
    shell:
        """
        gatk --java-options -Xmx{mem}G MergeBamAlignment \
        --ALIGNED {input.bam} \
        --UNMAPPED {input.filtered_sorted_unmapped_bam} \
        --OUTPUT {output.bam} \
        --REFERENCE_SEQUENCE {ref} \
        --SORT_ORDER coordinate \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --ALIGNED_READS_ONLY true \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --MAX_INSERTIONS_OR_DELETIONS -1
        """

###############################################################################
####                                Metrics                                ####
###############################################################################
rule CollectHsMetrics:
    input:
        bam = OUTPUT+"bam/{sample}.umi.consensus.aligned.bam",
        bai = OUTPUT+"bam/{sample}.umi.consensus.aligned.bai",
        interval = {interval_list}
    output:
        OUTPUT+"quality_control/{sample}_capture-targets_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectHsMetrics \
        --INPUT {input.bam} \
        --REFERENCE_SEQUENCE {ref} \
        --OUTPUT {output} \
        --BAIT_INTERVALS {input.interval} \
        --TARGET_INTERVALS 	{input.interval}
        """

rule CollectInsertSizeMetrics:
    input:
        bam = OUTPUT+"bam/{sample}.umi.consensus.aligned.bam"
    output:
        hist = OUTPUT+"quality_control/{sample}.insert_size_histogram.pdf",
        metrics = OUTPUT+"quality_control/{sample}.insert_size_metrics.txt"
    shell:
        """
        gatk --java-options -Xmx{mem}G CollectInsertSizeMetrics \
        -I {input.bam} \
        -H {output.hist} \
        -O {output.metrics}
        """


rule qualimap:
    input:
        bam = OUTPUT+"bam/{sample}.umi.consensus.aligned.bam"
    output:
        html = OUTPUT+"quality_control/{sample}/qualimapReport.html"
    params:
        outdir = OUTPUT+"quality_control/{sample}"
    resources: cpus=6, mem=30000
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt {resources.cpus} \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size={resources.mem}
        """

rule multiqc:
    input:
        expand("quality_control/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand(OUTPUT+"quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand(OUTPUT+"quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand(OUTPUT+"quality_control/{sample}_correctUMIs_metrics.txt", sample=SAMPLES)
    output:
        html = OUTPUT+"quality_control/takara_umi_multiqc_report.html"
    params:
        indir1 = "quality_control/fastqc/",
        indir2 = OUTPUT+"quality_control/",
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir1} {params.indir2} \
        -n {output.html} \
        -c {params.config}
        """
