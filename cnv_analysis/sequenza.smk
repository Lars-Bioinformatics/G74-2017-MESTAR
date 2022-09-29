from snakemake.utils import R

WORK = "/work/Data/GRCh38.dragen.consensus/"
INPUT_BAM = WORK+"bam/"
OUTPUT = WORK+"sequenza/"

#configfile: WORK+"samples-pairwise.yaml" # RUN 1 - Medexome
configfile: WORK+"samples-pairwise-DRAGEN-MedExTwist.yaml" # RUN 2 - TwistMedEx

ref_build = "hg38"

#bam_suffix = ".connor.recalibrated.bam" # RUN 1 - Medexome
bam_suffix = ".GRCh38.dragen.consensus.umi-min-1.bam" # RUN 2 - TwistMedEx

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
    wig = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.gc50Base.wig.gz"
    chr_prefix = ""
else:
    # GRCh38 paths
    ref = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa"
    wig = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.gc50Base.wig.gz"
    chr_prefix = "chr"


# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G56-sampleA1_nimblegen-medexome_HGL2LDSXX"]
print(SAMPLES)

window = "50"

onstart:
    shell("mkdir -p "+OUTPUT)

rule all:
    input:
        [expand(OUTPUT+"Sequenza_{tumor}_vs_{normal}_bin"+window+"/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf", 
            tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in SAMPLES],



###############################################################################
#### Sequenza
###############################################################################
rule gc_wiggle:
    input:
        ref = ref
    output:
        wig = wig
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils gc_wiggle --fasta {input.ref} -w {window} -o {output.wig}
        """

rule bam2seqz:
    input:
        normal=INPUT_BAM+"{normal}"+bam_suffix,
        tumor=INPUT_BAM+"{tumor}"+bam_suffix,
        wig = wig
    output:
        seqz=OUTPUT+"{tumor}_vs_{normal}.seqz.gz"
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} \
        -q 20 -N 20 -f "illumina" --fasta {ref} -gc {wig} -o {output}
        """

rule seqz_binning:
    input:
        seqz=OUTPUT+"{tumor}_vs_{normal}.seqz.gz"
    output:
        seqz_bin=OUTPUT+"{tumor}_vs_{normal}.seqz.bin{window}.gz"
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils seqz_binning --seqz {input.seqz} -w {window} -o {output.seqz_bin}
        """

rule sequenza_analysis:
    input:
        seqz_bin=OUTPUT+"{tumor}_vs_{normal}.seqz.bin{window}.gz"
    output:
        res=OUTPUT+"Sequenza_{tumor}_vs_{normal}_bin{window}/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf"
    conda: "envs/sequenza.yaml"
    shell:
        """
        Rscript -e '
        library(sequenza)
        seqzdata = sequenza.extract("{input.seqz_bin}", chromosome.list = paste0("{chr_prefix}",c(1:22,"X","Y")))

        CP = sequenza.fit(seqzdata)
        
        id = "Sequenza_{wildcards.tumor}_vs_{wildcards.normal}"
        out.dir = paste0("{OUTPUT}", id, "_bin{wildcards.window}")
        sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir)
        '
        """
