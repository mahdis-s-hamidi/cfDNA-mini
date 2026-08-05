# =============================================================
# cfDNA-mini  |  Snakemake Pipeline
# =============================================================
# Run:  snakemake --cores 4
# Env:  conda activate cfdna-mini
# =============================================================


SAMPLES = [
    "SRR2130005",
    "SRR2130052"
]

READS_TO_PROCESS = 500000


# Short codes are used only for final fragment-length outputs.
# Upstream tools require full SRA accession identifiers.

SAMPLE_CODE = {
    "SRR2130005": "0005",
    "SRR2130052": "0052",
}

CODE_TO_SAMPLE = {
    code: sample
    for sample, code in SAMPLE_CODE.items()
}

CODES = list(CODE_TO_SAMPLE.keys())


# =============================================================
# Final targets
# =============================================================

rule all:
    input:

        # QC outputs
        expand(
            "results/qc/{sample}.flagstat",
            sample=SAMPLES
        ),

        # Fragment length extraction
        expand(
            "results/fragment_lengths/fragment_lengths_{code}.txt",
            code=CODES
        ),

        # Summary table
        "results/fragment_lengths/fragment_summary.csv",

        # Figures
        expand(
            "results/figures/fragment_hist_{code}.png",
            code=CODES
        ),

        "results/figures/fragment_boxplot.png"



# =============================================================
# 1. Trim reads
# =============================================================

rule trim:

    input:
        r1 = "data_raw/{sample}_1.fastq.gz",
        r2 = "data_raw/{sample}_2.fastq.gz"

    output:
        r1 = "trimmed/{sample}_1.sub.fastq.gz",
        r2 = "trimmed/{sample}_2.sub.fastq.gz",
        json = "trimmed/{sample}_fastp.json",
        html = "trimmed/{sample}_fastp.html"

    threads:
        4

    log:
        "logs/fastp/{sample}.log"

    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --reads_to_process {READS_TO_PROCESS} \
            --detect_adapter_for_pe \
            --average_qual 20 \
            --length_required 30 \
            --thread {threads} \
            -j {output.json} \
            -h {output.html} \
            2> {log}
        """



# =============================================================
# 2. Alignment
# =============================================================

rule align:

    input:
        r1 = "trimmed/{sample}_1.sub.fastq.gz",
        r2 = "trimmed/{sample}_2.sub.fastq.gz",
        ref = "ref/chr1.fa"

    output:
        bam = "aligned/{sample}.sorted.bam",
        bai = "aligned/{sample}.sorted.bam.bai"

    threads:
        4

    log:
        "logs/bwa/{sample}.log"

    shell:
        """
        bwa mem \
            -t {threads} \
            {input.ref} \
            {input.r1} \
            {input.r2} \
            2> {log} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index {output.bam}
        """



# =============================================================
# 3. Mapping QC
# =============================================================

rule flagstat:

    input:
        bam = "aligned/{sample}.sorted.bam"

    output:
        "results/qc/{sample}.flagstat"

    shell:
        """
        samtools flagstat \
            {input.bam} \
            > {output}
        """



# =============================================================
# 4. Extract fragment lengths
# =============================================================

# Uses short codes for output filenames:
#
# fragment_lengths_0005.txt
# fragment_lengths_0052.txt
#
# BAM files still use full SRA accessions.

rule fragment_lengths:

    input:

        bam = lambda wildcards:
            f"aligned/{CODE_TO_SAMPLE[wildcards.code]}.sorted.bam",

        bai = lambda wildcards:
            f"aligned/{CODE_TO_SAMPLE[wildcards.code]}.sorted.bam.bai"


    output:

        "results/fragment_lengths/fragment_lengths_{code}.txt"


    shell:
        """
        samtools view -f 0x2 {input.bam} \
        | awk '$9 > 0 && $9 < 1000 {{print $9}}' \
        > {output}
        """



# =============================================================
# 5. Fragment length analysis and visualization
# =============================================================

# The R script consumes both fragment length files
# and generates:
#
# results/figures/
#   fragment_hist_0005.png
#   fragment_hist_0052.png
#   fragment_boxplot.png
#
# results/fragment_lengths/
#   fragment_summary.csv


rule plot:

    input:

        expand(
            "results/fragment_lengths/fragment_lengths_{code}.txt",
            code=CODES
        )


    output:

        hist_0005 =
            "results/figures/fragment_hist_0005.png",

        hist_0052 =
            "results/figures/fragment_hist_0052.png",

        box =
            "results/figures/fragment_boxplot.png",

        stats =
            "results/fragment_lengths/fragment_summary.csv"


    log:

        "logs/R/plot_fragment.log"


    shell:

        """
        Rscript scripts/plot_fragment.R \
        2> {log}
        """
