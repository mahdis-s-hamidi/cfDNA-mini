# =============================================================
# cfDNA-mini  |  Snakemake Pipeline
# =============================================================
# Run:  snakemake --cores 4
# Env:  conda activate cfdna-mini
# =============================================================

SAMPLES = ["SRR2130005", "SRR2130052"]
READS_TO_PROCESS = 500000

rule all:
    input:
        expand("results/qc/{sample}.flagstat",              sample=SAMPLES),
        expand("results/fragment_lengths/{sample}_fragment_lengths.txt", sample=SAMPLES),
        "results/fragment_lengths/summary_stats.csv",
        "results/figures/fragment_histogram.png",
        "results/figures/fragment_boxplot.png"


# ── 1. Trim ───────────────────────────────────────────────────
rule trim:
    input:
        r1 = "data_raw/{sample}_1.fastq.gz",
        r2 = "data_raw/{sample}_2.fastq.gz"
    output:
        r1 = "trimmed/{sample}_1.sub.fastq.gz",
        r2 = "trimmed/{sample}_2.sub.fastq.gz",
        json = "trimmed/{sample}_fastp.json",
        html = "trimmed/{sample}_fastp.html"
    threads: 4
    log: "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --reads_to_process {READS_TO_PROCESS} \
            --detect_adapter_for_pe \
            --average_qual 20 \
            --length_required 30 \
            --thread {threads} \
            -j {output.json} -h {output.html} \
            2> {log}
        """


# ── 2. Align ──────────────────────────────────────────────────
rule align:
    input:
        r1 = "trimmed/{sample}_1.sub.fastq.gz",
        r2 = "trimmed/{sample}_2.sub.fastq.gz",
        ref = "ref/chr1.fa"
    output:
        bam = "aligned/{sample}.sorted.bam",
        bai = "aligned/{sample}.sorted.bam.bai"
    threads: 4
    log: "logs/bwa/{sample}.log"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


# ── 3. Flagstat ───────────────────────────────────────────────
rule flagstat:
    input:
        bam = "aligned/{sample}.sorted.bam"
    output:
        "results/qc/{sample}.flagstat"
    shell:
        "samtools flagstat {input.bam} > {output}"


# ── 4. Extract fragment lengths ───────────────────────────────
rule fragment_lengths:
    input:
        bam = "aligned/{sample}.sorted.bam",
        bai = "aligned/{sample}.sorted.bam.bai"
    output:
        "results/fragment_lengths/{sample}_fragment_lengths.txt"
    shell:
        """
        samtools view -f 0x2 {input.bam} \
            | awk '$9 > 0 && $9 < 1000 {{print $9}}' \
            > {output}
        """


# ── 5. R analysis & plots ─────────────────────────────────────
rule plot:
    input:
        expand("results/fragment_lengths/{sample}_fragment_lengths.txt",
               sample=SAMPLES)
    output:
        hist   = "results/figures/fragment_histogram.png",
        box    = "results/figures/fragment_boxplot.png",
        stats  = "results/fragment_lengths/summary_stats.csv"
    log: "logs/R/plot_fragment.log"
    shell:
        "Rscript scripts/plot_fragment.R 2> {log}"
