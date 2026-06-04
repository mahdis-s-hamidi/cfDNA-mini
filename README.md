# cfDNA-mini

> **A reproducible educational pipeline for cell-free DNA (cfDNA) fragment-length analysis**  
> Portfolio project · Mahdis Saffar-Hamidi · [LinkedIn](https://linkedin.com/in/mahdis-saffar-hamidi) · mahdiisshamidi79@gmail.com

---

## Overview

Cell-free DNA (cfDNA) circulates in plasma as short nucleosome-protected fragments, typically ~167 bp (mononucleosome) and ~340 bp (dinucleosome). Deviations from these characteristic lengths are a hallmark of cancer and other diseases, making cfDNA fragment-length profiling a key method in liquid biopsy research.

This repository demonstrates a **complete, reproducible processing workflow** from raw paired-end sequencing reads to fragment-length visualisation. It is intentionally scoped as an educational project — using chromosome 1 only and subsetted reads — so that every step can be run on a laptop without a computing cluster.

> **Scope disclaimer:** This is a portfolio/learning project, not a publication. Subsetted reads and single-chromosome alignment are used for computational feasibility. Results are illustrative, not scientifically conclusive.

---

## Samples

| Accession | BioProject | Status | Notes |
|---|---|---|---|
| SRR2130005 | PRJNA291561 | ✅ Used | Healthy donor cfDNA |
| SRR2130052 | PRJNA291561 | ✅ Used | Healthy donor cfDNA |
| SRR2130040 | PRJNA291561 | ❌ Discarded | Excluded due to low mapping rate on chr1 subset (<40%); likely a data quality issue amplified by subsetting |

Raw data are available from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra?term=PRJNA291561). They are **not** included in this repository to keep it lightweight.

---

## Pipeline Steps

```
FASTQ (SRA)  →  fastp (trim)  →  bwa mem (align, chr1)  →  samtools (BAM)  →  R (fragment analysis)
```

### 1. Read Trimming — `fastp`
Adapter sequences and low-quality bases (Q < 20) are removed. A read subset (`--reads_to_process 500000`) is used to reduce runtime.

```bash
fastp \
  -i data_raw/${SAMPLE}_1.fastq.gz \
  -I data_raw/${SAMPLE}_2.fastq.gz \
  -o trimmed/${SAMPLE}_1.sub.fastq.gz \
  -O trimmed/${SAMPLE}_2.sub.fastq.gz \
  --reads_to_process 500000 \
  --detect_adapter_for_pe \
  --thread 4 \
  -j trimmed/${SAMPLE}_fastp.json \
  -h trimmed/${SAMPLE}_fastp.html
```

### 2. Alignment — `bwa mem`
Reads are aligned to chromosome 1 of hg38. Aligning to a single chromosome keeps reference and output files small.

```bash
bwa mem -t 4 ref/chr1.fa \
  trimmed/${SAMPLE}_1.sub.fastq.gz \
  trimmed/${SAMPLE}_2.sub.fastq.gz \
  > aligned/${SAMPLE}.sam
```

### 3. BAM Processing — `samtools`

```bash
samtools view  -bS aligned/${SAMPLE}.sam | \
samtools sort  -o aligned/${SAMPLE}.sorted.bam
samtools index     aligned/${SAMPLE}.sorted.bam
```

### 4. QC — `samtools flagstat`

```bash
samtools flagstat aligned/${SAMPLE}.sorted.bam \
  > results/qc/${SAMPLE}.flagstat
```

### 5. Fragment Length Extraction

Insert sizes are taken from the TLEN field of properly paired (`-f 0x2`) reads:

```bash
samtools view -f 0x2 aligned/${SAMPLE}.sorted.bam \
  | awk '$9 > 0 && $9 < 1000 {print $9}' \
  > results/fragment_lengths/${SAMPLE}_fragment_lengths.txt
```

### 6. Statistical Analysis & Visualisation — R

```r
Rscript scripts/plot_fragment.R
```

Produces:
- `results/figures/fragment_histogram.png`
- `results/figures/fragment_boxplot.png`
- `results/fragment_lengths/summary_stats.csv`

See `scripts/plot_fragment.R` for full annotated source. The script reads real pipeline output files when present, and falls back to clearly labelled simulated data for display purposes when BAM-derived files are absent.

---

## Repository Structure

```
cfDNA-mini/
├── .gitignore
├── LICENSE
├── README.md
├── METHOD.md
├── Snakefile
├── environment.yml
├── plot_fragment.R
└── results/
    ├── qc/
    │   ├── SRR2130005.flagstat
    │   └── SRR2130052.flagstat
    ├── fragment_lengths/
    │   ├── SRR2130005_fragment_lengths.txt
    │   ├── SRR2130052_fragment_lengths.txt
    │   └── summary_stats.csv
    └── figures/
        ├── fragment_hist_0005.png
        ├── fragment_hist_0052.png
        └── fragment_boxplot.png

---

## Reproducing the Analysis

### Requirements

- Linux or macOS
- [conda](https://docs.conda.io/en/latest/miniconda.html) ≥ 23.x  
- ~2 GB disk space (reference + trimmed reads + BAMs for chr1)

### Setup

```bash
git clone https://github.com/mahdis-s-hamidi/cfDNA-mini.git
cd cfDNA-mini
conda env create -f envs/environment.yml
conda activate cfdna-mini
```

### Download raw data

```bash
mkdir -p data_raw
fastq-dump --split-files --gzip SRR2130005 -O data_raw/
fastq-dump --split-files --gzip SRR2130052 -O data_raw/
```

### Download reference (chr1 only)

```bash
mkdir -p ref
wget -O ref/chr1.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
gunzip ref/chr1.fa.gz
bwa index ref/chr1.fa
```

### Run the full pipeline

```bash
snakemake --cores 4
```

Or run step-by-step using the commands in each pipeline section above.

---

## Key Results

| Sample | n fragments | Mean (bp) | Median (bp) | SD (bp) |
|---|---|---|---|---|
| SRR2130005 | 195,889 | 142.8 | 151 | 57.6 |
| SRR2130052 | 97,293 | 145.7 | 153 | 55.3 |

Both samples show median fragment lengths close to the canonical mononucleosomal 
peak (~167 bp), consistent with healthy donor cfDNA. The distributions include a 
short-fragment population (<100 bp) likely representing adapter-derived or chimeric 
reads not fully removed by trimming — a known artefact in cfDNA subset analyses and 
a documented limitation of this pipeline. Mean values are pulled below the median 
by this short-fragment tail.

> Fragment lengths were extracted from chr1-aligned, properly paired reads (TLEN field).
> n = number of fragments per sample after filtering TLEN 1–999 bp.

### Fragment Length Distributions

![Fragment Size Distribution - SRR2130005](results/figures/fragment_hist_0005.png)
![Fragment Size Distribution - SRR2130052](results/figures/fragment_hist_0052.png)
![Fragment Length Boxplot Comparison](results/figures/fragment_boxplot.png)

### Biological Interpretation

The median fragment lengths of ~151–153 bp in both samples fall below the 
canonical 167 bp mononucleosomal peak. This is consistent with published reports 
of cfDNA being slightly shorter than the full nucleosome-protected length, due to 
exonucleolytic trimming at fragment ends during apoptosis (Snyder et al., 2016). 
The broad distribution and short-fragment tail reflect the chr1-subset and 
read-depth limitations of this educational pipeline.

## Background

cfDNA is released into circulation during apoptosis, where chromosomal DNA is cleaved at linker regions between nucleosomes. This produces a characteristic ladder of fragment sizes: ~167 bp (mono-), ~340 bp (di-), ~510 bp (tri-nucleosome). In cancer, this pattern shifts due to altered chromatin architecture and different cell-of-origin contributions.

Key references:
- Snyder et al. (2016). Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell*, 164(1-2), 57–68.
- Cristiano et al. (2019). Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature*, 570, 385–389.

---

## Limitations

- Only chromosome 1 was used for alignment; genome-wide patterns may differ.
- Read subsets (500 k reads per sample) reduce statistical power.
- This is a demonstration pipeline; results should not be interpreted clinically or scientifically without full-genome, full-depth processing.

---

## License

MIT License. See `LICENSE`.

