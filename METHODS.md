# Methods

## Data Source

Paired-end cfDNA sequencing data were downloaded from NCBI SRA (BioProject PRJNA291561; Snyder et al., 2016, *Cell*). Two samples were retained for analysis: SRR2130005 and SRR2130052. A third sample (SRR2130040) was excluded due to anomalously low mapping rate on the chromosome 1 subset (<40%), likely a quality artefact amplified by subsetting. Raw FASTQ files are not included in this repository.

## Computational Environment

All steps were executed using the conda environment specified in `envs/environment.yml`. Key tool versions:

| Tool | Version | Purpose |
|---|---|---|
| fastp | 0.23.4 | Adapter trimming and QC |
| bwa | 0.7.17 | Short-read alignment |
| samtools | 1.18 | BAM processing and QC |
| R | 4.3.x | Statistical analysis |
| ggplot2 | 3.5.x | Visualisation |
| dplyr | 1.1.x | Data manipulation |
| Snakemake | 7.32.x | Workflow management |

## Reference Genome

Chromosome 1 of the human genome assembly GRCh38/hg38 was used as the alignment reference:

```
Source : UCSC Genome Browser
URL    : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
MD5    : (run md5sum ref/chr1.fa after download to verify)
```

Alignment to chr1 only was chosen to make the pipeline fully reproducible on a standard laptop without a computing cluster. This is a deliberate educational constraint and does not reflect a genome-wide analysis.

## Read Subsetting

To further reduce computational requirements, each sample was processed using the first 500,000 read pairs (`fastp --reads_to_process 500,000`). This subsetting is applied at the trimming stage and propagates through all downstream steps. All downstream results therefore reflect a subsetted, chr1-limited dataset.

## Adapter Trimming

Adapters were auto-detected and removed using `fastp` in paired-end mode. Reads with mean quality < Q20 or length < 30 bp after trimming were discarded. Per-sample HTML and JSON QC reports are written to `trimmed/`.

```
fastp --detect_adapter_for_pe --average_qual 20 --length_required 30 \
      --reads_to_process 500000 --thread 4
```

## Alignment

Trimmed reads were aligned to `ref/chr1.fa` using `bwa mem` with default parameters and 4 threads. Alignment output was written in SAM format and piped directly into `samtools sort` to produce a coordinate-sorted BAM. The BAM was indexed with `samtools index`.

## Alignment Quality Control

Alignment statistics were generated with `samtools flagstat` for each sample. Output files are committed to `results/qc/`. Mapping rate, properly paired rate, and total read counts are the primary QC metrics reviewed.

## Fragment Length Extraction

Fragment lengths were derived from the TLEN (template length) field of the BAM file, which encodes the inferred insert size for paired-end reads. Only properly paired reads (SAM flag `0x2`) were retained, and only positive TLEN values in the range 1–999 bp were used, consistent with the expected cfDNA fragment size range.

```bash
samtools view -f 0x2 aligned/${SAMPLE}.sorted.bam \
  | awk '$9 > 0 && $9 < 1000 {print $9}' \
  > results/fragment_lengths/${SAMPLE}_fragment_lengths.txt
```

## Statistical Analysis and Visualisation

Fragment length distributions were analysed in R. For each sample, the following summary statistics were calculated: total fragment count, mean, median, standard deviation, and proportion of fragments in the mononucleosomal range (120–200 bp). Results are written to `results/fragment_lengths/summary_stats.csv`.

Two figures were generated:

1. **Histogram** (`fragment_histogram.png`): overlaid fragment length distributions for both samples, with a dashed reference line at 167 bp (canonical mononucleosome peak).
2. **Boxplot** (`fragment_boxplot.png`): side-by-side boxplots enabling comparison of distribution spread and median fragment length between samples.

All visualisation code is in `scripts/plot_fragment.R`, which reads the fragment length text files produced in the extraction step above. A clearly labelled simulation fallback is included for cases where pipeline outputs are unavailable (e.g. CI environments), but all committed figures are derived from real pipeline output.

## Reproducibility

The full pipeline is encoded in `Snakefile` and can be re-executed with `snakemake --cores 4` after setting up the conda environment and downloading the reference and raw data (see README). All intermediate files (BAMs, SAMs, FASTQs) are excluded from version control via `.gitignore`. Only small text results files and final figures are committed.

## References

Snyder, M. W., Kircher, M., Hill, A. J., Daza, R. M., & Shendure, J. (2016). Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell*, 164(1-2), 57–68. https://doi.org/10.1016/j.cell.2015.11.050
