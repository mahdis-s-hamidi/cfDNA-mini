# Methods

## Data Source

Paired-end cfDNA sequencing data were obtained from NCBI SRA (BioProject PRJNA291561; Snyder et al., 2016, *Cell*).

Two samples were retained for analysis:

* SRR2130005
* SRR2130052

A third sample (SRR2130040) was excluded due to anomalously low mapping rate on the chromosome 1 subset (<40%), likely reflecting a data quality issue amplified by single-chromosome subsetting.

Raw FASTQ files are not included in this repository to keep the project lightweight.

---

## Computational Environment

All analyses were performed using the conda environment defined in:

```text
envs/environment.yml
```

Key software versions:

| Tool      | Version | Purpose                      |
| --------- | ------- | ---------------------------- |
| fastp     | 0.23.4  | Adapter trimming and read QC |
| bwa       | 0.7.17  | Short-read alignment         |
| samtools  | 1.18    | BAM processing and QC        |
| R         | 4.3.x   | Statistical analysis         |
| ggplot2   | 3.5.x   | Visualisation                |
| dplyr     | 1.1.x   | Data manipulation            |
| Snakemake | 7.32.x  | Workflow management          |

---

## Reference Genome

Chromosome 1 from the human genome assembly GRCh38/hg38 was used as the alignment reference.

```text
Source : UCSC Genome Browser
URL    : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
MD5    : (run md5sum ref/chr1.fa after download to verify)
```

Alignment was restricted to chromosome 1 to allow complete execution on a standard laptop without high-performance computing resources. This represents an educational simplification and does not reflect a genome-wide cfDNA analysis.

```

Source : UCSC Genome Browser
URL    : [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz)
MD5    : (run md5sum ref/chr1.fa after download to verify)

```

Alignment was restricted to chromosome 1 to allow complete execution on a standard laptop without high-performance computing resources. This represents an educational simplification and does not reflect a genome-wide cfDNA analysis.

---

## Read Subsetting

To reduce computational requirements, each sample was processed using the first 500,000 read pairs:

```

fastp --reads_to_process 500000

````

Subsetting was performed during trimming and propagated through all downstream steps. Therefore, all reported results represent a subsetted, chr1-limited dataset.

---

## Adapter Trimming

Adapter detection and removal were performed using `fastp` in paired-end mode.

Reads were filtered using:

- Mean quality threshold: Q20
- Minimum read length: 30 bp

Per-sample HTML and JSON QC reports are generated during the trimming step.

Example command:

```bash
fastp \
  --detect_adapter_for_pe \
  --average_qual 20 \
  --length_required 30 \
  --reads_to_process 500000 \
  --thread 4
````

---

## Alignment

Trimmed reads were aligned to `ref/chr1.fa` using `bwa mem` with 4 threads.

The resulting SAM output was converted to coordinate-sorted BAM format using `samtools sort`, followed by indexing with:

```bash
samtools index aligned/${SAMPLE}.sorted.bam
```

---

## Alignment Quality Control

Alignment statistics were generated using:

```bash
samtools flagstat aligned/${SAMPLE}.sorted.bam
```

Output files are stored in:

```
results/qc/
```

The primary QC metrics reviewed include:

* Total aligned reads
* Mapping rate
* Properly paired read percentage

---

## Fragment Length Extraction

Fragment lengths were calculated from the TLEN (template length) field of paired-end BAM alignments.

Only properly paired reads (`SAM flag 0x2`) were retained. Positive TLEN values between 1 and 999 bp were extracted to represent the expected cfDNA fragment-length range.

Example command:

```bash
samtools view -f 0x2 aligned/${SAMPLE}.sorted.bam \
  | awk '$9 > 0 && $9 < 1000 {print $9}' \
  > results/fragment_lengths/fragment_lengths_${CODE}.txt
```

Output filenames use the final four digits of each SRA accession:

```
SRR2130005 → fragment_lengths_0005.txt
SRR2130052 → fragment_lengths_0052.txt
```

---

## Statistical Analysis and Visualisation

Fragment length distributions were analysed using:

```
scripts/plot_fragment.R
```

The following summary statistics were calculated for each sample:

* Total fragment count
* Mean fragment length
* Median fragment length
* Standard deviation
* Percentage of fragments within the mononucleosomal range (120–200 bp)

Results are written to:

```
results/fragment_lengths/fragment_summary.csv
```

Three figures are generated:

1. **Per-sample fragment length histograms**

```
results/figures/fragment_hist_0005.png
results/figures/fragment_hist_0052.png
```

Each histogram displays the fragment-length distribution for one sample with a reference line at 167 bp, representing the canonical mononucleosome-associated fragment size.

2. **Cross-sample boxplot comparison**

```
results/figures/fragment_boxplot.png
```

The boxplot compares fragment-length distributions and median values between samples.

All visualisation code is contained in:

```
scripts/plot_fragment.R
```

The script reads fragment-length files generated by the upstream pipeline and stops with an explicit error if required inputs are unavailable. No simulated data are used for committed results.

---

## Reproducibility

The complete workflow is implemented in:

```
Snakefile
```

After creating the conda environment and downloading the required reference and raw sequencing data, the pipeline can be executed with:

```bash
snakemake --cores 4
```

Large intermediate files, including FASTQ files, BAM files, and alignment intermediates, are excluded from version control using `.gitignore`.

Committed repository outputs consist of:

* Fragment-length text files
* Summary statistics
* Final visualisations
* Pipeline and analysis scripts

---

## References

Snyder, M. W., Kircher, M., Hill, A. J., Daza, R. M., & Shendure, J. (2016).
Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin.
*Cell*, 164(1-2), 57–68.
[https://doi.org/10.1016/j.cell.2015.11.050](https://doi.org/10.1016/j.cell.2015.11.050)

```
```
