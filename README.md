```markdown
# cfDNA-mini

This repository is a mini project for processing cfDNA sequencing data. The goal is to provide a step-by-step pipeline for cfDNA data analysis and a reproducible example on GitHub.

## Project Structure

```

cfDNA-mini/
├── data_trimmed/
│   └── subset/               # subset fastq files
├── refs/                     # reference files
├── bam/                      # BAM outputs and sorted files
├── scripts/                  # scripts used (if any)
└── README.md

```

## Pipeline Steps

1. **Trimming**  
   - Using `fastp` to remove adapters and low-quality bases.  
   - Output files: `*_1.sub.fastq.gz` and `*_2.sub.fastq.gz`

2. **Alignment**  
   - Using `bwa mem` to align reads to the reference (e.g., chr1).  
   - Outputs: `.sam` and `.sorted.bam`

3. **BAM Processing**  
   - Convert `.sam` to `.bam` with `samtools view`  
   - Sort with `samtools sort`  
   - Index with `samtools index`

4. **QC & Summary**  
   - Check mapping and pairing with `samtools flagstat`  
   - Extract fragment lengths for insert size distribution analysis

## Notes

- This is a mini project for learning and creating a reproducible example on GitHub.  
- Subsets of the data were used for faster processing.  
- All steps are fully reproducible and can be applied to other samples.
```


> Contact: mahdiisshamidi79@gmail.com | [LinkedIn](https://linkedin.com/in/mahdis-saffar-hamidi)
