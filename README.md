# cfDNA Fragmentation Mini-Project (Self-Learning)

**Mahdis Saffar Hamidi** | MSc Genetics | November 2025

> **Goal**: Practice and learning project to explore cfDNA fragmentation patterns, inspired by ongoing research in molecular diagnostics (e.g., cfDNA studies at LMU Klinikum).  
> **Dataset**:  [GSE100684 – cfDNA in cancer patients](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100684)   
> **Tools**: R (ggplot2), FastQC, TrimGalore, BWA  
> **Status**: Simulation completed | Real data analysis planned (Dec 2025)

---

## Background

Cell-free DNA (cfDNA) refers to DNA fragments circulating freely in the bloodstream. Its fragmentation patterns reflect nucleosome positioning and can inform tissue-of-origin, making cfDNA a powerful tool in non-invasive diagnostics. This mini-project simulates cfDNA fragmentation and prepares a pipeline for analyzing real sequencing data from cancer patients.

---

## Results (Simulated)
- Mean fragment length: ~167 bp  
- Nucleosome peak: 166–170 bp  

![Fragment Size Distribution](fragment_plot.png)

---

## R Script
```r
# plot_fragment.R - Simulated cfDNA fragmentation
set.seed(123)
lengths <- c(rnorm(8000, 167, 10), rnorm(2000, 340, 15))
library(ggplot2)
ggplot(data.frame(length = lengths), aes(x = length)) +
  geom_histogram(bins = 60, fill = "#1f77b4", alpha = 0.8) +
  labs(title = "Simulated cfDNA Fragmentation", x = "Length (bp)", y = "Count") +
  theme_minimal()
ggsave("fragment_plot.png", width = 8, height = 5, dpi = 300)
```

---

## Pipeline (In Progress)
```bash
fastqc SRR10990278.fastq
trim_galore SRR10990278.fastq
bwa mem hg38.fa SRR10990278_trimmed.fq > aligned.sam
```

> Contact: mahdiisshamidi79@gmail.com | [LinkedIn](https://linkedin.com/in/mahdis-saffar-hamidi)
