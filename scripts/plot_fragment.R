# =============================================================
# cfDNA-mini | Fragment Length Analysis
# Author : Mahdis Saffar-Hamidi
# Contact: mahdiisshamidi79@gmail.com
# =============================================================
# Purpose:
#   Read insert-size summary files produced by samtools view
#   (properly paired reads, TLEN column) and generate:
#     1. Histogram comparison of fragment length distributions
#     2. Boxplot comparison across samples
#
# Input  : results/fragment_lengths/<sample>_fragment_lengths.txt
#          One integer (TLEN) per line, positive values only.
#
# Output : results/figures/fragment_histogram.png
#          results/figures/fragment_boxplot.png
#          results/fragment_lengths/summary_stats.csv
# =============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# ── 0.  Paths ─────────────────────────────────────────────────
frag_dir  <- "results/fragment_lengths"
fig_dir   <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("SRR2130005", "SRR2130052")

# ── 1.  Load data ─────────────────────────────────────────────
#
# NOTE ON DATA SOURCE
# -------------------
# Fragment lengths were extracted from properly paired alignments
# using the TLEN field of the BAM file:
#
#   samtools view -f 0x2 aligned/<sample>.sorted.bam \
#     | awk '$9 > 0 && $9 < 1000 {print $9}' \
#     > results/fragment_lengths/<sample>_fragment_lengths.txt
#
# The code below reads those real pipeline outputs.
# If the files are absent (e.g. in a GitHub preview / CI run
# without BAM data), a clearly labelled SIMULATION is used
# so the plotting logic can still be inspected.  All simulated
# panels are annotated with "[SIMULATED — not real data]".

load_or_simulate <- function(sample) {
  path <- file.path(frag_dir, paste0(sample, "_fragment_lengths.txt"))
  if (file.exists(path)) {
    lengths <- scan(path, quiet = TRUE)
    lengths <- lengths[lengths > 0 & lengths < 1000]
    message("Loaded real data for ", sample, " (n = ", length(lengths), ")")
    data.frame(sample = sample, length = lengths, source = "real")
  } else {
    warning("File not found: ", path, "\nUsing SIMULATED data for display only.")
    set.seed(42)
    sim <- c(rnorm(8000, mean = 167, sd = 10),
             rnorm(2000, mean = 340, sd = 15))
    sim <- round(sim[sim > 0 & sim < 1000])
    data.frame(sample = paste0(sample, " [SIMULATED]"),
               length = sim,
               source = "simulated")
  }
}

df <- bind_rows(lapply(samples, load_or_simulate))

has_simulated <- any(df$source == "simulated")

# ── 2.  Summary statistics ────────────────────────────────────
stats <- df %>%
  group_by(sample) %>%
  summarise(
    n        = n(),
    mean_bp  = round(mean(length), 1),
    median_bp = median(length),
    sd_bp    = round(sd(length), 1),
    pct_mono = round(mean(length >= 120 & length <= 200) * 100, 1),
    .groups  = "drop"
  )

print(stats)
write.csv(stats,
          file = file.path(frag_dir, "summary_stats.csv"),
          row.names = FALSE)

# ── 3.  Colour palette (colour-blind friendly) ────────────────
pal <- c("#0072B2", "#E69F00")
names(pal) <- unique(df$sample)

sim_note <- if (has_simulated)
  "\nNote: samples labelled [SIMULATED] do not reflect real BAM data."
else ""

# ── 4.  Histogram ─────────────────────────────────────────────
p_hist <- ggplot(df, aes(x = length, fill = sample)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity",
                 colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 167, colour = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 172, y = Inf, label = "167 bp\n(mono-nucleosome)",
           hjust = 0, vjust = 1.4, size = 3.2, colour = "grey30") +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(50, 600), breaks = seq(50, 600, 50)) +
  labs(
    title    = "cfDNA Fragment Length Distribution (chr1 subset)",
    subtitle = paste0("SRR2130005 vs SRR2130052  |  bwa mem → samtools TLEN", sim_note),
    x        = "Fragment Length (bp)",
    y        = "Read count",
    fill     = "Sample",
    caption  = "Data: SRA (PRJNA291561). chr1 only; subsetted reads."
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        plot.caption = element_text(colour = "grey50", size = 9))

ggsave(file.path(fig_dir, "fragment_histogram.png"),
       p_hist, width = 9, height = 5, dpi = 300)

# ── 5.  Boxplot ───────────────────────────────────────────────
p_box <- ggplot(df, aes(x = sample, y = length, fill = sample)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
               linewidth = 0.6, width = 0.5) +
  geom_hline(yintercept = 167, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(50, 600), breaks = seq(50, 600, 50)) +
  labs(
    title   = "Fragment Length Comparison (chr1 subset)",
    subtitle = paste0("Dashed line = 167 bp mono-nucleosome peak", sim_note),
    x       = NULL,
    y       = "Fragment Length (bp)",
    caption = "Data: SRA (PRJNA291561). chr1 only; subsetted reads."
  ) +
  theme_classic(base_size = 13) +
  theme(plot.caption = element_text(colour = "grey50", size = 9))

ggsave(file.path(fig_dir, "fragment_boxplot.png"),
       p_box, width = 7, height = 5, dpi = 300)

message("\nDone. Figures written to: ", fig_dir)
message("Summary stats written to: ", frag_dir, "/summary_stats.csv")
