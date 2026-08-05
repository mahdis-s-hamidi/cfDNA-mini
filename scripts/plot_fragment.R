# =============================================================
# cfDNA-mini | Fragment Length Analysis
# Author : Mahdis Saffar-Hamidi
# Contact: mahdiisshamidi79@gmail.com
# =============================================================
#
# Purpose:
#   Read fragment-length files generated from cfDNA BAM alignments
#   using the TLEN field from properly paired reads and generate:
#
#   1. Per-sample fragment length histograms
#   2. Cross-sample boxplot comparison
#   3. Summary statistics table
#
# Input:
#   results/fragment_lengths/fragment_lengths_0005.txt
#   results/fragment_lengths/fragment_lengths_0052.txt
#
# Output:
#   results/figures/
#       fragment_hist_0005.png
#       fragment_hist_0052.png
#       fragment_boxplot.png
#
#   results/fragment_lengths/
#       fragment_summary.csv
#
# Note on failure behaviour:
#   This script does NOT fall back to simulated data. If an expected
#   input file is missing, it stops with an explicit error pointing
#   to the upstream pipeline step that should have produced it. This
#   is intentional: a silently-generated fake plot is worse than a
#   clear failure for a reproducibility-focused portfolio project.
# =============================================================


library(dplyr)
library(ggplot2)


# -------------------------------------------------------------
# 0. Paths
# -------------------------------------------------------------

frag_dir <- "results/fragment_lengths"
fig_dir  <- "results/figures"

dir.create(fig_dir,
           recursive = TRUE,
           showWarnings = FALSE)


samples <- c("0005", "0052")


# -------------------------------------------------------------
# 1. Load fragment length data
# -------------------------------------------------------------

load_fragment_lengths <- function(sample_id) {

  path <- file.path(
    frag_dir,
    paste0("fragment_lengths_", sample_id, ".txt")
  )

  if (!file.exists(path)) {

    stop(
      "Missing fragment length file: ",
      path,
      "\nPlease run the upstream pipeline first ",
      "(see README: 'Run the full pipeline')."
    )

  }

  lengths <- scan(path, quiet = TRUE)

  lengths <- lengths[
    lengths > 0 &
    lengths < 1000
  ]

  message(
    "Loaded ",
    path,
    " (n = ",
    length(lengths),
    ")"
  )

  # Reconstruct the full SRA accession from the 4-digit code for use
  # in plot titles, axis labels, and the summary table. NOTE: the
  # accession prefix is "SRR213" (not "SRR21300") -- e.g. code
  # "0005" -> "SRR2130005", not "SRR213000005".
  data.frame(
    sample = paste0("SRR213", sample_id),
    length = lengths
  )

}


df <- bind_rows(
  lapply(samples, load_fragment_lengths)
)


# -------------------------------------------------------------
# 2. Summary statistics
# -------------------------------------------------------------

stats <- df %>%
  group_by(sample) %>%
  summarise(
    n = n(),
    mean_bp = round(mean(length), 1),
    median_bp = median(length),
    sd_bp = round(sd(length), 1),
    pct_mono_120_200bp =
      round(
        mean(length >= 120 &
             length <= 200) * 100,
        1
      ),
    .groups = "drop"
  )


print(stats)


write.csv(
  stats,
  file.path(
    frag_dir,
    "fragment_summary.csv"
  ),
  row.names = FALSE
)


# -------------------------------------------------------------
# 3. Per-sample histograms
# -------------------------------------------------------------

for (s in unique(df$sample)) {

  sample_df <- df %>%
    filter(sample == s)

  p <- ggplot(
    sample_df,
    aes(x = length)
  ) +

    geom_histogram(
      bins = 80,
      fill = "#0072B2",
      alpha = 0.8
    ) +

    geom_vline(
      xintercept = 167,
      linetype = "dashed",
      linewidth = 0.8
    ) +

    annotate(
      "text",
      x = 180,
      y = Inf,
      label = "167 bp\nmono-nucleosome",
      vjust = 1.5,
      size = 3.5
    ) +

    scale_x_continuous(
      limits = c(50, 600),
      breaks = seq(50, 600, 50)
    ) +

    labs(
      title =
        paste0(
          "cfDNA Fragment Length Distribution - ",
          s
        ),
      subtitle =
        "chr1 subset | TLEN from properly paired reads",
      x = "Fragment Length (bp)",
      y = "Read count",
      caption =
        "Data: SRA PRJNA291561"
    ) +

    theme_classic(
      base_size = 13
    )

  # Output filename uses the trailing 4-digit code, matching the
  # committed results/figures/ naming convention.
  output_name <- paste0(
    "fragment_hist_",
    substr(
      s,
      nchar(s) - 3,
      nchar(s)
    ),
    ".png"
  )

  ggsave(
    filename =
      file.path(
        fig_dir,
        output_name
      ),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )

}


# -------------------------------------------------------------
# 4. Boxplot comparison
# -------------------------------------------------------------

p_box <- ggplot(
  df,
  aes(
    x = sample,
    y = length
  )
) +

  geom_boxplot(
    outlier.size = 0.4,
    outlier.alpha = 0.3,
    width = 0.5
  ) +

  geom_hline(
    yintercept = 167,
    linetype = "dashed",
    linewidth = 0.8
  ) +

  scale_y_continuous(
    limits = c(50, 600),
    breaks = seq(50, 600, 50)
  ) +

  labs(
    title =
      "cfDNA Fragment Length Comparison",
    subtitle =
      "Dashed line represents canonical ~167 bp nucleosomal peak",
    x = NULL,
    y = "Fragment Length (bp)",
    caption =
      "Data: SRA PRJNA291561 | chr1 subset"
  ) +

  theme_classic(
    base_size = 13
  )


ggsave(
  filename =
    file.path(
      fig_dir,
      "fragment_boxplot.png"
    ),
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 300
)


message(
  "\nAnalysis completed successfully."
)

message(
  "Figures saved in: ",
  fig_dir
)

message(
  "Summary saved in: ",
  file.path(
    frag_dir,
    "fragment_summary.csv"
  )
)
