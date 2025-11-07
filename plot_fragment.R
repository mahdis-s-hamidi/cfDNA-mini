# Simulated cfDNA Fragment Size Distribution
# Learning exercise inspired by cfDNA studies in clinical genomics
set.seed(123)
fragment_lengths <- c(rnorm(8000, mean = 167, sd = 10),
                      rnorm(2000, mean = 340, sd = 15))

library(ggplot2)
df <- data.frame(length = fragment_lengths)

ggplot(df, aes(x = length)) +
  geom_histogram(bins = 60, fill = "#1f77b4", alpha = 0.8, color = "black") +
  geom_vline(xintercept = 167, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Simulated cfDNA Fragmentation",
       subtitle = "Peak at 167 bp (mononucleosome)",
       x = "Fragment Length (bp)", y = "Count") +
  theme_minimal()

ggsave("fragment_plot.png", width = 8, height = 5, dpi = 300)
