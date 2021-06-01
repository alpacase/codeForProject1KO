library(tidyverse)

G1G2_1KO <- read_csv("G1G2_1KO.csv")

# data processing
d <- G1G2_1KO %>%
  mutate(
    bins = gsub(bins, pattern = "1KO_", replacement = ""),
    peak = gsub(peak, pattern = "\\..*Peak", replacement = ""),
    peak = gsub(peak, pattern = "atac.*", replacement = "ATAC-seq")
  ) %>%
  pivot_longer(cols = where(is.numeric), names_to = "pos") %>%
  separate(
    col = bins,
    sep = "_",
    remove = T,
    into = c("mod", "phs")
  ) %>%
  pivot_wider(names_from = phs, values_from = value) %>%
  mutate(
    pos = as.numeric(pos),
    log = log2(G1 / G2),
    mod = factor(mod, levels = c("5mC", "5hmC")),
    peak = factor(
      peak,
      levels = c("ATAC-seq", "H3K9ac", "H3K27ac", "H3K9me3", "H3K27me3")
    )
  )

# plot
p <- ggplot(d) +
  aes(x = pos, y = log, colour = mod) +
  geom_line(size = 0.5) +
  scale_x_continuous(breaks = c(0, 25, 50),
                     labels = c("-2.5", "0", "2.5")) +
  scale_y_continuous(
    limits = c(-0.1, 0.3),
    breaks = c(-0.1, 0, 0.1, 0.2, 0.3),
    n.breaks = 5
  ) +
  scale_color_manual(values = c("#00aa00", "#ff0000")) +
  labs(x = "distance of peak center (kb)", y = "log2(G1/G2)", color = "") +
  theme_user() +
  theme(aspect.ratio = 5 / 4) +
  facet_wrap(vars(peak), nrow = 1L, scales = "free_y")

ggsave(
  plot = p,
  filename = "fig2c.pdf",
  units = "mm",
  width = 120,
  height = 50
)
