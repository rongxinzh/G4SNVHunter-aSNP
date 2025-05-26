set.seed(1)

library(G4SNVHunter)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(scales)

seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
              names = paste0("chr", 1:22))
G4_gr <- G4HunterDetect(seq)
rm(seq)
gc()

all_var_files <- list.files("./Sprime_results/mendeley_data/",
                            "*.ND_match", full.names = TRUE)

variants_list <- list()
for (i in seq_along(all_var_files)) {
  tmp_file <- all_var_files[i]
  message(tmp_file)
  tmp_intvar <- fread(tmp_file, sep = "\t", header = TRUE) %>% data.frame()
  variants_list[[i]] <- tmp_intvar
}
variants <- do.call(rbind, variants_list)
variants <- variants %>% filter(if_any(c(NMATCH, DMATCH), ~ . == "match"), ALLELE == 1) %>% 
  select(CHROM, POS, REF, ALT, NMATCH, DMATCH ) %>% data.frame() %>% unique() 

variants <- GRanges(seqnames = paste0("chr", variants$CHROM),
                    ranges = IRanges(start = variants$POS, width = nchar(variants$REF)),
                    strand = "*",
                    ref = variants$REF,
                    alt = variants$ALT,
                    neandertal = variants$NMATCH,
                    denisovan = variants$DMATCH)

result <- G4VarImpact(G4 = G4_gr,
                      variants = variants,
                      ref_col = "ref",
                      alt_col = "alt")

filtered_mutG4s <- filterVarImpact(result, mut_score_threshold = 1.2)

uniq_mutG4 <- data.frame(
  wild_score = filtered_mutG4s$G4.info.max_score,
  mut_score = filtered_mutG4s$mutated.max_score,
  ref = filtered_mutG4s$variant.info.ref,
  introgressed = ifelse((filtered_mutG4s$variant.info.neandertal == "match") & (filtered_mutG4s$variant.info.denisovan != "match"),
                        "neandertal", ifelse((filtered_mutG4s$variant.info.neandertal != "match") & (filtered_mutG4s$variant.info.denisovan == "match"),
                        "denisovan", "both"))
)

uniq_mutG4_long <- uniq_mutG4 %>% mutate(id = row_number()) %>%
  pivot_longer(cols = c(wild_score, mut_score),
               names_to = "type",
               values_to = "score") %>%
  mutate(score = abs(score),
         introgressed = factor(introgressed, levels = c("neandertal", "denisovan", "both")),
         x_group = paste(introgressed, type, sep = "_"))

uniq_mutG4_long$x_group <- factor(uniq_mutG4_long$x_group, levels = c(
  "neandertal_wild_score", "neandertal_mut_score",
  "denisovan_wild_score", "denisovan_mut_score",
  "both_wild_score", "both_mut_score"
))

group_labels <- data.frame(
  x = c(1.5, 3.5, 5.5),
  label = c("Neandertal only", "Denisovan only", "Both Archaic Matches")
)

uniq_mutG4_points <- uniq_mutG4_long %>%
  group_by(introgressed, type, score) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(x_group = paste(introgressed, type, sep = "_"),
         x_group = factor(x_group, levels = levels(uniq_mutG4_long$x_group)))

fig1 <- ggplot() +geom_line(data = uniq_mutG4_long, aes(x = x_group, y = score, group = id),
                            color = "#8a9097", linewidth = 0.25, alpha = 0.3) +
  geom_point(data = uniq_mutG4_points, aes(x = x_group, y = score, size = count, color = type)) +
  scale_x_discrete(labels = rep(c("Wild-G4s", "Mutant-G4s"), 3)) +
  scale_color_manual(values = c(wild_score = "#cc1b00", mut_score = "#197ebf"), guide = "none") +
  scale_size_continuous(range = c(2, 6), breaks = c(1, 3, 5, 7)) +
  labs(x = NULL, y = expression("Absolute G4Hunter max scores"), size = "Count") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) + geom_text(data = group_labels,
                aes(x = x, y = max(uniq_mutG4_long$score) + 0.05, label = label),
                inherit.aes = FALSE, size = 4, fontface = "bold")

ggsave("G4_maxscore_changes.pdf", fig1, width = 8, height = 5)

ref_counts <- table(filtered_mutG4s$variant.info.ref)
subset_C <- variants[variants$ref == "C", ]
subset_G <- variants[variants$ref == "G", ]
sample_C <- subset_C[sample(length(subset_C), ref_counts["C"])]
sample_G <- subset_G[sample(length(subset_G), ref_counts["G"])]

random_ref_matched <- c(sample_C, sample_G)

fwrite(data.frame(paste0(seqnames(filtered_mutG4s), ":", filtered_mutG4s$variant.info.ranges.start-1, "-", filtered_mutG4s$variant.info.ranges.start)), 
       "sprime_var.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

fwrite(data.frame(paste0(seqnames(random_ref_matched), ":", start(random_ref_matched)-1, "-", end(random_ref_matched))), 
       "rand_var.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# for vep
filtered_mutG4s_df <- data.frame(filtered_mutG4s)[, c("seqnames", "variant.info.ranges.start", "variant.info.ranges.end", "variant.info.ref", "variant.info.alt")]
filtered_mutG4s_df$seqnames <- str_split(filtered_mutG4s_df$seqnames, "chr", simplify = TRUE)[, 2]
filtered_mutG4s_df$variant.info.ref <- paste0(filtered_mutG4s_df$variant.info.ref, "/", filtered_mutG4s_df$variant.info.alt)
filtered_mutG4s_df$variant.info.alt <- 1
fwrite(filtered_mutG4s_df, "./filtered_mutG4s_df.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#https://jamesotto852.github.io/ggdensity/
ggsave("./all_mutG4s.pdf", plotVarImpact(result), width = 10, height = 4)
ggsave("./filtered_mutG4s.pdf", plotVarImpact(filtered_mutG4s), width = 10, height = 4)

save(list=ls(), file = "./sprime_G4_var_hg19.RData")
