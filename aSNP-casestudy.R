
library(G4SNVHunter)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(stringr)
library(gridExtra)

seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,
              names = paste0("chr", 1:22))
G4_gr <- G4HunterDetect(seq)

aSNP <- fread("./aSNP_Aut.txt", sep = "\t", header = TRUE) %>%
  data.frame()

aSNP_gr <- GRanges(seqnames = aSNP$Chromosome,
                   ranges = IRanges(start = aSNP$Position, end = aSNP$Position),
                   ref = aSNP$Reference..hg38.,
                   alt = aSNP$Alternative,
                   ancestry = aSNP$Haplotype.ancestry)

# G4 formation
asnp_eff <- SNVImpactG4(G4_gr, aSNP_gr, alt_col = "alt")

high_conf <- filterSNVImpact(asnp_eff,
                             raw_score_threshold = 1.5,
                             mut_score_threshold = 1.2)

low_conf <- filterSNVImpact(asnp_eff,
                            mut_score_threshold = 1.2,
                            score_diff_threshold = -0.25)

low_conf <- subsetByOverlaps(low_conf, high_conf, invert = TRUE)

eff_snp_dt <- rbind(data.frame(high_conf, group = "high_confidence"),
                    data.frame(low_conf, group = "low_confidence"))

fwrite(eff_snp_dt, "./aSNP-disruptG4.txt", sep = "\t", col.names = TRUE, 
  row.names = FALSE, quote = FALSE)

# plot changes in G4Hunter scores

plotSNVImpact(high_conf)

eff_snp_dt$ID <- paste0("ID_", 1:nrow(eff_snp_dt)[1])
setDT(eff_snp_dt)
eff_snp_long <- melt(eff_snp_dt[, c("G4.info.score", "mut.score",
                                    "SNV.info.ancestry", "group", "ID")],
                  id.vars = c("SNV.info.ancestry", "group", "ID"),
                  variable.name = "Score_Type", value.name = "Score") %>%
  mutate(Score_Type = recode(Score_Type,
                             "G4.info.score" = "Raw", "mut.score" = "Mut")) %>%
  mutate(group = recode(group,
                             "high_confidence" = "high confidence",
                             "low_confidence" = "low confidence"))

paired_plot <- function(data, title, legend = "none") {
  ggplot(data, aes(y = abs(Score), x = Score_Type,
                   fill = Score_Type, color = group)) +
    geom_line(aes(group = ID), color = "black", linewidth = 0.5, alpha = 0.7) +
    geom_point(aes(shape = Score_Type), size = 3,
               position = position_dodge(width = 0.3), show.legend = TRUE) +
    facet_wrap(~ group, nrow = 1, strip.position = "top") +
    labs(x = NULL, y = NULL, title = title) +
    scale_color_manual(values = c("lightcoral", "darkcyan")) +
    scale_shape_manual(values = c(16, 17)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.title = element_blank(),
      legend.position = legend,
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

p1 <- paired_plot(eff_snp_long %>%
                    filter(SNV.info.ancestry == "Denisovan-like set 1&2"),
                  "Denisovan-like set 1&2")
p2 <- paired_plot(eff_snp_long %>%
                    filter(SNV.info.ancestry == "Denisovan-like set 2"),
                  "Denisovan-like set 2")
p3 <- paired_plot(eff_snp_long %>%
                    filter(SNV.info.ancestry == "Neandertal-like"),
                  "Neandertal-like")

grid.arrange(p1, p2, p3, ncol = 3, widths = c(1, 1, 1))


