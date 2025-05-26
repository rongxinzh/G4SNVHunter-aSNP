set.seed(1)

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(scales)

mut_G4var_reg <- fread("./sprime_var_regulomescore.txt", sep = "\t", header = FALSE) %>% data.frame()
shuffled_var_reg <- fread("./rand_var_regulomescore.txt", sep = "\t", header = FALSE) %>% data.frame()

group_names <- c("G4-impacting\nintrogressed variants", "Base-matched\nintrogressed variants")

regulome_score_plt <- data.frame(
  group = factor(rep(group_names, c(length(mut_G4var_reg$V5), 
                                    length(shuffled_var_reg$V5))), levels = group_names),
  score = c(mut_G4var_reg$V5, shuffled_var_reg$V5)
)

pval <- wilcox.test(mut_G4var_reg$V5, shuffled_var_reg$V5, alternative = "greater")$p.value
p_label <- ifelse(pval < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(pval, 2)))

fig1 <- ggplot(regulome_score_plt, aes(x = group, y = score, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black") +
  annotate("text", x = 1.5, y = max(regulome_score_plt$score, na.rm = TRUE), 
           label = p_label, size = 4.5, fontface = "italic") +
  scale_fill_manual(values = setNames(c("#D55E00", "#0072B2"), group_names)) +
  labs(x = NULL, y = "RegulomeDB Scores") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("regulome_scores.pdf", fig1, width = 5, height = 4)
