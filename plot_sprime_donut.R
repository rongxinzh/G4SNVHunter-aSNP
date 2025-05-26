library(data.table)
library(dplyr)
library(ggpie)
library(ggplot2)

data <- fread("rQ9m3566QAN1KGJl.txt", sep = "\t", header = TRUE) %>% data.frame()

target_types <- c(
  "intergenic_variant", "intron_variant", "upstream_gene_variant",
  "downstream_gene_variant", "regulatory_region_variant"
)
data$group <- ifelse(data$Consequence %in% target_types, data$Consequence, "others")
  
fig1 <- ggdonut(data = data, group_key = "group", count_type = "full",
                label_info = "ratio", label_type = "horizon",
                label_size = 5, label_pos = "in", label_threshold = 10)

ggsave("./sprime_anno.pdf", fig1, width = 9, height = 5.5)