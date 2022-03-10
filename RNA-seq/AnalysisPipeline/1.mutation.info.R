#' @description: plot the mutation information

#### Based on company pannel, aggregate detected mutations
T1 <- data.frame(mutation = c("VHL"),
                abundance = c(35.88), 
                Type = c("Frameshift indel"),
                patient = rep("T1", 1),stringsAsFactors = F)
T2 <- data.frame(mutation = c("VHL", "ARID1B"),
                abundance = c(37.78, 27.28), 
                Type = c("Missense", "Frameshift indel"),
                patient = rep("T2", 2), stringsAsFactors = F)
T3 <- data.frame(mutation = c("VHL", "TERT", "PTEN", "BARD1"),
                abundance = c(20.0, 17.59, 1.61, 1.35), 
                Type = c("Frameshift indel", "Missense", "Missense", "Missense"),
                patient = rep("T3", 4), stringsAsFactors = F)
#Because there are many mutations detected in T4, the mutation frequency is retained >2%; at the same time, VEGFA/CCND1 is added due to multiple mutations at multiple sites
T4 <- data.frame(mutation = c("VHL", "BAP1", "NF1", "MTOR", "COL5A3", "VEGFA"),
                abundance = c(54.80, 58.90, 3.80, 69.10, 18.20, 20.7), 
                Type = c("Frameshift indel", "Frameshift indel", "Missense", "Missense", "Missense", "Missense"),
                patient = rep("T4", 6), stringsAsFactors = F)
all.mutation <- list(T1, T2, T3, T4)
all.mutation <- do.call(rbind, all.mutation)

#Data transformation: length to width
library(tidyr)
all.mutation.abundance <- spread(all.mutation[,-c(3)], patient, abundance)
rownames(all.mutation.abundance) <- all.mutation.abundance$mutation
all.mutation.abundance <- all.mutation.abundance[,-1]
all.mutation.Type <- spread(all.mutation[,-c(2)], patient, Type)
rownames(all.mutation.Type) <- all.mutation.Type$mutation
all.mutation.Type <- all.mutation.Type[,-1]

all.mutation.abundance <- all.mutation.abundance[c("VHL", "BAP1", "MTOR", "ARID1B",  "VEGFA", "COL5A3", "TERT", "NF1", "BARD1", "PTEN"),]
all.mutation.Type <- all.mutation.Type[c("VHL", "BAP1", "MTOR", "ARID1B", "VEGFA", "COL5A3", "TERT", "NF1", "BARD1", "PTEN"),]

library(ComplexHeatmap)
library(circlize)
pdf("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scRNA/mutation.pdf")
column_split <- c("T1", "T2", "T3", "T4")
row_split <- factor(rownames(all.mutation.Type), levels = rownames(all.mutation.Type))
Heatmap(all.mutation.abundance, name = "Mutation abundance", column_split = column_split, width = unit(6, "cm"), height = unit(8, "cm"), cluster_rows = F, cluster_columns = F, col = c("white", "red"), na_col = "white", row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(all.mutation.abundance[i, j])) grid.text(sprintf("%.1f", all.mutation.abundance[i, j]), x, y, gp = gpar(fontsize = 8))})
#colors <- c("#CC62B0", "#297D6D", "#58C6C1")
colors = structure(2:4, names = c("Frameshift indel", "Missense", "Nonsense"))
Heatmap(all.mutation.Type, name = "Type", row_split = row_split, column_split = column_split, border_gp = gpar(col = "grey"), width = unit(6, "cm"), height = unit(8, "cm"), cluster_rows = F, cluster_columns = F, col = colors, na_col = "white",row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()