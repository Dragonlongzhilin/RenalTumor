#' @description: Based on the article, build a new ligand-receptor interaction model by yourself
# Require comma separated interaction

ligand_receptor1 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand-receptor.LiepingChen.2013.NatureReviewImmunology.txt", header = T, sep = "\t", stringsAsFactors = F)
#168
idx1 <- "LiepingChen.2013.NatureReviewImmunology" 
ligand_receptor2 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand-receptor.Ramilowski.2015.NatureComm.literatureSupport.txt", header = T, sep = "\t", stringsAsFactors = F)
#1894
idx2 <- "Ramilowski.2015.NatureComm" 
ligand_receptor3 <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/human_lr_pair.XinShao.2020.BriefinginBioinfomatics.txt", header = T, sep = "\t", stringsAsFactors = F)
#3398
idx3 <- "XinShao.2020.BriefinginBioinfomatics" 

#Required Format
#id_cp_interaction partner_a partner_b protein_name_a protein_name_b annotation_strategy source
#The missing can be left blank, it is NULL
#Because partner_a is required to be ligand and receptor_b is to acceptor, it is directly limited

ligand_receptor.data1 <- data.frame(id_cp_interaction = paste0(idx1, "_", 1:nrow(ligand_receptor1)),
                                    partner_a = ligand_receptor1[,1],
                                    partner_b = ligand_receptor1[,2],
                                    protein_name_a = rep("", nrow(ligand_receptor1)),
                                    protein_name_b = rep("", nrow(ligand_receptor1)),
                                    annotation_strategy = rep("curated", nrow(ligand_receptor1)),
                                    source = rep(idx1, nrow(ligand_receptor1)), stringsAsFactors = F)

ligand_receptor.data2 <- data.frame(id_cp_interaction = paste0(idx2, "_", 1:nrow(ligand_receptor2)),
                                    partner_a = ligand_receptor2[,2],
                                    partner_b = ligand_receptor2[,3],
                                    protein_name_a = rep("", nrow(ligand_receptor2)),
                                    protein_name_b = rep("", nrow(ligand_receptor2)),
                                    annotation_strategy = ligand_receptor2[,5],
                                    source = rep(idx2, nrow(ligand_receptor2)), stringsAsFactors = F)

ligand_receptor.data3 <- data.frame(id_cp_interaction = paste0(idx3, "_", 1:nrow(ligand_receptor3)),
                                    partner_a = ligand_receptor3[,2],
                                    partner_b = ligand_receptor3[,3],
                                    protein_name_a = rep("", nrow(ligand_receptor3)),
                                    protein_name_b = rep("", nrow(ligand_receptor3)),
                                    annotation_strategy = gsub(",", ";", ligand_receptor3$evidence),
                                    source = rep(idx3, nrow(ligand_receptor3)), stringsAsFactors = F)
ligand_receptor.data <- rbind(ligand_receptor.data1, ligand_receptor.data2) 
ligand_receptor.data <- rbind(ligand_receptor.data, ligand_receptor.data3) #5460

library(dplyr)
newdata <- ligand_receptor.data %>% distinct(partner_a, partner_b, .keep_all = TRUE) #3534
saveRDS(newdata, file = "/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.rds")
# Need to be converted to Uniprot id
# Use online tools to convert https://www.uniprot.org/uploadlists/
gene.id <- data.frame(gene_name = unique(c(newdata[,2], newdata[,3])))
write.csv(gene.id, file = "/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/gene.id.csv", quote = F, row.names = F)
write.csv(newdata, file = "/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.csv", quote = F,row.names = F)

#### After manual verification
# Different data resources have more redundancy, and there are many gene aliases that cause duplication
# And cellphoneDB usually accepts the identity of uniport id, so it is converted
ligand_receptor.data <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.csv", header = T, stringsAsFactors = F, sep = ",")
gene.info <- read.table("/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/gene.info.txt", header = F, sep = "\t", stringsAsFactors = F)

ligand_receptor.data.new <- apply(ligand_receptor.data, 1, function(x){
    index1 <- which(gene.info[,1] == x[2])
    index2 <- which(gene.info[,1] == x[3])
    if(length(index1)==1 & length(index2)==1){
        x[2] <- gene.info[index1, 2]
        x[3] <- gene.info[index2, 2]
    }else{
        x[2] <- NA
        x[3] <- NA
    }
    return(x)
})
ligand_receptor.data.new <- as.data.frame(t(ligand_receptor.data.new))
ligand_receptor.data.new <- ligand_receptor.data.new %>% distinct(partner_a, partner_b, .keep_all = TRUE) #3472
write.csv(ligand_receptor.data.new, file = "/data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.new.csv", quote = F,row.names = F)

## Build a custom interaction_input
cd /data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/cellphoneDB
cellphonedb database generate --user-interactions /data/ExtraDisk/sdd/longzhilin/Data/signatureGeneSet/Immune/ligand_receptor/ligand_receptor.data.new.csv