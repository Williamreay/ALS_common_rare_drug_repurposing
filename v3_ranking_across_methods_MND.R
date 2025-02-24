##########################################

## MND gene-based results - ranking across different input methods

## V3 - May 2024

## Dr William Reay (2024)

## MAGMA P
## mBAT P
## PEC TWAS Z multiplied by COLOC H4
## SMR MetaBrain Z weighted by HEIDI Z, bounded by max level at 1.96 (nominal significance)
## SMR eQTLgen blood Z weighted by HEIDI Z

## Meaian across all (including those with missing), scaled rank by number of non-missing
## Annotate ranks to TCRD TDL levels

#########################################

library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(viridis)
library(ggrepel)
library(httr)
library(gprofiler2)
library(GGally)
library(corrplot)

setwd("~/Desktop/Zac_MND_project/")

## MAGMA
MAGMA <- fread("MAGMA_string_network/MND_EUR_cons_boundaries_commonSNPs.genes.out")
MAGMA <- MAGMA %>% select(GENE, P)
MAGMA <- rename(MAGMA, "ID"="GENE", "MAGMA_P"="P")

## mBAT
mBAT <- fread("Zac_2024_als_genes/als_mbat.csv")
mBAT <- mBAT %>% select(gene_name, P_mBATcombo)
mBAT <- rename(mBAT, "ID"="gene_name", "mBAT_P"="P_mBATcombo")

## TWAS PEC - Z * COLOC H4 and TWAS Z
TWAS_PEC <- fread("Zac_2024_als_genes/als_twas_pec.csv")
TWAS_PEC$Weighted_Z <- TWAS_PEC$TWAS.Z*TWAS_PEC$COLOC.PP4

## Output in GSEA preranked form
## Correlate Raw TWAS Z with coloc weighted
ggplot(TWAS_PEC, aes(x=TWAS.Z, y=Weighted_Z)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("TWAS Z (raw)") + ylab("TWAS Z (coloc weighted)")

Common_GSEApreranked <- TWAS_PEC %>% select(ID, Weighted_Z)
Common_GSEApreranked <- Common_GSEApreranked %>% filter(!duplicated(ID))
Common_GSEApreranked <- rename(Common_GSEApreranked, "RNK"="Weighted_Z")
Common_GSEApreranked <- Common_GSEApreranked %>% arrange(-RNK)
write.table(Common_GSEApreranked, file="Ranking_of_gene_based_results/Common_var_ranks_final/TWAS_PEC_coloc_weighted_all_genes.rnk",
            sep = "\t", row.names = F, quote = F)

TWAS_PEC$TWAS_Abs_Weighted_Z <- abs(TWAS_PEC$Weighted_Z)
TWAS_PEC$Abs_TWAS_Z <- abs(TWAS_PEC$TWAS.Z)
TWAS_PEC <- TWAS_PEC %>% select(ID, TWAS_Abs_Weighted_Z, Abs_TWAS_Z)

## SMR MetaBrain - Z/scaled HEIDI Z, with max of 1.96, raw Z (SMR)
SMR_Brain <- fread("Zac_2024_als_genes/als_smr_metabrain.csv")
SMR_Brain$SMR_Z <- SMR_Brain$b_SMR/SMR_Brain$se_SMR
SMR_Brain <- SMR_Brain %>% filter(!is.na(p_HEIDI))
SMR_Brain$Abs_SMR_Z <- abs(SMR_Brain$SMR_Z)
SMR_Brain$HEIDI_Z <- abs(qnorm((SMR_Brain$p_HEIDI/2), F))
SMR_Brain$HEIDI_Z <- ifelse(SMR_Brain$HEIDI_Z < 1.959964, 1, SMR_Brain$HEIDI_Z)
SMR_Brain$SMR_Z_HEIDI_weighted <- SMR_Brain$SMR_Z/SMR_Brain$HEIDI_Z

## Output in GSEA preranked form
## Correlate Raw TWAS Z with coloc weighted
ggplot(SMR_Brain, aes(x=SMR_Z, y=SMR_Z_HEIDI_weighted)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("SMR Z (raw)") + ylab("SMR Z (HEIDI weighted)")

Common_GSEApreranked_SMR <- SMR_Brain %>% select(gene_name, SMR_Z_HEIDI_weighted)
Common_GSEApreranked_SMR <- Common_GSEApreranked_SMR %>% filter(!duplicated(gene_name))
Common_GSEApreranked_SMR <- rename(Common_GSEApreranked_SMR, "RNK"="SMR_Z_HEIDI_weighted", "ID"="gene_name")
Common_GSEApreranked_SMR <- Common_GSEApreranked_SMR %>% arrange(-RNK)
write.table(Common_GSEApreranked_SMR, file="Ranking_of_gene_based_results/Common_var_ranks_final/SMR_brain_HEIDI_weighted_all_genes.rnk",
            sep = "\t", row.names = F, quote = F)


SMR_Brain$Abs_SMR_Z_HEIDI_weighted <- SMR_Brain$Abs_SMR_Z/SMR_Brain$HEIDI_Z
SMR_Brain <- rename(SMR_Brain, "ID"="gene_name")
SMR_Brain <- SMR_Brain %>% select(ID, Abs_SMR_Z, Abs_SMR_Z_HEIDI_weighted)

## SMR eQTLgen blood - Z/scaled HEIDI Z, with max of 1.96, and raw Z (SMR)
SMR_Blood <- fread("Zac_2024_als_genes/als_smr_eqtlgen.csv")
SMR_Blood$SMR_Z <- SMR_Blood$b_SMR/SMR_Blood$se_SMR
SMR_Blood <- SMR_Blood %>% filter(!is.na(p_HEIDI))
SMR_Blood$Abs_SMR_Z_blood <- abs(SMR_Blood$SMR_Z)
SMR_Blood$HEIDI_Z <- abs(qnorm((SMR_Blood$p_HEIDI/2), F))
SMR_Blood$HEIDI_Z <- ifelse(SMR_Blood$HEIDI_Z < 1.959964, 1, SMR_Blood$HEIDI_Z)

SMR_Blood$SMR_Z_HEIDI_weighted <- SMR_Blood$SMR_Z/SMR_Blood$HEIDI_Z

## Output in GSEA preranked form
## Correlate Raw TWAS Z with coloc weighted
ggplot(SMR_Blood, aes(x=SMR_Z, y=SMR_Z_HEIDI_weighted)) + geom_point() + geom_smooth() +
  theme_bw() + xlab("SMR Z (raw)") + ylab("SMR Z (HEIDI weighted)")

Common_GSEApreranked_SMR_Blood <- SMR_Blood %>% select(gene_name, SMR_Z_HEIDI_weighted)
Common_GSEApreranked_SMR_Blood <- Common_GSEApreranked_SMR_Blood %>% filter(!duplicated(gene_name))
Common_GSEApreranked_SMR_Blood <- rename(Common_GSEApreranked_SMR_Blood, "RNK"="SMR_Z_HEIDI_weighted", "ID"="gene_name")
Common_GSEApreranked_SMR_Blood <- Common_GSEApreranked_SMR_Blood %>% arrange(-RNK)
write.table(Common_GSEApreranked_SMR_Blood, file="Ranking_of_gene_based_results/Common_var_ranks_final/SMR_Blood_Blood_HEIDI_weighted_all_genes.rnk",
            sep = "\t", row.names = F, quote = F)

SMR_Blood$Abs_SMR_Z_HEIDI_weighted_blood <- SMR_Blood$Abs_SMR_Z/SMR_Blood$HEIDI_Z
SMR_Blood <- rename(SMR_Blood, "ID"="gene_name")
SMR_Blood <- SMR_Blood %>% select(ID, Abs_SMR_Z_blood, Abs_SMR_Z_HEIDI_weighted_blood)

## Merge all the above - keep all of the above

Ranking_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all=T),
                         list(MAGMA, mBAT, TWAS_PEC, SMR_Blood, SMR_Brain))

Ranking_merged <- Ranking_merged %>% filter(!duplicated(ID))

## Rank
Ranking_merged <- Ranking_merged %>% mutate(MAGMA_rank = row_number(MAGMA_P))
Ranking_merged <- Ranking_merged %>% mutate(mBAT_rank = row_number(mBAT_P))
Ranking_merged <- Ranking_merged %>% mutate(Raw_TWAS_rank = row_number(-Abs_TWAS_Z))
Ranking_merged <- Ranking_merged %>% mutate(Coloc_TWAS_rank = row_number(-TWAS_Abs_Weighted_Z))
Ranking_merged <- Ranking_merged %>% mutate(Raw_SMR_blood_rank = row_number(-Abs_SMR_Z_blood))
Ranking_merged <- Ranking_merged %>% mutate(HEIDI_SMR_blood_rank = row_number(-Abs_SMR_Z_HEIDI_weighted_blood))
Ranking_merged <- Ranking_merged %>% mutate(Raw_SMR_brain_rank = row_number(-Abs_SMR_Z))
Ranking_merged <- Ranking_merged %>% mutate(HEIDI_SMR_brain_rank = row_number(-Abs_SMR_Z_HEIDI_weighted))


## Correlation between 8 different ranks

Named_ranking_merged <- Ranking_merged
Named_ranking_merged <- rename(Named_ranking_merged, "MAGMA"="MAGMA_rank", "mBAT"="mBAT_rank",
                               "SMR Blood (HEIDI weighted)"="HEIDI_SMR_blood_rank", "SMR Brain"="Raw_SMR_brain_rank",
                               "SMR Brain (HEIDI weighted)"="HEIDI_SMR_brain_rank", "TWAS"="Raw_TWAS_rank", 
                               "TWAS (coloc weighted)"="Coloc_TWAS_rank", "SMR Blood"="Raw_SMR_blood_rank")

MAT <- cor(Named_ranking_merged[,10:17], method = "spearman", use = "pairwise.complete.obs")

## Load function to calculate corr p-values on matrix

load("Cor_mat_pval.rda")

p.mat <- cor.mtest(Named_ranking_merged[,10:17])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

v2_corr <- corrplot(MAT, method="color", col=col(200), order="hclust", type ="upper",addCoef.col = "black", number.cex = 0.7, 
         tl.col="black", tl.srt=45, tl.cex = 0.7,
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         diag=F)



## Median raw rank

Ranking_merged$Median_rank_raw <- apply(Ranking_merged[,c(10:17)], 1, median, na.rm = T)

## Count number of non-missing ranks per gene

Ranking_merged$Number_non_missing_rnks <- apply(Ranking_merged[,c(10:17)], 1, function(x){ length(which(!is.na(x))) })

hist(Ranking_merged$Number_non_missing_rnks, xlab = "Non missing ranks", main = "Histogram of gene-wise non-missing ranks")

## Scale median rank by number of non-missing ranks

Ranking_merged$Median_rank_scaled <- (Ranking_merged$Median_rank_raw/Ranking_merged$Number_non_missing_rnks)

## sample estimates: rho 0.832221 

cbPalette <- c("#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#CC79A7", "#D2B48C", "#8B4513", "#800080", "#4682B4", "#FFD700")

All <- ggplot(Ranking_merged, aes(x=Median_rank_raw, y=Median_rank_scaled, colour=as.factor(Number_non_missing_rnks))) +
  geom_point(size = 1) +
  theme_bw() + xlab("Median rank (raw)") + ylab("Median rank (scaled by # of non-missing ranks)") +
  scale_color_manual(values = cbPalette) +
  labs(colour="# of non-missing ranks") +
  theme(legend.position = "bottom") +
  ggtitle("All genes")

Subset <- Ranking_merged %>% filter(Median_rank_raw < 50)

S <- ggplot(Subset, aes(x=Median_rank_raw, y=Median_rank_scaled, fill=as.factor(Number_non_missing_rnks))) +
  geom_point(size = 1) +
  theme_bw() + xlab("Median rank (raw)") + ylab("Median rank (scaled by # of non-missing ranks)") +
  scale_color_manual(values = cbPalette) +
  labs(colour="# of non-missing ranks") +
  theme(legend.position = "bottom") +
  ggtitle("Top 50 genes (raw)")

P1 <- ggarrange(All, S, common.legend = T, legend = "right")

## Plot per predictor result of top 10 genes - median scaled

TE <- Named_ranking_merged[,c(1,10:18,20)]
TE <- TE %>% arrange(Median_rank_scaled)
TE <- rename(TE, "Median (raw)"="Median_rank_raw", "Median (scaled)"="Median_rank_scaled")

Top_10_dis <- TE[1:10, ]

Melted_top_10 <-  melt(Top_10_dis, id = c("ID"))
Melted_top_10$Rank <- Melted_top_10$variable

ME <- ggplot(Melted_top_10, aes(x=ID, y=Rank, size=value)) +
  geom_point(colour="black") +
  scale_size(name = "Rank",
             breaks = c(10, 50, 100, 1000, 5000, 10000, 15000)) +
  #scale_color_viridis(discrete = T, option = "A") +
  guides(colour="none") +
  theme_bw() +
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Top scaled median ranks")

OU <- ggarrange(P1, ME, nrow = 2)



ggsave(OU, device = "svg", height = 10, width = 12,
       file="Ranking_of_gene_based_results/NEW_Common_var_ranks_final.svg")


write.table(Ranking_merged, file="Ranking_of_gene_based_results/Common_var_ranks_final/Common_var_raw_median_and_scaled_median_rank.txt",
            sep = "\t", row.names = F, quote = F)
