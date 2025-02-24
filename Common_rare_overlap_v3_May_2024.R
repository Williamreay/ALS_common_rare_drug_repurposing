##########################################

## MND gene-based results - ranking across different input methods

## Common and rare overlap

## V3 - May 2024

## Dr William Reay (2024)

##########################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(corrplot)
library(forcats)
library(gprofiler2)

setwd("~/Desktop/Zac_MND_project/Ranking_of_gene_based_results/")

Rare_rank <- fread("Rare_var_ranks_final/Rare_var_raw_median_and_scaled_median_rank.txt")
Rare_rank <- rename(Rare_rank, "Disruptive (MAF < 0.01)"="Rank_dis_0.01",
                    "Disruptive (MAF < 0.005)"="Rank_dis_0.005", "Damaging missense (MAF < 0.01)"="Rank_dam_0.01",
                    "Damaging missense (MAF < 0.005)"="Rank_dam_0.005", "Disruptive & Damaging (MAF < 0.01)"="Rank_dis_dam_0.01",
                    "Disruptive & Damaging (MAF < 0.005)"="Rank_dis_dam_0.005", "Disruptive & Damaging & Non-damaging missense (MAF < 0.01)"="Rank_dis_dam_mis_0.01",
                    "Disruptive & Damaging & Non-damaging missense (MAF < 0.005)"="Rank_dis_dam_mis_0.005",                              
                    "Rare - Median (raw)"="Median_rank_raw", "Rare - Median (scaled)"="Median_rank_scaled")

Common_rank <- fread("Common_var_ranks_final/Common_var_raw_median_and_scaled_median_rank.txt")
Common_rank <- rename(Common_rank, "Gene"="ID")
Common_rank <- rename(Common_rank, "MAGMA"="MAGMA_rank", "mBAT"="mBAT_rank",
                      "SMR Blood (HEIDI weighted)"="HEIDI_SMR_blood_rank", "SMR Brain"="Raw_SMR_brain_rank",
                      "SMR Brain (HEIDI weighted)"="HEIDI_SMR_brain_rank", "TWAS"="Raw_TWAS_rank", 
                      "TWAS (coloc weighted)"="Coloc_TWAS_rank", "SMR Blood"="Raw_SMR_blood_rank",
                      "Common - Median (raw)"="Median_rank_raw", "Common - Median (scaled)"="Median_rank_scaled")

## Merge

Merged_com_rar <- merge(Common_rank, Rare_rank, by="Gene")

## Top 5% of scaled common var rank

quantile(Merged_com_rar$`Common - Median (scaled)`, 0.05)

Top_5_per_cmr <- Merged_com_rar %>% filter(`Common - Median (scaled)` < 180.6)

CM <- gost(Top_5_per_cmr$Gene, sources = c("GO:MF", "GO:BP", "REAC", "KEGG", "WP"))

quantile(Merged_com_rar$`Rare - Median (scaled)`, 0.05)

Top_5_per_rmr <- Merged_com_rar %>% filter(`Common - Median (scaled)` < 157.5625)

RM <- gost(Top_5_per_rmr$Gene, sources = c("GO:MF", "GO:BP", "REAC", "KEGG", "WP"))

## Nothing much here

## Correlation plot

Corr_in <- Merged_com_rar %>% select(`Common - Median (raw)`, `Rare - Median (raw)`,
                                     `Common - Median (scaled)`, `Rare - Median (scaled)`)

## Corr plot

## Common vs rare: raw = 0.02138688, p = 0.005

## Common vs rare: scaled = 0.052, p = 1.10e-11

cor.test(Corr_in$`Common - Median (scaled)`, Corr_in$`Rare - Median (scaled)`, method = "spearman")

MAT <- cor(Corr_in, method="spearman")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

v2_corr <- corrplot(MAT, method="color", col=col(200), order="hclust", type ="full",addCoef.col = "black", number.cex = 0.3, 
                    tl.col="black", tl.srt=45, tl.cex = 0.3,
                    p.mat = p.mat, sig.level = 0.05, insig = "blank",
                    diag=T)

Merged_com_rar$Common_quant <- ntile(Merged_com_rar$`Common - Median (scaled)`, 10)




## Look at lowly ranked genes annotated Tchem or Tclin

TCRD <- fread("../Target_annot_pharmacological/TCRDv6.1.0_ALLexp.csv", header = T)
TCRD <- TCRD %>% select(`HGNC Sym`, TDL)
TCRD$Gene <- TCRD$`HGNC Sym`

TCRD_merge <- merge(TCRD, Merged_com_rar, by = "Gene")

TCRD_merge$TDL <- factor(TCRD_merge$TDL,                                    # Change ordering manually
                         levels = c("Tdark", "Tbio", "Tchem", "Tclin"))


## FWER correction TWAS PEC look for genes that are consistent

Tclin <- TCRD_merge %>% filter(TDL == "Tclin") %>% arrange(`Common - Median (raw)`)


## Arrange common and rare by Tclin

head(Tclin %>% arrange(`Rare - Median (scaled)`))

head(Tclin %>% arrange(`Common - Median (scaled)`))

## Find top 5% of scaled distribution for each

quantile(TCRD_merge$`Common - Median (scaled)`, 0.05)
## 169.9
quantile(TCRD_merge$`Rare - Median (scaled)`, 0.05)
## 158.86

## Just common

TCRD_merge_cmr <- TCRD_merge %>% filter(`Common - Median (scaled)` < 170)

TCRD_merge_rmr <- TCRD_merge %>% filter(`Rare - Median (scaled)` < 159)


TCRD_merge_both <- TCRD_merge %>% filter(`Common - Median (scaled)` < 170 & `Rare - Median (scaled)` < 159)

## Lowly ranked for both

library(viridis)

ggplot(TCRD_merge_both, aes(x=`Common - Median (scaled)`, y=`Rare - Median (scaled)`, colour=TDL)) +
  geom_point() +
  scale_color_viridis(discrete = T, option = "turbo") +
  theme_bw() +
  ggtitle("Median (scaled) - rare vs common") +
  facet_wrap(~TDL) +
  xlab("Common - Median (scaled rank)") +
  ylab("Rare - Median (scaled rank)")


TCRD_merge_3 <- TCRD_merge %>% arrange(`Common - Median (scaled)`)
TCRD_merge_3 <- TCRD_merge_3[1:500, ]

## Modal

ggplot(TCRD_merge_3, aes(x=`Common - Median (scaled)`, y=`Rare - Median (scaled)`, colour=TDL)) +
  geom_point() +
  scale_color_nejm() +
  theme_bw() +
  geom_text(data=subset(TCRD_merge_3, MAGMA_P == 4.2180e-05),
            label="NEK1", colour="black") +
  ggtitle("Median (scaled) - rare vs common")

ggarrange(P2, P2, common.legend = T)
