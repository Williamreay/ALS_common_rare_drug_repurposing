##########################################

## MND gene-based results - ranking across different input methods

## No exclusion of genes with missing rankings - v3 May 2024

## Rare variant gene based

## Dr William Reay (2024)

## Different genic masks and the relationship of their rank to that using synonymous variation

##########################################

setwd("~/Desktop/Zac_MND_project/Ranking_of_gene_based_results/")

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
set.seed(47347)

## Disruptive” variants were those variants classified as frame-shift, splice-site, exon loss, stop gained, start loss and transcription ablation by SnpEff71. 
## “Damaging” variants were missense variants predicted to be damaging by six prediction algorithms (SIFT, Polyphen-2, LRT, MutationTaster2, Mutations Assessor, and PROVEAN). 
## “Missense” variants are those missense variants that did not meet the “damaging” criteria

## Use Firth logistic regression p values

## Combine transcript level P-values into a single P-value per gene using the Cauchy approach
source("~/OneDrive - University of Tasmania/OneDrive_for_UTAS/cloudstor/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Scripts/ACAT_function.R")

## Define function for the following
## 1) Remove NA Firth P-values
## 2) Cauchy combination across transcripts
## 3) Rank by combined P-value (ascending)
## Return Gene name, Combined P-value, and rank

Rare_var_process <- function(input_df, ranking_name) {
  DF <- input_df
  ## Remove NA p-val
  DF <- DF %>% filter(!is.na(p.firth))
  
  ## Cauchy comb across all transcripts per gene ID, genes with only one tested transcript, the same p-val is returned by ACAT function
  genes <- sort(unique(DF$Name))
  gene_ps <- list()
  for (i in genes) { gene_ps[[i]] <- c(DF[DF$Name==i,][,5])}
  genes_acat <- list()
  for (i in genes) { genes_acat[[i]] <- ACAT(unlist(gene_ps[[i]], use.names=FALSE))}
  af_acat <- data.frame(cbind(as.vector(genes), as.vector(unlist(genes_acat))))
  af_acat <- rename(af_acat, c("Gene"="X1", "Cauchy.P"="X2"))
  af_acat$Cauchy.P <- as.numeric(af_acat$Cauchy.P)
  
  ## Ascending rank
  #af_acat <- af_acat %>% mutate(Rank = row_number(Cauchy.P))
  ## Name with correct mask for ranking
  names(af_acat)[names(af_acat) == "Cauchy.P"] <- paste("", ranking_name, sep="")
  #af_acat <- af_acat %>% select(-Cauchy.P)
  return(af_acat)
}

## Set 1 - disruptive, MAF 0.01
Dis_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.disruptive.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_0.01 <- Rare_var_process(Dis_0.01, "Dis_0.01")

## Make GSEA preranked output for disruptive variants

Prernk <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.disruptive.results.proteinCoding.withGeneNames.txt.gz")
Prernk <- Prernk %>% filter(!is.na(p.firth)) %>% arrange(p.firth)
Prernk_filt <- Prernk[!duplicated(Prernk$Name), ]
Prernk_filt$RNK <- Prernk_filt$beta.firth/Prernk_filt$se.firth
Prernk_filt <- Prernk_filt %>% select(Name, RNK)
Prernk_filt <- rename(Prernk_filt, "ID"="Name")
Rare_GSEApreranked_disruptive <- Prernk_filt %>% arrange(-RNK)
write.table(Rare_GSEApreranked_disruptive, file="Rare_var_ranks_final/Rare_variant_disruptive_0.01_MAF.rnk",
            sep = "\t", row.names = F, quote = F)

## Set 2 - disruptive, MAF 0.005
Dis_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.disruptive.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_0.005 <- Rare_var_process(Dis_0.005, "Dis_0.005")

## Set 3 - damaging, MAF 0.01
Dam_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.damaging.results.proteinCoding.withGeneNames.txt.gz")
Rank_dam_0.01 <- Rare_var_process(Dam_0.01, "Dam_0.01")

## Set 4 - damaging, MAF 0.005
Dam_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.damaging.results.proteinCoding.withGeneNames.txt.gz")
Rank_dam_0.005 <- Rare_var_process(Dam_0.005, "Dam_0.005")

## Set 5 - missense, MAF 0.01
Mis_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.missense.results.proteinCoding.withGeneNames.txt.gz")
Rank_mis_0.01 <- Rare_var_process(Mis_0.01, "Mis_0.01")

## Set 6 - missense, MAF 0.005
Mis_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.missense.results.proteinCoding.withGeneNames.txt.gz")
Rank_mis_0.005 <- Rare_var_process(Mis_0.005, "Mis_0.005")

## Set 7 - disruptive & damaging, MAF 0.01
Dis_dam_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.dis.dam.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_dam_0.01 <- Rare_var_process(Dis_dam_0.01, "Dis_dam_0.01")

## Set 8 - disruptive & damaging, MAF 0.005
Dis_dam_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.dis.dam.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_dam_0.005 <- Rare_var_process(Dis_dam_0.005, "Dis_dam_0.005")

## Set 9 - disruptive, damaging, & missense, MAF 0.01
Dis_dam_mis_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.dis.dam.mis.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_dam_mis_0.01 <- Rare_var_process(Dis_dam_mis_0.01, "Dis_dam_mis_0.01")

## Set 10 - disruptive, damaging, & missense MAF 0.005
Dis_dam_mis_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.dis.dam.mis.results.proteinCoding.withGeneNames.txt.gz")
Rank_dis_dam_mis_0.005 <- Rare_var_process(Dis_dam_mis_0.005, "Dis_dam_mis_0.005")

## Set 11 - negative control, synonymous MAF 0.01
Syn_0.01 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.01.synonymous.results.proteinCoding.withGeneNames.txt.gz")
Rank_syn_0.01 <- Rare_var_process(Syn_0.01, "Syn_0.01")

## Set 12 - negative control, synonymous MAF 0.005
Syn_0.005 <- fread("../GWAS_dat/rare_variant_burden/all.step7.exonic.maf0.005.synonymous.results.proteinCoding.withGeneNames.txt.gz")
Rank_syn_0.005 <- Rare_var_process(Syn_0.005, "Syn_0.005")

## Combine all above 12 masks

Ranking_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Gene", all=T),
                         list(Rank_dis_0.01, Rank_dis_0.005,
                              Rank_dam_0.01, Rank_dam_0.005,
                              Rank_mis_0.01, Rank_mis_0.005,
                              Rank_dis_dam_0.01, Rank_dis_dam_0.005,
                              Rank_dis_dam_mis_0.01, Rank_dis_dam_mis_0.005,
                              Rank_syn_0.01, Rank_syn_0.005))

Ranking_merged <- Ranking_merged %>% filter(!duplicated(Gene))

## Apply ranks
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_0.01 = row_number(Dis_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_0.005 = row_number(Dis_0.005))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dam_0.01 = row_number(Dam_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dam_0.005 = row_number(Dam_0.005))
Ranking_merged <- Ranking_merged %>% mutate(Rank_mis_0.01 = row_number(Mis_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_mis_0.005 = row_number(Mis_0.005))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_dam_0.01 = row_number(Dis_dam_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_dam_0.005 = row_number(Dis_dam_0.005))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_dam_mis_0.01 = row_number(Dis_dam_mis_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_dis_dam_mis_0.005 = row_number(Dis_dam_mis_0.005))
Ranking_merged <- Ranking_merged %>% mutate(Rank_syn_0.01 = row_number(Syn_0.01))
Ranking_merged <- Ranking_merged %>% mutate(Rank_syn_0.005 = row_number(Syn_0.005))


## Rename and correlate

Named_ranking_merged <- rename(Ranking_merged, "Disruptive (MAF < 0.01)"="Rank_dis_0.01",
                               "Disruptive (MAF < 0.005)"="Rank_dis_0.005", "Damaging missense (MAF < 0.01)"="Rank_dam_0.01",
                               "Damaging missense (MAF < 0.005)"="Rank_dam_0.005", "Non-damaging missense (MAF < 0.01)"="Rank_mis_0.01",
                               "Non-damaging missense (MAF < 0.005)"="Rank_mis_0.005", "Disruptive & Damaging (MAF < 0.01)"="Rank_dis_dam_0.01",
                               "Disruptive & Damaging (MAF < 0.005)"="Rank_dis_dam_0.005", "Disruptive & Damaging & Non-damaging missense (MAF < 0.01)"="Rank_dis_dam_mis_0.01",
                               "Disruptive & Damaging & Non-damaging missense (MAF < 0.005)"="Rank_dis_dam_mis_0.005",
                               "Synonymous (MAF < 0.01)"="Rank_syn_0.01", "Synonymous (MAF < 0.005)"="Rank_syn_0.005")



## 14:25
library(Hmisc)

MAT <- cor(Named_ranking_merged[,14:25], method = "spearman", use = "pairwise.complete.obs")


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

v2_corr <- corrplot(MAT, method="color", col=col(200), order="hclust", type ="upper",addCoef.col = "black", number.cex = 0.7, 
                    tl.col="black", tl.srt=45, tl.cex = 0.5, 
                    diag=F)


## Combined ranks - removing synonymous and non-damaging missense by themselves

Ranking_merged_v2 <- Ranking_merged[,c(1,14:17,20:23)]

## Median raw rank

Ranking_merged_v2$Median_rank_raw <- apply(Ranking_merged_v2[,c(2:9)], 1, median, na.rm = T)

## Count number of non-missing ranks per gene

Ranking_merged_v2$Number_non_missing_rnks <- apply(Ranking_merged_v2[,c(2:9)], 1, function(x){ length(which(!is.na(x))) })

Ranking_merged_v2 <- Ranking_merged_v2 %>% filter(Number_non_missing_rnks > 0)
hist(Ranking_merged_v2$Number_non_missing_rnks, xlab = "Non missing ranks", main = "Histogram of gene-wise non-missing ranks")

## Scale median rank by number of non-missing ranks

Ranking_merged_v2$Median_rank_scaled <- (Ranking_merged_v2$Median_rank_raw/Ranking_merged_v2$Number_non_missing_rnks)

## sample estimates: rho 0.9541308
cbPalette <- c("#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#CC79A7", "#D2B48C", "#8B4513", "#800080", "#4682B4", "#FFD700")
All <- ggplot(Ranking_merged_v2, aes(x=Median_rank_raw, y=Median_rank_scaled, colour=as.factor(Number_non_missing_rnks))) +
  geom_point(size = 1) +
  theme_bw() + xlab("Median rank (raw)") + ylab("Median rank (scaled by # of non-missing ranks)") +
  scale_color_manual(values = cbPalette) +
  labs(colour="# of non-missing ranks") +
  theme(legend.position = "bottom") +
  ggtitle("All genes")

Subset <- Ranking_merged_v2 %>% filter(Median_rank_raw < 50)

S <- ggplot(Subset, aes(x=Median_rank_raw, y=Median_rank_scaled, colour=as.factor(Number_non_missing_rnks))) +
  geom_point(size = 1) +
  theme_bw() + xlab("Median rank (raw)") + ylab("Median rank (scaled by # of non-missing ranks)") +
  scale_color_manual(values = c("#d95f02", "#e6ab02", "#666666")) +
  labs(colour="# of non-missing ranks") +
  theme(legend.position = "bottom") +
  ggtitle("Top 50 genes (raw)")

P1 <- ggarrange(All, S, common.legend = T, legend = "right")

## Plot per predictor result of top 10 genes (Median < 10)

Named_ranking_merged_v2 <- rename(Ranking_merged_v2, "Disruptive (MAF < 0.01)"="Rank_dis_0.01",
                               "Disruptive (MAF < 0.005)"="Rank_dis_0.005", "Damaging missense (MAF < 0.01)"="Rank_dam_0.01",
                               "Damaging missense (MAF < 0.005)"="Rank_dam_0.005", "Disruptive & Damaging (MAF < 0.01)"="Rank_dis_dam_0.01",
                               "Disruptive & Damaging (MAF < 0.005)"="Rank_dis_dam_0.005", "Disruptive & Damaging & Non-damaging missense (MAF < 0.01)"="Rank_dis_dam_mis_0.01",
                               "Disruptive & Damaging & Non-damaging missense (MAF < 0.005)"="Rank_dis_dam_mis_0.005")

TE <- Named_ranking_merged_v2
TE <- rename(TE, "Median (raw)"="Median_rank_raw", "Median (scaled)"="Median_rank_scaled")
TE <- TE %>% select(-Number_non_missing_rnks)

Top_10_dis <- TE %>% filter(`Median (raw)` < 11)

Melted_top_10 <-  melt(Top_10_dis, id = c("Gene"))
Melted_top_10$Rank <- Melted_top_10$variable

ME_1 <- ggplot(Melted_top_10, aes(x=Gene, y=Rank, size=value)) +
  geom_point(colour="black") +
  scale_size(name = "Rank",
             breaks = c(10, 50, 100, 1000, 5000, 10000, 15000)) +
  scale_color_aaas() +
  guides(colour="none") +
  theme_bw() +
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Median rank (raw) < 10")




write.table(Ranking_merged_v2, file="Rare_var_ranks_final/Rare_var_raw_median_and_scaled_median_rank.txt",
            sep = "\t", row.names = F, quote = F)

## Scaled median < 10

Top_10_sc <- TE %>% filter(`Median (scaled)` < 11)

Melted_top_10 <-  melt(Top_10_sc, id = c("Gene"))
Melted_top_10$Rank <- Melted_top_10$variable

ME <- ggplot(Melted_top_10, aes(x=Gene, y=Rank, size=value)) +
  geom_point(colour="black") +
  scale_size(name = "Rank",
             breaks = c(10, 50, 100, 1000, 5000, 10000, 15000)) +
  guides(colour="none") +
  theme_bw() +
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Median rank (scaled) < 10")

OU <- ggarrange(P1, ME_1, ME, nrow = 3)

ggsave(OU, device = "svg", height = 10, width = 12,
       file="Rare_var_ranks_final/UPDATED_2025_Rare_variant_raw_median_rank_less_than_10.svg")

ggsave(OU, device = "tiff", height = 10, width = 12,
       file="Rare_var_ranks_final/UPDATED_2025_Rare_variant_raw_median_rank_less_than_10.tiff")
