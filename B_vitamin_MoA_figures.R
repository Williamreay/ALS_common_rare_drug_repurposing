######################################

## B vitamin MoA
## Dr William Reay (2024)

######################################

setwd("~/Desktop/Zac_MND_project/Ranking_of_gene_based_results/GSEApreranked/Zif_et_al_ipsmn_data/")

library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
## Read in rare variant ranks

Rmr <- fread("../../Rare_var_ranks_final/Rare_var_raw_median_and_scaled_median_rank.txt")

## Read in B vitamin MoA

B_MoA <- fread("B_vitamin_MoA.txt")

## Merge

B_merge <- merge(B_MoA, Rmr, by = "Gene")

## Dot plot of ranks per MoA - 

R_B <- ggplot(B_merge, aes(x=B_vitamin, y=Median_rank_scaled)) +
  geom_jitter(aes(fill=B_vitamin), 
             colour="black",pch=21, size=2.5, alpha = 0.7, width = 0.15) + 
  geom_hline(yintercept = 100, lty = "dashed", colour="firebrick2") +
  geom_hline(yintercept = 679, lty = "dashed") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() + xlab(" ") +  theme(legend.position = "none") + ylab("Rare scaled median rank") +
  geom_text_repel(data=subset(B_merge, Gene=="DGAT2"),
            label="DGAT2", colour="black") +
  geom_text_repel(data=subset(B_merge, Gene=="HCAR3"),
                  label="HCAR3", colour="black")


## Volcano plot of MoA rare variant results for panel a

MoA_GSEA <- read_excel("../Rare_MoA_GSEA.xlsx")

Panel_A <- ggplot(MoA_GSEA, aes(y=-log10(`FWER p-val`), x=NES)) +
  geom_point() +
  ylim(0, 1.8) +
  geom_text(data=subset(MoA_GSEA, Term=="vitamin b"),
                  label="Vitamin B", colour="dodgerblue3", nudge_x = 0.3) +
  theme_bw() +
  geom_hline(yintercept = 1.3, lty="dashed") +
  ylab("-log10 (FWER)") +
  geom_vline(xintercept = 0, lty="dashed")

ggarrange(Panel_A, R_B, ncol = 2)


  
       