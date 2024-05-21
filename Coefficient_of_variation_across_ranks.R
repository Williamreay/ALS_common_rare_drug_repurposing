##########################################

## MND gene-based results - ranking across different input methods

## V3 - May 2024


## Coefficient of variation across ranks
## Common and Rare
## Genes with at least 5 lines of evidence

###########################################

library(EnvStats)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/Zac_MND_project/Ranking_of_gene_based_results/")

## Common ##

Common <- fread("Common_var_ranks_final/Common_var_raw_median_and_scaled_median_rank.txt")
## > 4 ranks
Filt_common <- Common %>% filter(Number_non_missing_rnks > 4)

## 8764 remaining

Filt_common$Row_mean <- apply(Filt_common[,10:17],1,function(x) mean(x, na.rm = TRUE))
Filt_common$Row_sd <- apply(Filt_common[,10:17],1,function(x) sd(x, na.rm = TRUE))
Filt_common$CV <- (Filt_common$Row_sd/Filt_common$Row_mean)*100

P1 <- ggplot(Filt_common, aes(x=CV)) +
  geom_density(fill="steelblue", alpha=0.7) +
  theme_bw() +
  xlab("CV (%)") +
  ylab("Density") +
  ggtitle("Common variant ranks") +
  geom_vline(xintercept = 59, lty="dashed")


## Rare

Rare <- fread("Rare_var_ranks_final/Rare_var_raw_median_and_scaled_median_rank.txt")
## > 4 ranks
Filt_Rare <- Rare %>% filter(Number_non_missing_rnks > 4)

## 9365 remaining

Filt_Rare$Row_mean <- apply(Filt_Rare[,2:9],1,function(x) mean(x, na.rm = TRUE))
Filt_Rare$Row_sd <- apply(Filt_Rare[,2:9],1,function(x) sd(x, na.rm = TRUE))
Filt_Rare$CV <- (Filt_Rare$Row_sd/Filt_Rare$Row_mean)*100

P2 <- ggplot(Filt_Rare, aes(x=CV)) +
  geom_density(fill="firebrick3", alpha=0.7) +
  theme_bw() +
  xlab("CV (%)") +
  ylab("Density") +
  ggtitle("Rare variant ranks") +
  geom_vline(xintercept = 58.5, lty="dashed")

ggarrange(P1, P2)

## Test corr between ranks
Filt_common$Gene <- Filt_common$ID
Merged_CV <- merge(Filt_common, Filt_Rare, by = "Gene")

cor.test(Merged_CV$CV.x, Merged_CV$CV.y)

## No strong correlation between variability in ranks between common and rare

ggplot(Merged_CV, aes(x=CV.x, y=CV.y)) +
  geom_point() +
  geom_smooth(method="lm")

Merged_CV %>% filter(CV.y < 25 & CV.x < 25) %>% select(Median_rank_raw.x, Median_rank_raw.y)
