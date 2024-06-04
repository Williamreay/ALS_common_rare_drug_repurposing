#####################################

## GSEApreranked results - comparing common vs rare scaled

## Dr William Reay (2024)

####################################

library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

setwd("~/Desktop/Zac_MND_project/Ranking_of_gene_based_results/GSEApreranked/Formatting_and_comparing_results/")


## Read in ranks

Common_lv2 <- read_excel("Common_median_scaled_level_2.xlsx")
Common_lv3 <- read_excel("Common_median_scaled_level_3.xlsx")
Common_lv4 <- read_excel("Common_median_scaled_level_4.xlsx")

Rare_lv2 <- read_excel("Rare_median_scaled_level_2.xlsx")
Rare_lv3 <- read_excel("Rare_median_scaled_level_3.xlsx")
Rare_lv4 <- read_excel("Rare_median_scaled_level_4.xlsx")


## Correlate NES at each level using kendel's tau as a small sample

Lvl2 <- merge(Common_lv2, Rare_lv2, by = "Term")
cor.test(Lvl2$NES.x, Lvl2$NES.y, method = "kendall")
## tau = 0.1115, p = 0.1168

Lvl3 <- merge(Common_lv3, Rare_lv3, by = "Term")
cor.test(Lvl3$NES.x, Lvl3$NES.y, method="kendall")
## tau = 0.17, p = 0.00072

Lvl4 <- merge(Common_lv4, Rare_lv4, by = "Term")
cor.test(Lvl4$NES.x, Lvl4$NES.y, method="kendall")
## tau = 0.097, p = 0.0061

## Plot that visualises common vs rare NES - highlight and label drug classes significant after multiple-testing correction
Lvl2$Threshold <- ifelse(Lvl2$Term == "H02", "1", "0")
SP_LVL2 <- ggplot(Lvl2, aes(x=NES.x, y=NES.y, colour=Threshold)) +
  geom_point() +
  scale_color_manual(name="Threshold", values = c("black","firebrick3")) +
  theme_bw() +
  xlab("NES - common variants (scaled)") +
  ylab("NES - rare variants (scaled)") +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_vline(xintercept = 0, lty="dashed") +
  geom_text(data=subset(Lvl2, Term=="H02"),
            label="Corticosteroids (systemic)", nudge_y = 0.15, nudge_x = 0.1, colour="firebrick3") +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  ggtitle("ATC level 2") +
  theme(legend.position = "none")
  #annotate("rect", xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf, fill = "blue", alpha = .07, color = NA) +
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = -1, fill = "orange", alpha = .07, color = NA)
  
SP_LVL3 <- ggplot(Lvl3, aes(x=NES.x, y=NES.y)) +
  geom_point() +
  theme_bw() +
  xlab("NES - common variants (scaled)") +
  ylab("NES - rare variants (scaled)") +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_vline(xintercept = 0, lty="dashed") +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  ggtitle("ATC level 3") +
  theme(legend.position = "none")

Lvl4$Threshold <- ifelse(Lvl4$Term == "L01EC", "1", "0")

SP_LVL4 <- ggplot(Lvl4, aes(x=NES.x, y=NES.y, colour=Threshold)) +
  geom_point() +
  scale_color_manual(name="Threshold", values = c("black","firebrick3")) +
  theme_bw() +
  xlab("NES - common variants (scaled)") +
  ylab("NES - rare variants (scaled)") +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_vline(xintercept = 0, lty="dashed") +
  geom_text(data=subset(Lvl4, Term=="L01EC"),
            label="BRAF inhibitors", colour="firebrick3", nudge_x = 0.3, nudge_y = 0.15) +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  ggtitle("ATC level 4") +
  theme(legend.position = "none")
  

Panel_1 <- ggarrange(SP_LVL2, SP_LVL3, SP_LVL4, nrow = 1)
  

## L01EC boxplot or KDE - include comparing gene-size

L01EC_genes <- fread("L01EC_genes.txt")

## Extract L01EC genes out of all genes

## Read in ranks

Comm_ranks <- fread("../../Data_for_Zac_30_04_24/Common_var_raw_median_and_scaled_median_rank.txt")
Comm_ranks <- Comm_ranks %>% select(ID, Median_rank_scaled)

Plot_input_L01EC <- merge(L01EC_genes, Comm_ranks, by = "ID")
Plot_input_L01EC$Category <- "L01EC"
T <- Comm_ranks[!L01EC_genes, on=.(ID)]
T$Category = "All other genes"
PI <- rbind(Plot_input_L01EC, T)

L01EC_plot <- ggplot(PI, aes(y=log(Median_rank_scaled), x=Category)) +
  geom_boxplot(aes(fill=Category)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab(" ") + ylab("Common variant rank (log scale)") +
  scale_fill_manual(aes(fill=Category), values = c("grey","deepskyblue1")) +
  ggtitle("L01EC: Common variant ranks") +
  geom_text(data=subset(PI, ID=="LIMK2"),
            label="LIMK2", nudge_x = 0.12, nudge_y = -0.2) +
  geom_text_repel(data=subset(PI, ID=="STK36"),
            label="STK36", nudge_x = 0.12, nudge_y = 0.25)
  

## Same as above for H02 with rare

H02_genes <- fread("../../Data_for_Zac_30_04_24/H02_ATC.txt")

Rare_ranks <- fread("../../Data_for_Zac_30_04_24/Rare_var_raw_median_and_scaled_median_rank.txt")
Rare_ranks <- Rare_ranks %>% select(Gene, Median_rank_scaled)
Rare_ranks <- rename(Rare_ranks, "ID"="Gene")

Plot_input_H02 <- merge(H02_genes, Rare_ranks, by = "ID")
Plot_input_H02$Category <- "H02"
T2 <- Rare_ranks[!H02_genes, on=.(ID)]
T2$Category = "All other genes"
PI_2 <- rbind(Plot_input_H02, T2)

H02_plot <- ggplot(PI_2, aes(y=log(Median_rank_scaled), x=Category)) +
  geom_boxplot(aes(fill=Category)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab(" ") + ylab("Rare variant rank (log scale)") +
  scale_fill_manual(aes(fill=Category), values = c("grey","lightsalmon")) +
  ggtitle("H02: Rare variant ranks") +
  geom_text(data=subset(PI_2, ID=="NR3C2"),
            label="NR3C2", nudge_x = 0.12) +
  geom_text_repel(data=subset(PI_2, ID=="CYP11B2"),
            label="CYP11B2", nudge_x = 0.12)

BP <- ggarrange(L01EC_plot, H02_plot)

ggarrange(Panel_1, BP, nrow = 2)
