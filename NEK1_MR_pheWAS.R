##############################

## NEK1 IVs - MR-pheWAS

## FinnGen r10

## Dr William Reay (2024)

##############################

library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(TwoSampleMR)

setwd("~/Desktop/Zac_MND_project/NEK1_MR_pheWAS/")

## Common variant eQTL proxy (brain) - T allele corr with increased expression. rs10520157 (G:T) 
## beta (brain - GTEx cerbellum, highest PIP) =  0.14, SE = 0.024, PIP = 0.86 (finemapping derived)


## FinnGen MR-pheWAS updated function for release 10 - same allele config as eQTL

FinnGen_MR_PheWAS <- function(SNP, IV_beta, df) {
  
  df <- df %>% filter(beta != "NA")
  
  df$Z <- sign(df$beta) * sqrt(qchisq(df$pval, 1, lower =F))
  df$SE <- abs(df$beta/df$Z)
  
  ## MR - Wald ratio
  
  df$MR_beta <- df$beta/IV_beta
  df$MR_SE <- abs(df$SE/IV_beta)
  df$MR_pval <- pnorm(abs(df$MR_beta)/df$MR_SE, lower.tail=FALSE) * 2
  
  ## FWER and FDR correction
  
  #df$FWER <- p.adjust(df$MR_pval, method="bonferroni")
  #df$FDR <- p.adjust(df$MR_pval, method="fdr")
  
  return(df)
}


## Read in FinnGen release 10 outcome data - Bonferroni for number of tests 2.076412e-05



Out_1 <- fread("4_169612585_G_T_phenotype_associations.tsv")
Out_1$SNP <- "rs10520157"
Out_1$ref <- "G"
Out_1$alt <- "T"

Out_1 <- Out_1 %>% filter(n_case > 1000)

NEK1_eQTL_MR <- FinnGen_MR_PheWAS("rs10520157", 0.14, Out_1)

## Flip sign to be in the direction of per SD reduction in NEK1

NEK1_eQTL_MR$MR_beta <- (NEK1_eQTL_MR$MR_beta)*-1

write.table(NEK1_eQTL_MR, file="NEK1_eQTL_MR_pheWAS_direction_decreased_NEK1.txt", sep = "\t", row.names = F, quote = F)

## Volcano plot of results

NEK1_eQTL_MR$category <- as.factor(NEK1_eQTL_MR$category)

levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Certain infectious and parasitic diseases (AB1_)"] <- "Infectious disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the digestive system (K11_)"] <- "Diseases of digestive system"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Neurological endpoints"] <-"Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Mental and behavioural disorders (F5_)"] <-"Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Alcohol related diseases"] <-"Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Asthma and related endpoints"] <-"Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of Asthma"] <-"Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of Neurological endpoints"] <-"Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases marked as autimmune origin"] <-"Inflammatory/Automimmune related"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Endocrine, nutritional and metabolic diseases (E4_)"] <-"Metabolic/Endocrine disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Miscellaneous, not yet classified endpoints"] <-"Miscellaneous"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Neoplasms, from cancer register (ICD-O-3)"] <-"Neoplasms"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of COPD"] <-"Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Neoplasms from hospital discharges (CD2_)"] <-"Neoplasms"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Gastrointestinal endpoints"] <- "Diseases of digestive system"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Cardiometabolic endpoints"] <-"Metabolic/Endocrine disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="COPD and related endpoints"] <- "Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Interstitial lung disease endpoints"] <-"Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of Diabetes"] <- "Metabolic/Endocrine disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)"] <-"Inflammatory/Automimmune related"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Common endpoint"] <-"Miscellaneous"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Dental endpoints"] <- "Dental"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diabetes endpoints"] <- "Metabolic/Endocrine disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the eye and adnexa (H7_)"] <- "Opthalamological"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the nervous system (G6_)"] <- "Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the circulatory system (I9_)"] <-"Cardiovascular Disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Rheuma endpoints"] <- "Rheumatological endpoints"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the ear and mastoid process (H8_)"] <-"Ear disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the respiratory system (J10_)"] <-"Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the musculoskeletal system and connective tissue (M13_)"] <- "Musculoskeletal/Connective Tissue disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Psychiatric endpoints from Katri Räikkönen"] <- "Neurological and Psychiatric"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the skin and subcutaneous tissue (L12_)"] <-"Skin/subcutaneous disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Other, not yet classified endpoints (same as #MISC)"] <-"Miscellaneous"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Diseases of the genitourinary system (N14_)"] <-"Genitourinary disease"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of Gastrointestinal endpoints"] <-"Diseases of digestive system"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Pregnancy, childbirth and the puerperium (O15_)"] <-"Pregnancy, childbirth and the puerperium"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Congenital malformations, deformations and chromosomal abnormalities (Q17)"] <- "Congential Malformations"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)"] <- "Miscellaneous"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Drug purchase endpoints"] <-"Drug purchase endpoints"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Comorbidities of Interstitial lung disease endpoints"] <- "Respiratory disorders"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Injury, poisoning and certain other consequences of external causes (ST19_)"] <-"Injury/poisoning"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Codes for special purposes (U22_)"] <-"Miscellaneous"
levels(NEK1_eQTL_MR$category)[levels(NEK1_eQTL_MR$category)=="Factors influencing health status and contact with health services (Z21_)"] <-"Miscellaneous"


NEK1_eQTL_MR <- rename(NEK1_eQTL_MR, "Category"="category")

library(viridis)

ggplot(data = NEK1_eQTL_MR,
           aes(x=MR_beta, y = -log10(MR_pval),
               colour = Category)) +
  geom_point(alpha = 0.8, size = 1.75) +
  geom_hline(yintercept = -log10(4.14e-5), lty = "dashed") +
  labs(x = expression("Log odds per SD decrease in brain NEK1 expression"), y = expression(paste("-log"[10], " P-value"))) +
  theme_bw() +
  theme(legend.position ="right", legend.text = element_text(size=12)) +
  guides(colour=guide_legend(ncol=1)) +
  geom_vline(xintercept = 0, lty = "dashed") +
  ggtitle("NEK1 downregulation pheWAS (nominally significant traits)") +
  geom_text(data=subset(NEK1_eQTL_MR, phenotype=="Soft tissue disorders"),
            label = "Soft tissue disorder, MR-pvalue = 0.0005", colour="black", nudge_x = 0.5, nudge_y = 0.1) +
  geom_text(data=subset(NEK1_eQTL_MR, phenotype=="Pain in limb"),
            label = "Pain in limb, MR-pvalue = 0.001", colour="black", nudge_x = 0.35, nudge_y = 0.1) +
  geom_text(data=subset(NEK1_eQTL_MR, phenotype=="Polycythaemia vera"),
            label = "Polycythaemia vera, MR-pvalue = 0.003", colour="black", nudge_x = 0.7, nudge_y=0.1) +
  scale_colour_viridis(discrete = T, option = "B")
  
  
