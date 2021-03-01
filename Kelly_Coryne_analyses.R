# Botswana Infant Microbiome Study - 16S Sequencing Data Analyses
# Matthew Kelly, MD, MPH 

remove(list=ls())
setwd("LIST WORKING DIRECTORY") 

set.seed(1234)

library(phyloseq)
library(tidyverse)
library(dplyr)
library(plyr)
library(data.table)
library(gridExtra)
library(metagenomeSeq)
library(httr)
library(reshape2)
library(ggplot2)
library(extrafont)
library(cowplot)
library(DataCombine)
library(vegan)
library(ggpubr)
library(R2admb)
library(glmmADMB)
library(bbmle)
library(performance)

earthy_cols_11 <- c("gray20", "indianred4", "darkslateblue", "thistle", "peru", "mistyrose4", "darkolivegreen4", "mediumorchid4",
                   "coral3", "dodgerblue3", "navajowhite3")
seq_cols_11 <- c("indianred4", "#481567FF", "#453781FF", "#39568CFF", "#2D708EFF", "#238A8DFF", "#20A387FF", "#3CBB75FF", "#73D055FF",
                 "#B8DE29FF", "#FDE725FF")

phy.np.pruned <- read_rds("phy.np.pruned.020121.rds")
nsamples(phy.np.pruned)
ntaxa(phy.np.pruned)
metadata_all <- data.frame(sample_data(phy.np.pruned))
metadata_all$month <- as.factor(metadata_all$month)
metadata_all$SampleType <- as.factor(metadata_all$SampleType)
sample_data(phy.np.pruned) <- metadata_all

#***********************************************************************************************************************
# INFANT CHARACTERISTICS
#***********************************************************************************************************************

metadata_inf <- subset(metadata_all, Source=="INF")
patients <- metadata_inf %>% group_by(study_id) %>% filter(row_number(month) == 1)
infants <- droplevels(patients$study_id)

table(patients$sex)
prop.table(table(patients$sex))
summary(patients$bw)
table(patients$mat_hiv)
prop.table(table(patients$mat_hiv))
mat_hiv <- subset(patients, mat_hiv=="Y")
table(mat_hiv$mat_art_reg, useNA="always")
prop.table(table(mat_hiv$mat_art_reg, useNA="always"))
summary(mat_hiv$mat_art_mo)
summary(mat_hiv$enr_cd4_num)
table(mat_hiv$enr_vl_num, useNA="always")
prop.table(table(mat_hiv$enr_vl_num))
table(patients$residence)
prop.table(table(patients$residence))
table(patients$electric)
prop.table(table(patients$electric))
table(patients$wood)
prop.table(table(patients$wood))
table(patients$mat_educ)
prop.table(table(patients$mat_educ))
summary(patients$num_adults)
summary(patients$num_adol)
summary(patients$num_kids)
table(patients$season)
prop.table(table(patients$season))
remove(mat_hiv)

pcv <- metadata_inf
pcv$pcv[!is.na(pcv$pcv1) | !is.na(pcv$pcv2) | !is.na(pcv$pcv3)] <- 1
pcv$pcv[is.na(pcv$pcv1) & !is.na(pcv$pcv2) & !is.na(pcv$pcv3)] <- 0
pcv$pcv <- as.numeric(pcv$pcv)
pcv_any <- pcv %>% group_by(study_id) %>% tally(pcv)
pcv_any$n <- as.numeric(pcv_any$n)
pcv_any$pcv_any <- NA
pcv_any$pcv_any[pcv_any$n>0] <- "Y"
pcv_any$pcv_any[pcv_any$n==0] <- "N"
table(pcv_any$pcv_any, useNA="always")
prop.table(table(pcv_any$pcv_any))
remove(pcv, pcv_any)

bm <- metadata_inf
bm$bm[bm$breastmilk=="Y"] <- 1
bm$bm[bm$breastmilk=="N"] <- 0
bm$bm <- as.numeric(bm$bm)
bm_any <- bm %>% group_by(study_id) %>% tally(bm)
bm_any$n <- as.numeric(bm_any$n)
bm_any$bm_any <- NA
bm_any$bm_any[bm_any$n>0] <- "Y"
bm_any$bm_any[bm_any$n==0] <- "N"
table(bm_any$bm_any, useNA="always")
prop.table(table(bm_any$bm_any))
remove(bm, bm_any)

abx <- metadata_inf
abx$abx[abx$inf_abx_any=="Y"] <- 1
abx$abx[abx$inf_abx_any=="N"] <- 0
abx$abx <- as.numeric(abx$abx)
abx_any <- abx %>% group_by(study_id) %>% tally(abx)
abx_any$n <- as.numeric(abx_any$n)
abx_any$abx_any <- NA
abx_any$abx_any[abx_any$n>0] <- "Y"
abx_any$abx_any[abx_any$n==0] <- "N"
table(abx_any$abx_any, useNA="always")
prop.table(table(abx_any$abx_any))
abx$amox[abx$inf_abx_amox=="Y"] <- 1
abx$amox[abx$inf_abx_amox=="N"] <- 0
abx$amox <- as.numeric(abx$amox)
amox_any <- abx %>% group_by(study_id) %>% tally(amox)
amox_any$n <- as.numeric(amox_any$n)
amox_any$amox_any <- NA
amox_any$amox_any[amox_any$n>0] <- "Y"
amox_any$amox_any[amox_any$n==0] <- "N"
table(amox_any$amox_any, useNA="always")
prop.table(table(amox_any$amox_any))
abx$metro[abx$inf_abx_metro=="Y"] <- 1
abx$metro[abx$inf_abx_metro=="N"] <- 0
abx$metro <- as.numeric(abx$metro)
metro_any <- abx %>% group_by(study_id) %>% tally(metro)
metro_any$n <- as.numeric(metro_any$n)
metro_any$metro_any <- NA
metro_any$metro_any[metro_any$n>0] <- "Y"
metro_any$metro_any[metro_any$n==0] <- "N"
table(metro_any$metro_any, useNA="always")
prop.table(table(metro_any$metro_any))
abx$cotrim[abx$inf_abx_cotrim=="Y"] <- 1
abx$cotrim[abx$inf_abx_cotrim=="N"] <- 0
abx$cotrim <- as.numeric(abx$cotrim)
cotrim_any <- abx %>% group_by(study_id) %>% tally(cotrim)
cotrim_any$n <- as.numeric(cotrim_any$n)
cotrim_any$cotrim_any <- NA
cotrim_any$cotrim_any[cotrim_any$n>0] <- "Y"
cotrim_any$cotrim_any[cotrim_any$n==0] <- "N"
table(cotrim_any$cotrim_any, useNA="always")
prop.table(table(cotrim_any$cotrim_any))
remove(abx, abx_any, amox_any, metro_any, cotrim_any)

# **********************************************************************************************************************
# CLEANING UP TAXTABLE
# **********************************************************************************************************************

taxtable <- as.data.frame(as(tax_table(phy.np.pruned),"matrix"),stringsAsFactors=FALSE)
taxtable <- taxtable[,-c(7,8)]
table(taxtable$domain)
table(taxtable$phylum)
#table(taxtable$class)
#table(taxtable$order)
#table(taxtable$family)
#table(taxtable$genus)
taxtable$genus[taxtable$genus=="Corynebacterium_1"] <- "Corynebacterium"
taxtable$genus[taxtable$genus=="Prevotella_2"] <- "Prevotella"
taxtable$genus[taxtable$genus=="Prevotella_6"] <- "Prevotella"
taxtable$genus[taxtable$genus=="Prevotella_7"] <- "Prevotella"
taxtable$genus[taxtable$genus=="Prevotella_9"] <- "Prevotella"
taxtable$genus[taxtable$genus=="Ruminiclostridium_1"] <- "Ruminiclostridium"
taxtable$genus[taxtable$genus=="Ruminiclostridium_5"] <- "Ruminiclostridium"
taxtable$genus[taxtable$genus=="Ruminiclostridium_6"] <- "Ruminiclostridium"
taxtable$genus[taxtable$genus=="Ruminiclostridium_9"] <- "Ruminiclostridium"
taxtable$genus[taxtable$genus=="Ruminococcus_1"] <- "Ruminococcus"
taxtable$genus[taxtable$genus=="Ruminococcus_2"] <- "Ruminococcus"
taxtable$genus[taxtable$genus=="Tyzzerella_3"] <- "Tyzzerella"
taxtable$genus[taxtable$genus=="Tyzzerella_4"] <- "Tyzzerella"
taxtable$genus[taxtable$genus=="Treponema_2"] <- "Treponema"
taxtable$genus[taxtable$genus=="Selenomonas_3"] <- "Selenomonas"
taxtable$genus[taxtable$genus=="Lachnoclostridium_5"] <- "Lachnoclostridium"
taxtable$genus[taxtable$genus=="Lachnoclostridium_10"] <- "Lachnoclostridium"
taxtable$genus[taxtable$genus=="Coprococcus_1"] <- "Coprococcus"
taxtable$genus[taxtable$genus=="Coprococcus_2"] <- "Coprococcus"
taxtable$genus[taxtable$genus=="Coprococcus_3"] <- "Coprococcus"

# Add species for Corynebacterium ASVs based on BLAST searches
coryne_blast <- read.csv("blast_coryne_030121.csv")
row.names(coryne_blast) <- coryne_blast$ASV
coryne_blast <- coryne_blast[,-c(3,4)]
taxtable_coryne <- merge(taxtable, coryne_blast, by="row.names", all.x = TRUE)
row.names(taxtable_coryne) <- taxtable_coryne[,1]
taxtable_coryne <- taxtable_coryne[,-1]
remove(taxtable)

# Replace taxtable with modified taxtable
taxtable <- taxtable_coryne
taxtable_matrix <- as.matrix(taxtable)
tax_table(phy.np.pruned) <- taxtable_matrix
remove(taxtable_coryne, taxtable_matrix)

#***********************************************************************************************************************
# SAMPLE CLUSTERING
#***********************************************************************************************************************

ord_theme   <-  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), axis.text.y = element_text(size=12), 
                      axis.title.y = element_text(size=16), axis.text.x  = element_text(size=12, angle=90, hjust=1, vjust=0.4), 
                      axis.title.x = element_text(size=16), strip.text.x=element_text(size=18,angle=0),
                      strip.background=element_rect(fill="black"), legend.position = "right",
                      legend.text = element_text(size=19, family = "mono", color="grey20", face="bold"), 
                      legend.title = element_blank(), legend.key = element_blank())

phy.np.cluster <- phy.np.pruned
metadata_tmp <- data.frame(sample_data(phy.np.cluster))
metadata_tmp$month2[metadata_tmp$Source=="MOM"] <- "M0"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="0"] <- "I0"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="1"] <- "I1"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="2"] <- "I2"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="3"] <- "I3"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="4"] <- "I4"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="5"] <- "I5"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="6"] <- "I6"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="8"] <- "I8"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="10"] <- "I10"
metadata_tmp$month2[metadata_tmp$Source=="INF" & metadata_tmp$month=="12"] <- "I12"
table(metadata_tmp$month2)
metadata_tmp$month2 <- as.factor(metadata_tmp$month2)
metadata_tmp$month2 <- reorder(metadata_tmp$month2, new.order=c("M0", "I0", "I1", "I2", "I3", "I4", "I5", "I6", "I8", "I10", "I12"))
sample_data(phy.np.cluster) <- metadata_tmp

# These graphs look at PCoA by SAMPLE TYPE and AGE using Bray-Curtis distances
ord <- ordinate(phy.np.cluster, method="PCoA", distance="bray")
lines <- c("solid", "solid", "blank", "blank", "blank", "blank", "blank", "blank", "blank", "blank", "solid")
bray <- plot_ordination(phy.np.cluster, ord, color="month2", shape="Source") +
  geom_point(size=2) + ord_theme + scale_color_manual(values=seq_cols_11) +
  stat_ellipse(aes(linetype=month2), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) + 
  xlab("PC1 (22.5%)") + ylab("PC2 (12.3%)") +
  scale_linetype_manual(values=lines)
png(file="R Plots/Figure_2a.png", 
    width = 9.2, height = 8.3, units = 'in', res=600)
bray
dev.off()
remove(bray, ord, ord_theme, metadata_tmp)

# Comparison of overall microbiome composition by PERMANOVA
phy.I0_vs_I1 <- subset_samples(phy.np.cluster, month2=="I0" | month2=="I1")
metadata_I0_vs_I1 <- data.frame(sample_data(phy.I0_vs_I1))
adonis(distance(phy.I0_vs_I1, method="bray") ~ month2, strata=metadata_I0_vs_I1$study_id, data = metadata_I0_vs_I1)
phy.I1_vs_I2 <- subset_samples(phy.np.cluster, month2=="I1" | month2=="I2")
metadata_I1_vs_I2 <- data.frame(sample_data(phy.I1_vs_I2))
adonis(distance(phy.I1_vs_I2, method="bray") ~ month2, strata=metadata_I1_vs_I2$study_id, data = metadata_I1_vs_I2)
phy.I2_vs_I3 <- subset_samples(phy.np.cluster, month2=="I2" | month2=="I3")
metadata_I2_vs_I3 <- data.frame(sample_data(phy.I2_vs_I3))
adonis(distance(phy.I2_vs_I3, method="bray") ~ month2, strata=metadata_I2_vs_I3$study_id, data = metadata_I2_vs_I3)
phy.I3_vs_I4 <- subset_samples(phy.np.cluster, month2=="I3" | month2=="I4")
metadata_I3_vs_I4 <- data.frame(sample_data(phy.I3_vs_I4))
adonis(distance(phy.I3_vs_I4, method="bray") ~ month2, strata=metadata_I3_vs_I4$study_id, data = metadata_I3_vs_I4)
phy.I4_vs_I5 <- subset_samples(phy.np.cluster, month2=="I4" | month2=="I5")
metadata_I4_vs_I5 <- data.frame(sample_data(phy.I4_vs_I5))
adonis(distance(phy.I4_vs_I5, method="bray") ~ month2, strata=metadata_I4_vs_I5$study_id, data = metadata_I4_vs_I5)
phy.I5_vs_I6 <- subset_samples(phy.np.cluster, month2=="I5" | month2=="I6")
metadata_I5_vs_I6 <- data.frame(sample_data(phy.I5_vs_I6))
adonis(distance(phy.I5_vs_I6, method="bray") ~ month2, strata=metadata_I5_vs_I6$study_id, data = metadata_I5_vs_I6)
phy.I6_vs_I8 <- subset_samples(phy.np.cluster, month2=="I6" | month2=="I8")
metadata_I6_vs_I8 <- data.frame(sample_data(phy.I6_vs_I8))
adonis(distance(phy.I6_vs_I8, method="bray") ~ month2, strata=metadata_I6_vs_I8$study_id, data = metadata_I6_vs_I8)
phy.I8_vs_I10 <- subset_samples(phy.np.cluster, month2=="I8" | month2=="I10")
metadata_I8_vs_I10 <- data.frame(sample_data(phy.I8_vs_I10))
adonis(distance(phy.I8_vs_I10, method="bray") ~ month2, strata=metadata_I8_vs_I10$study_id, data = metadata_I8_vs_I10)
phy.I10_vs_I12 <- subset_samples(phy.np.cluster, month2=="I10" | month2=="I12")
metadata_I10_vs_I12 <- data.frame(sample_data(phy.I10_vs_I12))
adonis(distance(phy.I10_vs_I12, method="bray") ~ month2, strata=metadata_I10_vs_I12$study_id, data = metadata_I10_vs_I12)

phy.M0_vs_I0 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I0")
metadata_M0_vs_I0 <- data.frame(sample_data(phy.M0_vs_I0))
adonis(distance(phy.M0_vs_I0, method="bray") ~ Source, data = metadata_M0_vs_I0) 
phy.M0_vs_I1 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I1")
metadata_M0_vs_I1 <- data.frame(sample_data(phy.M0_vs_I1))
adonis(distance(phy.M0_vs_I1, method="bray") ~ Source, data = metadata_M0_vs_I1)
phy.M0_vs_I2 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I2")
metadata_M0_vs_I2 <- data.frame(sample_data(phy.M0_vs_I2))
adonis(distance(phy.M0_vs_I2, method="bray") ~ Source, data = metadata_M0_vs_I2)
phy.M0_vs_I3 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I3")
metadata_M0_vs_I3 <- data.frame(sample_data(phy.M0_vs_I3))
adonis(distance(phy.M0_vs_I3, method="bray") ~ Source, data = metadata_M0_vs_I3)
phy.M0_vs_I4 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I4")
metadata_M0_vs_I4 <- data.frame(sample_data(phy.M0_vs_I4))
adonis(distance(phy.M0_vs_I4, method="bray") ~ Source, data = metadata_M0_vs_I4)
phy.M0_vs_I5 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I5")
metadata_M0_vs_I5 <- data.frame(sample_data(phy.M0_vs_I5))
adonis(distance(phy.M0_vs_I5, method="bray") ~ Source, data = metadata_M0_vs_I5)
phy.M0_vs_I6 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I6")
metadata_M0_vs_I6 <- data.frame(sample_data(phy.M0_vs_I6))
adonis(distance(phy.M0_vs_I6, method="bray") ~ Source, data = metadata_M0_vs_I6)
phy.M0_vs_I8 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I8")
metadata_M0_vs_I8 <- data.frame(sample_data(phy.M0_vs_I8))
adonis(distance(phy.M0_vs_I8, method="bray") ~ Source, data = metadata_M0_vs_I8)
phy.M0_vs_I10 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I10")
metadata_M0_vs_I10 <- data.frame(sample_data(phy.M0_vs_I10))
adonis(distance(phy.M0_vs_I10, method="bray") ~ Source, data = metadata_M0_vs_I10)
phy.M0_vs_I12 <- subset_samples(phy.np.cluster, month2=="M0" | month2=="I12")
metadata_M0_vs_I12 <- data.frame(sample_data(phy.M0_vs_I12))
adonis(distance(phy.M0_vs_I12, method="bray") ~ Source, data = metadata_M0_vs_I12)

# Comparison of overall microbiome composition by PERMANOVA in mothers by HIV status
phy.M0 <- subset_samples(phy.np.cluster, month2=="M0")
metadata_M0 <- data.frame(sample_data(phy.M0))
adonis(distance(phy.M0, method="bray") ~ mat_hiv, data = metadata_M0)

remove(metadata_I0_vs_I1, metadata_I1_vs_I2, metadata_I2_vs_I3, metadata_I3_vs_I4, metadata_I4_vs_I5, metadata_I5_vs_I6, metadata_I6_vs_I8,
       metadata_I8_vs_I10, metadata_I10_vs_I12, metadata_M0_vs_I0, metadata_M0_vs_I1, metadata_M0_vs_I2, metadata_M0_vs_I3, metadata_M0_vs_I4,
       metadata_M0_vs_I5, metadata_M0_vs_I6, metadata_M0_vs_I8, metadata_M0_vs_I10, metadata_M0_vs_I12, metadata_M0)
remove(phy.I0_vs_I1, phy.I1_vs_I2, phy.I2_vs_I3, phy.I3_vs_I4, phy.I4_vs_I5, phy.I5_vs_I6, phy.I6_vs_I8, phy.I8_vs_I10, phy.I10_vs_I12,
       phy.M0_vs_I0, phy.M0_vs_I1, phy.M0_vs_I2, phy.M0_vs_I3, phy.M0_vs_I4, phy.M0_vs_I5, phy.M0_vs_I6, phy.M0_vs_I8, phy.M0_vs_I10,
       phy.M0_vs_I12, phy.M0, phy.np.cluster)

# **********************************************************************************************************************
# ALPHA DIVERSITY
# **********************************************************************************************************************

# Create dataframe with Shannon diversity index for each specimen
diversity_tmp <- estimate_richness(phy.np.pruned, measures = c("Shannon", "Chao1"))
setDT(diversity_tmp, keep.rownames = TRUE)[]
diversity_tmp$SampleID <- diversity_tmp$rn
diversity_tmp$rn <- NULL
diversity_tmp$SampleID <- gsub("[.]", "-", diversity_tmp$SampleID)
diversity_all <- merge(metadata_all, diversity_tmp, by="SampleID")
remove(diversity_tmp, metadata_all)

nsamples(phy.np.pruned)
sum(sample_sums(phy.np.pruned))
mean(sample_sums(phy.np.pruned))
summary(diversity_all$Shannon)
summary(diversity_all$Chao1)

phy.np.pruned.inf <- subset_samples(phy.np.pruned, Source=="INF")
nsamples(phy.np.pruned.inf)
sum(sample_sums(phy.np.pruned.inf))
mean(sample_sums(phy.np.pruned.inf))
diversity_inf <- subset(diversity_all, Source=="INF")
summary(diversity_inf$Shannon)
summary(diversity_inf$Chao1)

diversity_plot <- diversity_all
diversity_plot$month2[diversity_plot$Source=="MOM" & diversity_plot$month=="0"] <- "M0"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="0"] <- "I0"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="1"] <- "I1"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="2"] <- "I2"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="3"] <- "I3"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="4"] <- "I4"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="5"] <- "I5"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="6"] <- "I6"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="8"] <- "I8"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="10"] <- "I10"
diversity_plot$month2[diversity_plot$Source=="INF" & diversity_plot$month=="12"] <- "I12"
diversity_plot$month2 <- as.factor(diversity_plot$month2)
diversity_plot$month2 <- reorder(diversity_plot$month2, new.order=c("M0", "I0", "I1", "I2", "I3", "I4", "I5", "I6", "I8", "I10", "I12"))

sdi_np_age <- ggplot(diversity_plot, aes(x=month2, y=Shannon, color=Source)) + 
  geom_boxplot(aes(fill=Source, group=month2), color="black", 
               fill=c("indianred1", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue",
                      "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue")) +
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=15, family="mono", color="grey20", face="bold"), 
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=15, family="mono", color="grey20", face="bold"),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=16), legend.position="none") + ylab("Shannon diversity index") + xlab("Time (months)")
chao_np_age <- ggplot(diversity_plot, aes(x=month2, y=Chao1, color=Source)) + 
  geom_boxplot(aes(fill=Source, group=month2), color="black", 
               fill=c("indianred1", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue",
                      "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue")) +
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=15, family="mono", color="grey20", face="bold"), 
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=15, family="mono", color="grey20", face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background=element_rect(fill="white"), 
        strip.text.y=element_text(size=16), legend.position="none") + ylab("Chao1 index") + xlab("Time (months)") + ylim(0,400)

png(file="R Plots/Figure_1.png", 
    width = 12, height = 7, units = 'in', res = 300)
plot_grid(sdi_np_age, chao_np_age, labels="auto", label_size=22) 
dev.off()
remove(chao_np_age, sdi_np_age)

# Compare alpha diversity in maternal and infant samples

diversity_M0_I0 <- subset(diversity_plot, month2=="M0" | month2=="I0")
diversity_M0_I1 <- subset(diversity_plot, month2=="M0" | month2=="I1")
diversity_M0_I2 <- subset(diversity_plot, month2=="M0" | month2=="I2")
diversity_M0_I3 <- subset(diversity_plot, month2=="M0" | month2=="I3")
diversity_M0_I4 <- subset(diversity_plot, month2=="M0" | month2=="I4")
diversity_M0_I5 <- subset(diversity_plot, month2=="M0" | month2=="I5")
diversity_M0_I6 <- subset(diversity_plot, month2=="M0" | month2=="I6")
diversity_M0_I8 <- subset(diversity_plot, month2=="M0" | month2=="I8")
diversity_M0_I10 <- subset(diversity_plot, month2=="M0" | month2=="I10")
diversity_M0_I12 <- subset(diversity_plot, month2=="M0" | month2=="I12")

# Chao1 richness lower than maternal samples through 3 months of age
wilcox.test(diversity_M0_I0$Chao1 ~ diversity_M0_I0$Source)
wilcox.test(diversity_M0_I1$Chao1 ~ diversity_M0_I1$Source)
wilcox.test(diversity_M0_I2$Chao1 ~ diversity_M0_I2$Source)
wilcox.test(diversity_M0_I3$Chao1 ~ diversity_M0_I3$Source)
wilcox.test(diversity_M0_I4$Chao1 ~ diversity_M0_I4$Source)
wilcox.test(diversity_M0_I5$Chao1 ~ diversity_M0_I5$Source)
wilcox.test(diversity_M0_I6$Chao1 ~ diversity_M0_I6$Source)
wilcox.test(diversity_M0_I8$Chao1 ~ diversity_M0_I8$Source)
wilcox.test(diversity_M0_I10$Chao1 ~ diversity_M0_I10$Source)
wilcox.test(diversity_M0_I12$Chao1 ~ diversity_M0_I12$Source)

remove(diversity_M0_I0, diversity_M0_I1, diversity_M0_I2, diversity_M0_I3, diversity_M0_I4, diversity_M0_I5, diversity_M0_I6,
       diversity_M0_I8, diversity_M0_I10, diversity_M0_I12)

# Chao1 richness increases with age
ggdensity(diversity_inf$Chao1)
Chao1_model <- diversity_inf
Chao1_model$month <- as.numeric(Chao1_model$month)
Chao1_age <- glmmadmb(Chao1 ~ month (1|study_id), data=Chao1_model, family="nbinom")
summary(Chao1_age)
remove(Chao1_age, Chao1_model)

# Compare alpha diversity in birth samples compared with later timepoints

diversity_I0_I1 <- subset(diversity_plot, month2=="I0" | month2=="I1")
diversity_I0_I2 <- subset(diversity_plot, month2=="I0" | month2=="I2")
diversity_I0_I3 <- subset(diversity_plot, month2=="I0" | month2=="I3")
diversity_I0_I4 <- subset(diversity_plot, month2=="I0" | month2=="I4")
diversity_I0_I5 <- subset(diversity_plot, month2=="I0" | month2=="I5")
diversity_I0_I6 <- subset(diversity_plot, month2=="I0" | month2=="I6")
diversity_I0_I8 <- subset(diversity_plot, month2=="I0" | month2=="I8")
diversity_I0_I10 <- subset(diversity_plot, month2=="I0" | month2=="I10")
diversity_I0_I12 <- subset(diversity_plot, month2=="I0" | month2=="I12")

# Shannon diversity at birth differs from diversity at all other infant timepoints
wilcox.test(diversity_I0_I1$Shannon ~ diversity_I0_I1$month2)
wilcox.test(diversity_I0_I2$Shannon ~ diversity_I0_I2$month2)
wilcox.test(diversity_I0_I3$Shannon ~ diversity_I0_I3$month2)
wilcox.test(diversity_I0_I4$Shannon ~ diversity_I0_I4$month2)
wilcox.test(diversity_I0_I5$Shannon ~ diversity_I0_I5$month2)
wilcox.test(diversity_I0_I6$Shannon ~ diversity_I0_I6$month2)
wilcox.test(diversity_I0_I8$Shannon ~ diversity_I0_I8$month2)
wilcox.test(diversity_I0_I10$Shannon ~ diversity_I0_I10$month2)
wilcox.test(diversity_I0_I12$Shannon ~ diversity_I0_I12$month2)

remove(diversity_plot, diversity_I0_I1, diversity_I0_I2, diversity_I0_I3, diversity_I0_I4, diversity_I0_I5, diversity_I0_I6,
       diversity_I0_I8, diversity_I0_I10, diversity_I0_I12)

# **********************************************************************************************************************
# MICROBIOME COMPOSITION
# **********************************************************************************************************************

# ******************************
# NORMALIZATION IN METAGENOMESEQ
# ******************************

phy.np.normalized <- phy.np.pruned
mgs <- phyloseq_to_metagenomeSeq(phy.np.normalized)
phenoData(mgs)
featureData(mgs)

# Normalize to the cumulative 75%ile
mgs.norm = cumNorm(mgs, p = 0.75)
normFactor = normFactors(mgs.norm)
normFactor = log2(normFactor/median(normFactor) + 1)
normFactor <- normFactor[!is.na(names(normFactor))]

# Output ASV table with normalized counts to a dataframe and import this back into phyloseq as otu_table
asv1 <- MRcounts(mgs.norm, norm = TRUE, log = FALSE)
asv1 <- otu_table(asv1, taxa_are_rows=TRUE)
otu_table(phy.np.normalized) <- asv1
remove(mgs, mgs.norm, asv1, normFactor)

# *****************
# COMPOSITION PLOTS
# *****************

# Agglomerate taxa into genera for compositional comparisons
ntaxa(phy.np.normalized)
phy.np.agglom <- tax_glom(phy.np.normalized, taxrank = 'genus')
ntaxa(phy.np.agglom)

# Transform to relative abundances
phy.np.relative <- transform_sample_counts(phy.np.agglom, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.np.relative))  
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling
relative_np <- psmelt(phy.np.relative)

# Create dataframes with overall relative abundances of phyla and genera
relative_np$phylum <- as.character(relative_np$phylum)
phyla_abundances <- aggregate(relative_np$Abundance, by=list(phylum=relative_np$phylum), FUN=sum)
phyla_abundances$x <- (phyla_abundances$x)/(nsamples(phy.np.relative))
phyla_abundances <- rename(phyla_abundances, c("x"="phyla_Ab"))
sum(phyla_abundances$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances)            # Corresponds to # of unique phyla
relative_np$genus <- as.character(relative_np$genus)
genus_abundances <- aggregate(relative_np$Abundance, by=list(phylum=relative_np$phylum, genus=relative_np$genus,
                                                                 OTU=relative_np$OTU), FUN=mean)
genus_abundances <- aggregate(genus_abundances$x, by=list(phylum=genus_abundances$phylum,
                                                          genus=genus_abundances$genus), FUN=sum)
genus_abundances <- rename(genus_abundances, c("x"="genus_Ab"))
sum(genus_abundances$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genus_abundances)            # Corresponds to # of unique genera
abundances <- merge(genus_abundances, phyla_abundances, by="phylum")
genus_abundances <- arrange(genus_abundances, genus_Ab, decreasing=TRUE)  
TOPGenera <- unique(genus_abundances$genus[1:10])
genus_df <- genus_abundances[genus_abundances$genus %in% TOPGenera,]
genus_df$genus <- factor(genus_df$genus, levels = genus_df$genus[order(-genus_df$genus_Ab)])
phyla_abundances <- arrange(phyla_abundances, phyla_Ab, decreasing=TRUE)  
head(phyla_abundances, 5)
head(genus_abundances, 19)

# Rename genera other than top genera as "Other" in creating dataframe relative_np 
domination_tmp <- relative_np[,c("SampleID","genus","Abundance")]
relative_np$genus[!(relative_np$genus %in% TOPGenera)] <- "Other"
sum(relative_np$Abundance)
remove(abundances, genus_df)

relative_np$month2[relative_np$Source=="MOM"] <- "M0"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="0"] <- "I0"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="1"] <- "I1"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="2"] <- "I2"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="3"] <- "I3"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="4"] <- "I4"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="5"] <- "I5"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="6"] <- "I6"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="8"] <- "I8"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="10"] <- "I10"
relative_np$month2[relative_np$Source=="INF" & relative_np$month=="12"] <- "I12"
table(relative_np$month2)
relative_np$month2 <- as.factor(relative_np$month2)
relative_np$month2 <- reorder(relative_np$month2, new.order=c("M0", "I0", "I1", "I2", "I3", "I4", "I5", "I6", "I8", "I10", "I12"))

relative_np$genus <- as.factor(relative_np$genus)
relative_np$genus <- reorder(relative_np$genus, new.order=c("Acinetobacter", "Corynebacterium", "Dolosigranulum", "Gemella", "Haemophilus", 
                                                            "Lactobacillus", "Moraxella", "Staphylococcus", "Streptococcus", "Vibrio", "Other"))

plot_relative_np <- ggplot(arrange(relative_np, genus), aes(x=month2, y=Abundance, fill=genus)) +
  geom_bar(stat="identity", position="fill") +
  theme(panel.background = element_blank(), legend.text=element_text(size=16, face="italic"), legend.title=element_text(size=18), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=18, family="mono", color="grey20", face="bold"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=18, family="mono", color="grey20", face="bold"), plot.title = element_text(size=17, hjust=1.4), 
        strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) + 
  scale_fill_manual(values=alpha(earthy_cols_11, c(1,1,1,1,1,1,1,1,1,1,1)), 
                    labels=c("Acinetobacter", "Corynebacterium", "Dolosigranulum", "Gemella", "Haemophilus", 
                            "Lactobacillus", "Moraxella", "Staphylococcus", "Streptococcus", "Vibrio", "Other")) + 
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Genera")) +
  xlab("Time (months)") + ylab("Relative abundance") 
# Change alphas in scale_fill_manual to alter the transparency of specific colors

png(file="R Plots/Figure_2b.png", 
    width = 13, height = 8, units = 'in', res = 600)
print(plot_relative_np)
dev.off()
remove(plot_relative_np)

# Examine infant nasopharyngeal microbiome composition at the birth visit (I0)

phy.relative.I0_vs_I1 <- subset_samples(phy.np.relative, month=="0" & Source=="INF")
relative_np_I0 <- psmelt(phy.relative.I0_vs_I1)
relative_np_I0$phylum <- as.character(relative_np_I0$phylum)
phyla_abundances_I0 <- aggregate(relative_np_I0$Abundance, by=list(phylum=relative_np_I0$phylum), FUN=sum)
phyla_abundances_I0$x <- (phyla_abundances_I0$x)/(nsamples(phy.relative.I0_vs_I1))
phyla_abundances_I0 <- rename(phyla_abundances_I0, c("x"="phyla_Ab"))
sum(phyla_abundances_I0$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances_I0)            # Corresponds to # of unique phyla
relative_np_I0$genus <- as.character(relative_np_I0$genus)
genus_abundances_I0 <- aggregate(relative_np_I0$Abundance, by=list(phylum=relative_np_I0$phylum, genus=relative_np_I0$genus,
                                                             OTU=relative_np_I0$OTU), FUN=mean)
genus_abundances_I0 <- aggregate(genus_abundances_I0$x, by=list(phylum=genus_abundances_I0$phylum,
                                                          genus=genus_abundances_I0$genus), FUN=sum)
genus_abundances_I0 <- rename(genus_abundances_I0, c("x"="genus_Ab"))
sum(genus_abundances_I0$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genus_abundances_I0)            # Corresponds to # of unique genera
abundances_I0 <- merge(genus_abundances_I0, phyla_abundances_I0, by="phylum")
genus_abundances_I0 <- arrange(genus_abundances_I0, genus_Ab, decreasing=TRUE)  
phyla_abundances_I0 <- arrange(phyla_abundances_I0, phyla_Ab, decreasing=TRUE)  
head(phyla_abundances_I0, 5)
head(genus_abundances_I0, 20)

# Examine infant nasopharyngeal microbiome composition at the birth visit (I0)

phy.relative.I1_I12 <- subset_samples(phy.np.relative, month!="0" & Source=="INF")
relative_np_I1_I12 <- psmelt(phy.relative.I1_I12)
relative_np_I1_I12$phylum <- as.character(relative_np_I1_I12$phylum)
phyla_abundances_I1_I12 <- aggregate(relative_np_I1_I12$Abundance, by=list(phylum=relative_np_I1_I12$phylum), FUN=sum)
phyla_abundances_I1_I12$x <- (phyla_abundances_I1_I12$x)/(nsamples(phy.relative.I1_I12))
phyla_abundances_I1_I12 <- rename(phyla_abundances_I1_I12, c("x"="phyla_Ab"))
sum(phyla_abundances_I1_I12$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances_I1_I12)            # Corresponds to # of unique phyla
relative_np_I1_I12$genus <- as.character(relative_np_I1_I12$genus)
genus_abundances_I1_I12 <- aggregate(relative_np_I1_I12$Abundance, by=list(phylum=relative_np_I1_I12$phylum, genus=relative_np_I1_I12$genus,
                                                                   OTU=relative_np_I1_I12$OTU), FUN=mean)
genus_abundances_I1_I12 <- aggregate(genus_abundances_I1_I12$x, by=list(phylum=genus_abundances_I1_I12$phylum,
                                                                genus=genus_abundances_I1_I12$genus), FUN=sum)
genus_abundances_I1_I12 <- rename(genus_abundances_I1_I12, c("x"="genus_Ab"))
sum(genus_abundances_I1_I12$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genus_abundances_I1_I12)            # Corresponds to # of unique genera
abundances_I1_I12 <- merge(genus_abundances_I1_I12, phyla_abundances_I1_I12, by="phylum")
genus_abundances_I1_I12 <- arrange(genus_abundances_I1_I12, genus_Ab, decreasing=TRUE)  
phyla_abundances_I1_I12 <- arrange(phyla_abundances_I1_I12, phyla_Ab, decreasing=TRUE)  
head(phyla_abundances_I1_I12, 5)
head(genus_abundances_I1_I12, 20)

remove(phy.np.agglom, phy.np.relative, relative_np, phyla_abundances_I0, genus_abundances_I0, phyla_abundances_I1_I12, genus_abundances_I1_I12,
       abundances_I0, abundances_I1_I12, relative_np_I0, relative_np_I1_I12, phy.relative.I0_vs_I1, phy.relative.I1_I12)

# **********************************************************************************************************************
# CATEGORIZE NP SAMPLES BY BIOTYPE & CREATE CIRCOS PLOTS
# **********************************************************************************************************************

# Classify infant NP samples based on dominant genus
domination_tmp <- unique(domination_tmp[,c(1:3)])
domination_tmp2 <- subset(domination_tmp, Abundance>=0.50)
domination <- domination_tmp2[,-3]
biotypes <- merge(diversity_inf, domination, by="SampleID", all.x=TRUE)
names(biotypes)[names(biotypes)=="genus"] <- "biotype"
table(biotypes$biotype, useNA="always")
biotypes$biotype[is.na(biotypes$biotype)] <- "Biodiverse"
table(biotypes$biotype, useNA="always")
core_biotypes = c("Biodiverse", "Corynebacterium", "Dolosigranulum", "Haemophilus", "Moraxella", "Staphylococcus", "Streptococcus")
biotypes$biotype[!(biotypes$biotype %in% core_biotypes)] <- "Other"
table(biotypes$biotype, useNA="always")
prop.table(table(biotypes$biotype, useNA="always"))
remove(domination_tmp, domination_tmp2, domination)

# Shows domination of NP microbiota by month
prop.table(table(biotypes$month, biotypes$biotype), 1)

biotypes_legend <- ggplot(arrange(biotypes, biotype), aes(x=month, y=biotype, fill=biotype)) +
  geom_bar(stat="identity", position="fill") +
  theme(panel.background = element_blank(), legend.text=element_text(size=24, face="italic"), legend.title=element_text(size=24), 
        axis.title.y = element_text(angle=90, size=18), 
        axis.text.y = element_text(size=16, colour="black"), axis.title.x = element_text(size=18), 
        axis.text.x = element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust=1.4), 
        strip.text.x = element_text(size = 16), strip.background=element_rect(fill="white")) + 
  scale_fill_manual(values=c("Black", "darkslateblue", "indianred4", "peru", "darkolivegreen4", 
                             "mediumorchid4", "coral3", "navajowhite3"), 
                    labels=c("Biodiverse", "Corynebacterium", "Dolosigranulum", "Haemophilus", 
                             "Moraxella", "Staphylococcus", "Streptococcus", "Other")) + 
  guides(fill = guide_legend(ncol=1, byrow=TRUE, title="Biotype")) +
  xlab("Time (months)") + ylab("Relative abundance") 
png(file="R Plots/Biotypes_legend.png", 
    width = 20, height = 15, units = 'in', res = 600)
print(biotypes_legend)
dev.off()
remove(biotypes_legend)

library("circlize")
colors_circos = c("Corynebacterium" = "darkslateblue", "Dolosigranulum" = "indianred4", "Haemophilus" = "peru", 
                  "Moraxella" = "darkolivegreen4", "Staphylococcus" = "mediumorchid4", "Streptococcus" = "coral3", 
                  "Other" = "navajowhite3", "Biodiverse" = "Black")

# Create chord diagram that contains transitions from month 0 to 1
biotypes_0to1 <- subset(biotypes, month=="0" | month=="1")
biotypes_0to1 <- biotypes_0to1[,c("study_id", "month", "biotype")]
biotypes_0to1 <- spread(biotypes_0to1, month, biotype)
biotypes_0to1 <- subset(biotypes_0to1, !is.na(biotypes_0to1[,2]))
biotypes_0to1 <- subset(biotypes_0to1, !is.na(biotypes_0to1[,3]))
colnames(biotypes_0to1)[2] <- "Month_0"
colnames(biotypes_0to1)[3] <- "Month_1"
circos_0to1 <- as.matrix(table(biotypes_0to1$`Month_0`, biotypes_0to1$`Month_1`))
png(file="R Plots/Chord Diagrams/Month0_to_1.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_0to1, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_0to1, circos_0to1)

# Create chord diagram that contains transitions from month 1 to 2
biotypes_1to2 <- subset(biotypes, month=="1" | month=="2")
biotypes_1to2 <- biotypes_1to2[,c("study_id", "month", "biotype")]
biotypes_1to2 <- spread(biotypes_1to2, month, biotype)
biotypes_1to2 <- subset(biotypes_1to2, !is.na(biotypes_1to2[,2]))
biotypes_1to2 <- subset(biotypes_1to2, !is.na(biotypes_1to2[,3]))
colnames(biotypes_1to2)[2] <- "Month_1"
colnames(biotypes_1to2)[3] <- "Month_2"
circos_1to2 <- as.matrix(table(biotypes_1to2$`Month_1`, biotypes_1to2$`Month_2`))
png(file="R Plots/Chord Diagrams/Month1_to_2.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_1to2, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_1to2, circos_1to2)

# Create chord diagram that contains transitions from month 2 to 3
biotypes_2to3 <- subset(biotypes, month=="2" | month=="3")
biotypes_2to3 <- biotypes_2to3[,c("study_id", "month", "biotype")]
biotypes_2to3 <- spread(biotypes_2to3, month, biotype)
biotypes_2to3 <- subset(biotypes_2to3, !is.na(biotypes_2to3[,2]))
biotypes_2to3 <- subset(biotypes_2to3, !is.na(biotypes_2to3[,3]))
colnames(biotypes_2to3)[2] <- "Month_2"
colnames(biotypes_2to3)[3] <- "Month_3"
circos_2to3 <- as.matrix(table(biotypes_2to3$`Month_2`, biotypes_2to3$`Month_3`))
png(file="R Plots/Chord Diagrams/Month2_to_3.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_2to3, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_2to3, circos_2to3)

# Create chord diagram that contains transitions from month 3 to 4
biotypes_3to4 <- subset(biotypes, month=="3" | month=="4")
biotypes_3to4 <- biotypes_3to4[,c("study_id", "month", "biotype")]
biotypes_3to4 <- spread(biotypes_3to4, month, biotype)
biotypes_3to4 <- subset(biotypes_3to4, !is.na(biotypes_3to4[,2]))
biotypes_3to4 <- subset(biotypes_3to4, !is.na(biotypes_3to4[,3]))
colnames(biotypes_3to4)[2] <- "Month_3"
colnames(biotypes_3to4)[3] <- "Month_4"
circos_3to4 <- as.matrix(table(biotypes_3to4$`Month_3`, biotypes_3to4$`Month_4`))
png(file="R Plots/Chord Diagrams/Month3_to_4.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_3to4, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_3to4, circos_3to4)

# Create chord diagram that contains transitions from month 4 to 5
biotypes_4to5 <- subset(biotypes, month=="4" | month=="5")
biotypes_4to5 <- biotypes_4to5[,c("study_id", "month", "biotype")]
biotypes_4to5 <- spread(biotypes_4to5, month, biotype)
biotypes_4to5 <- subset(biotypes_4to5, !is.na(biotypes_4to5[,2]))
biotypes_4to5 <- subset(biotypes_4to5, !is.na(biotypes_4to5[,3]))
colnames(biotypes_4to5)[2] <- "Month_4"
colnames(biotypes_4to5)[3] <- "Month_5"
circos_4to5 <- as.matrix(table(biotypes_4to5$`Month_4`, biotypes_4to5$`Month_5`))
png(file="R Plots/Chord Diagrams/Month4_to_5.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_4to5, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_4to5, circos_4to5)

# Create matrix that contains transitions from month 5 to 6
biotypes_5to6 <- subset(biotypes, month=="5" | month=="6")
biotypes_5to6 <- biotypes_5to6[,c("study_id", "month", "biotype")]
biotypes_5to6 <- spread(biotypes_5to6, month, biotype)
biotypes_5to6 <- subset(biotypes_5to6, !is.na(biotypes_5to6[,2]))
biotypes_5to6 <- subset(biotypes_5to6, !is.na(biotypes_5to6[,3]))
colnames(biotypes_5to6)[2] <- "Month_5"
colnames(biotypes_5to6)[3] <- "Month_6"
circos_5to6 <- as.matrix(table(biotypes_5to6$`Month_5`, biotypes_5to6$`Month_6`))
png(file="R Plots/Chord Diagrams/Month5_to_6.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_5to6, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_5to6, circos_5to6)

# Create matrix that contains transitions from month 6 to 8
biotypes_6to8 <- subset(biotypes, month=="6" | month=="8")
biotypes_6to8 <- biotypes_6to8[,c("study_id", "month", "biotype")]
biotypes_6to8 <- spread(biotypes_6to8, month, biotype)
biotypes_6to8 <- subset(biotypes_6to8, !is.na(biotypes_6to8[,2]))
biotypes_6to8 <- subset(biotypes_6to8, !is.na(biotypes_6to8[,3]))
colnames(biotypes_6to8)[2] <- "Month_6"
colnames(biotypes_6to8)[3] <- "Month_8"
circos_6to8 <- as.matrix(table(biotypes_6to8$`Month_6`, biotypes_6to8$`Month_8`))
png(file="R Plots/Chord Diagrams/Month6_to_8.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_6to8, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_6to8, circos_6to8)

# Create matrix that contains transitions from month 8 to 10
biotypes_8to10 <- subset(biotypes, month=="8" | month=="10")
biotypes_8to10 <- biotypes_8to10[,c("study_id", "month", "biotype")]
biotypes_8to10 <- spread(biotypes_8to10, month, biotype)
biotypes_8to10 <- subset(biotypes_8to10, !is.na(biotypes_8to10[,2]))
biotypes_8to10 <- subset(biotypes_8to10, !is.na(biotypes_8to10[,3]))
colnames(biotypes_8to10)[2] <- "Month_8"
colnames(biotypes_8to10)[3] <- "Month_10"
circos_8to10 <- as.matrix(table(biotypes_8to10$`Month_8`, biotypes_8to10$`Month_10`))
png(file="R Plots/Chord Diagrams/Month8_to_10.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_8to10, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_8to10, circos_8to10)

# Create matrix that contains transitions from month 10 to 12
biotypes_10to12 <- subset(biotypes, month=="10" | month=="12")
biotypes_10to12 <- biotypes_10to12[,c("study_id", "month", "biotype")]
biotypes_10to12 <- spread(biotypes_10to12, month, biotype)
biotypes_10to12 <- subset(biotypes_10to12, !is.na(biotypes_10to12[,2]))
biotypes_10to12 <- subset(biotypes_10to12, !is.na(biotypes_10to12[,3]))
colnames(biotypes_10to12)[2] <- "Month_10"
colnames(biotypes_10to12)[3] <- "Month_12"
circos_10to12 <- as.matrix(table(biotypes_10to12$`Month_10`, biotypes_10to12$`Month_12`))
png(file="R Plots/Chord Diagrams/Month10_to_12.png", 
    width = 8, height = 8, units = 'in', res = 600)
chordDiagram(circos_10to12, grid.col=colors_circos, directional = 1, direction.type = "arrows",
             link.arr.length = 0.5, annotationTrack="grid")
dev.off()
remove(biotypes_10to12, circos_10to12)

# ********************************************************************************************************************** 
# CREATE FILE FOR BLAST SEARCHES FOR HIGHLY ABUNDANT ASVs
# ********************************************************************************************************************** 

# Create file with top 100 ASVs 
phy.np.normalized.inf <- subset_samples(phy.np.normalized, Source=="INF")
taxa_sums <- data.frame(sort(taxa_sums(phy.np.normalized.inf), TRUE)[1:100])
np_taxa <- merge(taxa_sums, taxtable, by="row.names", all.x=TRUE)
colnames(np_taxa)[2] <- "Abundance"
rownames(np_taxa) <- np_taxa[,1]
np_taxa <- np_taxa[,-1]
np_taxa_seqs <- data.frame(refseq(phy.np.normalized.inf))
top_taxa <- merge(np_taxa, np_taxa_seqs, by="row.names")
write.csv(top_taxa, "blast_top_taxa.csv")
remove(np_taxa_seqs, taxa_sums, np_taxa)

# ********************************************************************************************************************** 
# ADD RELATIVE ABUNDANCES OF CORYNEBACTERIUM SPECIES AND OTHER HIGHLY ABUNDANT GENERA
# ********************************************************************************************************************** 

# Transform to relative abundances
phy.np.inf.relative_asv <- transform_sample_counts(phy.np.normalized.inf, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.np.inf.relative_asv))
# This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling

# Corynebacterium spp.

phy.coryne <- subset_taxa(phy.np.inf.relative_asv, genus=="Corynebacterium")
ntaxa(phy.coryne)
coryne_refseq <- data.frame(refseq(phy.coryne))
coryne_refseq <- tibble::rownames_to_column(coryne_refseq, "ASV")
write.csv(coryne_refseq, "blast_coryne.csv")
remove(coryne_refseq)
coryne_all <- data.frame(otu_table(phy.coryne))
coryne_all <- data.frame(t(coryne_all), check.names = FALSE)
coryne_all$coryne_ALL <- rowSums(coryne_all)
coryne_all$SampleID <- row.names(coryne_all)
coryne_all <- coryne_all[,c("SampleID", "coryne_ALL")]

phy.coryne1 <- subset_taxa(phy.coryne, blast_species=="pseudodiphtheriticum_propinquum")
ntaxa(phy.coryne1)
coryne1 <- data.frame(otu_table(phy.coryne1))
coryne1 <- data.frame(t(coryne1), check.names = FALSE)
coryne1$pseudodiphtheriticum_propinquum <- rowSums(coryne1)
coryne1$SampleID <- row.names(coryne1)
coryne1 <- coryne1[,c("SampleID", "pseudodiphtheriticum_propinquum")]
summary(coryne1$pseudodiphtheriticum_propinquum)
(length(which(coryne1$pseudodiphtheriticum_propinquum!=0))/nrow(coryne1))

phy.coryne2 <- subset_taxa(phy.coryne, blast_species=="accolens_macginleyi")
ntaxa(phy.coryne2)
coryne2 <- data.frame(otu_table(phy.coryne2))
coryne2 <- data.frame(t(coryne2), check.names = FALSE)
coryne2$accolens_macginleyi <- rowSums(coryne2)
coryne2$SampleID <- row.names(coryne2)
coryne2 <- coryne2[,c("SampleID", "accolens_macginleyi")]
summary(coryne2$accolens_macginleyi)
(length(which(coryne2$accolens_macginleyi!=0))/nrow(coryne2))

phy.coryne3 <- subset_taxa(phy.coryne, blast_species=="tuberculostearicum")
ntaxa(phy.coryne3)
coryne3 <- data.frame(otu_table(phy.coryne3))
coryne3 <- data.frame(t(coryne3), check.names = FALSE)
coryne3$tuberculostearicum <- rowSums(coryne3)
coryne3$SampleID <- row.names(coryne3)
coryne3 <- coryne3[,c("SampleID", "tuberculostearicum")]
summary(coryne3$tuberculostearicum)
(length(which(coryne3$tuberculostearicum!=0))/nrow(coryne3))

phy.coryne4 <- subset_taxa(phy.coryne, blast_species=="aurimucosum")
ntaxa(phy.coryne5)
coryne4 <- data.frame(otu_table(phy.coryne4))
coryne4 <- data.frame(t(coryne4), check.names = FALSE)
coryne4$aurimucosum <- rowSums(coryne4)
coryne4$SampleID <- row.names(coryne4)
coryne4 <- coryne4[,c("SampleID", "aurimucosum")]
summary(coryne4$aurimucosum)
(length(which(coryne4$aurimucosum!=0))/nrow(coryne4))

phy.coryne5 <- subset_taxa(phy.coryne, blast_species=="striatum_simulans")
ntaxa(phy.coryne5)
coryne5 <- data.frame(otu_table(phy.coryne5))
coryne5 <- data.frame(t(coryne5), check.names = FALSE)
coryne5$striatum_simulans <- rowSums(coryne5)
coryne5$SampleID <- row.names(coryne5)
coryne5 <- coryne5[,c("SampleID", "striatum_simulans")]
summary(coryne5$striatum_simulans)
(length(which(coryne5$striatum_simulans!=0))/nrow(coryne5))

phy.coryne6 <- subset_taxa(phy.coryne, blast_species=="amycolatum_lactis")
ntaxa(phy.coryne6)
coryne6 <- data.frame(otu_table(phy.coryne6))
coryne6 <- data.frame(t(coryne6), check.names = FALSE)
coryne6$amycolatum_lactis <- rowSums(coryne6)
coryne6$SampleID <- row.names(coryne6)
coryne6 <- coryne6[,c("SampleID", "amycolatum_lactis")]
summary(coryne6$amycolatum_lactis)
(length(which(coryne6$amycolatum_lactis!=0))/nrow(coryne6))

phy.coryne7 <- subset_taxa(phy.coryne, blast_species=="kroppenstedtii")
ntaxa(phy.coryne7)
coryne7 <- data.frame(otu_table(phy.coryne7))
coryne7 <- data.frame(t(coryne7), check.names = FALSE)
coryne7$kroppenstedtii <- rowSums(coryne7)
coryne7$SampleID <- row.names(coryne7)
coryne7 <- coryne7[,c("SampleID", "kroppenstedtii")]
summary(coryne7$kroppenstedtii)
(length(which(coryne7$kroppenstedtii!=0))/nrow(coryne7))

phy.coryne8 <- subset_taxa(phy.coryne, blast_species=="ihumii_pilbarense_mucifaciens")
ntaxa(phy.coryne8)
coryne8 <- data.frame(otu_table(phy.coryne8))
coryne8 <- data.frame(t(coryne8), check.names = FALSE)
coryne8$ihumii_pilbarense_mucifaciens <- rowSums(coryne8)
coryne8$SampleID <- row.names(coryne8)
coryne8 <- coryne8[,c("SampleID", "ihumii_pilbarense_mucifaciens")]
summary(coryne8$ihumii_pilbarense_mucifaciens)
(length(which(coryne8$ihumii_pilbarense_mucifaciens!=0))/nrow(coryne8))

phy.coryne9 <- subset_taxa(phy.coryne, blast_species=="xerosis_freneyi")
ntaxa(phy.coryne9)
coryne9 <- data.frame(otu_table(phy.coryne9))
coryne9 <- data.frame(t(coryne9), check.names = FALSE)
coryne9$xerosis_freneyi <- rowSums(coryne9)
coryne9$SampleID <- row.names(coryne9)
coryne9 <- coryne9[,c("SampleID", "xerosis_freneyi")]
summary(coryne9$xerosis_freneyi)
(length(which(coryne9$xerosis_freneyi!=0))/nrow(coryne9))

phy.coryne10 <- subset_taxa(phy.coryne, blast_species=="minutissimum_spheniscorum_singulare")
ntaxa(phy.coryne10)
coryne10 <- data.frame(otu_table(phy.coryne10))
coryne10 <- data.frame(t(coryne10), check.names = FALSE)
coryne10$minutissimum_spheniscorum_singulare <- rowSums(coryne10)
coryne10$SampleID <- row.names(coryne10)
coryne10 <- coryne10[,c("SampleID", "minutissimum_spheniscorum_singulare")]
summary(coryne10$minutissimum_spheniscorum_singulare)
(length(which(coryne10$minutissimum_spheniscorum_singulare!=0))/nrow(coryne10))

phy.coryne11 <- subset_taxa(phy.coryne, blast_species=="coyleae")
ntaxa(phy.coryne11)
coryne11 <- data.frame(otu_table(phy.coryne11))
coryne11 <- data.frame(t(coryne11), check.names = FALSE)
coryne11$coyleae <- rowSums(coryne11)
coryne11$SampleID <- row.names(coryne11)
coryne11 <- coryne11[,c("SampleID", "coyleae")]
summary(coryne11$coyleae)
(length(which(coryne11$coyleae!=0))/nrow(coryne11))

phy.coryne12 <- subset_taxa(phy.coryne, blast_species=="appendicis_aquatimens")
ntaxa(phy.coryne12)
coryne12 <- data.frame(otu_table(phy.coryne12))
coryne12 <- data.frame(t(coryne12), check.names = FALSE)
coryne12$appendicis_aquatimens <- rowSums(coryne12)
coryne12$SampleID <- row.names(coryne12)
coryne12 <- coryne12[,c("SampleID", "appendicis_aquatimens")]
summary(coryne12$appendicis_aquatimens)
(length(which(coryne12$appendicis_aquatimens!=0))/nrow(coryne12))

phy.coryne13 <- subset_taxa(phy.coryne, blast_species=="jeikeium")
ntaxa(phy.coryne13)
coryne13 <- data.frame(otu_table(phy.coryne13))
coryne13 <- data.frame(t(coryne13), check.names = FALSE)
coryne13$jeikeium <- rowSums(coryne13)
coryne13$SampleID <- row.names(coryne13)
coryne13 <- coryne13[,c("SampleID", "jeikeium")]
summary(coryne13$jeikeium)
(length(which(coryne13$jeikeium!=0))/nrow(coryne13))

phy.coryne14 <- subset_taxa(phy.coryne, blast_species=="lowii_mastitidis")
ntaxa(phy.coryne14)
coryne14 <- data.frame(otu_table(phy.coryne14))
coryne14 <- data.frame(t(coryne14), check.names = FALSE)
coryne14$lowii_mastitidis <- rowSums(coryne14)
coryne14$SampleID <- row.names(coryne14)
coryne14 <- coryne14[,c("SampleID", "lowii_mastitidis")]
summary(coryne14$lowii_mastitidis)
(length(which(coryne14$lowii_mastitidis!=0))/nrow(coryne14))

phy.coryne15 <- subset_taxa(phy.coryne, blast_species=="diphtheriae_vitaeruminis")
ntaxa(phy.coryne15)
coryne15 <- data.frame(otu_table(phy.coryne15))
coryne15 <- data.frame(t(coryne15), check.names = FALSE)
coryne15$diphtheriae_vitaeruminis <- rowSums(coryne15)
coryne15$SampleID <- row.names(coryne15)
coryne15 <- coryne15[,c("SampleID", "diphtheriae_vitaeruminis")]
summary(coryne15$diphtheriae_vitaeruminis)
(length(which(coryne15$diphtheriae_vitaeruminis!=0))/nrow(coryne15))

phy.coryne16 <- subset_taxa(phy.coryne, blast_species=="oculi_ciconiae")
ntaxa(phy.coryne16)
coryne16 <- data.frame(otu_table(phy.coryne16))
coryne16 <- data.frame(t(coryne16), check.names = FALSE)
coryne16$oculi_ciconiae <- rowSums(coryne16)
coryne16$SampleID <- row.names(coryne16)
coryne16 <- coryne16[,c("SampleID", "oculi_ciconiae")]
summary(coryne16$oculi_ciconiae)
(length(which(coryne16$oculi_ciconiae!=0))/nrow(coryne16))

phy.coryne17 <- subset_taxa(phy.coryne, blast_species=="urealyticum")
ntaxa(phy.coryne17)
coryne17 <- data.frame(otu_table(phy.coryne17))
coryne17 <- data.frame(t(coryne17), check.names = FALSE)
coryne17$urealyticum <- rowSums(coryne17)
coryne17$SampleID <- row.names(coryne17)
coryne17 <- coryne17[,c("SampleID", "urealyticum")]
summary(coryne17$urealyticum)
(length(which(coryne17$urealyticum!=0))/nrow(coryne17))

phy.coryne18 <- subset_taxa(phy.coryne, blast_species=="jeddahense")
ntaxa(phy.coryne18)
coryne18 <- data.frame(otu_table(phy.coryne18))
coryne18 <- data.frame(t(coryne18), check.names = FALSE)
coryne18$jeddahense <- rowSums(coryne18)
coryne18$SampleID <- row.names(coryne18)
coryne18 <- coryne18[,c("SampleID", "jeddahense")]
summary(coryne18$jeddahense)
(length(which(coryne18$jeddahense!=0))/nrow(coryne18))

phy.coryne19 <- subset_taxa(phy.coryne, blast_species=="riegelii")
ntaxa(phy.coryne19)
coryne19 <- data.frame(otu_table(phy.coryne19))
coryne19 <- data.frame(t(coryne19), check.names = FALSE)
coryne19$riegelii <- rowSums(coryne19)
coryne19$SampleID <- row.names(coryne19)
coryne19 <- coryne19[,c("SampleID", "riegelii")]
summary(coryne19$riegelii)
(length(which(coryne19$riegelii!=0))/nrow(coryne19))

phy.coryne20 <- subset_taxa(phy.coryne, blast_species=="thomssenii_sundsvallense")
ntaxa(phy.coryne20)
coryne20 <- data.frame(otu_table(phy.coryne20))
coryne20 <- data.frame(t(coryne20), check.names = FALSE)
coryne20$thomssenii_sundsvallense <- rowSums(coryne20)
coryne20$SampleID <- row.names(coryne20)
coryne20 <- coryne20[,c("SampleID", "thomssenii_sundsvallense")]
summary(coryne20$thomssenii_sundsvallense)
(length(which(coryne20$thomssenii_sundsvallense!=0))/nrow(coryne20))

phy.coryne21 <- subset_taxa(phy.coryne, blast_species=="gottingense_imitans")
ntaxa(phy.coryne21)
coryne21 <- data.frame(otu_table(phy.coryne21))
coryne21 <- data.frame(t(coryne21), check.names = FALSE)
coryne21$gottingense_imitans <- rowSums(coryne21)
coryne21$SampleID <- row.names(coryne21)
coryne21 <- coryne21[,c("SampleID", "gottingense_imitans")]
summary(coryne21$gottingense_imitans)
(length(which(coryne21$gottingense_imitans!=0))/nrow(coryne21))

phy.coryne22 <- subset_taxa(phy.coryne, blast_species=="lipophiloflavum")
ntaxa(phy.coryne22)
coryne22 <- data.frame(otu_table(phy.coryne22))
coryne22 <- data.frame(t(coryne22), check.names = FALSE)
coryne22$lipophiloflavum <- rowSums(coryne22)
coryne22$SampleID <- row.names(coryne22)
coryne22 <- coryne22[,c("SampleID", "lipophiloflavum")]
summary(coryne22$lipophiloflavum)
(length(which(coryne22$lipophiloflavum!=0))/nrow(coryne22))

phy.coryne23 <- subset_taxa(phy.coryne, blast_species=="dorum")
ntaxa(phy.coryne23)
coryne23 <- data.frame(otu_table(phy.coryne23))
coryne23 <- data.frame(t(coryne23), check.names = FALSE)
coryne23$dorum <- rowSums(coryne23)
coryne23$SampleID <- row.names(coryne23)
coryne23 <- coryne23[,c("SampleID", "dorum")]
summary(coryne23$dorum)
(length(which(coryne23$dorum!=0))/nrow(coryne23))

phy.coryne24 <- subset_taxa(phy.coryne, blast_species=="glutamicum_efficiens_faecale")
ntaxa(phy.coryne24)
coryne24 <- data.frame(otu_table(phy.coryne24))
coryne24 <- data.frame(t(coryne24), check.names = FALSE)
coryne24$glutamicum_efficiens_faecale <- rowSums(coryne24)
coryne24$SampleID <- row.names(coryne24)
coryne24 <- coryne24[,c("SampleID", "glutamicum_efficiens_faecale")]
summary(coryne24$glutamicum_efficiens_faecale)
(length(which(coryne24$glutamicum_efficiens_faecale!=0))/nrow(coryne24))

phy.coryne25 <- subset_taxa(phy.coryne, blast_species=="massiliense")
ntaxa(phy.coryne25)
coryne25 <- data.frame(otu_table(phy.coryne25))
coryne25 <- data.frame(t(coryne25), check.names = FALSE)
coryne25$massiliense <- rowSums(coryne25)
coryne25$SampleID <- row.names(coryne25)
coryne25 <- coryne25[,c("SampleID", "massiliense")]
summary(coryne25$massiliense)
(length(which(coryne25$massiliense!=0))/nrow(coryne25))

phy.coryne26 <- subset_taxa(phy.coryne, blast_species=="terpenotabidum_variabile_freiburgense")
ntaxa(phy.coryne26)
coryne26 <- data.frame(otu_table(phy.coryne26))
coryne26 <- data.frame(t(coryne26), check.names = FALSE)
coryne26$terpenotabidum_variabile_freiburgense <- rowSums(coryne26)
coryne26$SampleID <- row.names(coryne26)
coryne26 <- coryne26[,c("SampleID", "terpenotabidum_variabile_freiburgense")]
summary(coryne26$terpenotabidum_variabile_freiburgense)
(length(which(coryne26$terpenotabidum_variabile_freiburgense!=0))/nrow(coryne26))

phy.coryne27 <- subset_taxa(phy.coryne, blast_species=="matruchotii")
ntaxa(phy.coryne27)
coryne27 <- data.frame(otu_table(phy.coryne27))
coryne27 <- data.frame(t(coryne27), check.names = FALSE)
coryne27$matruchotii <- rowSums(coryne27)
coryne27$SampleID <- row.names(coryne27)
coryne27 <- coryne27[,c("SampleID", "matruchotii")]
summary(coryne27$matruchotii)
(length(which(coryne27$matruchotii!=0))/nrow(coryne27))

phy.coryne28 <- subset_taxa(phy.coryne, blast_species=="ammoniagenes_casei_stationis")
ntaxa(phy.coryne28)
coryne28 <- data.frame(otu_table(phy.coryne28))
coryne28 <- data.frame(t(coryne28), check.names = FALSE)
coryne28$ammoniagenes_casei_stationis <- rowSums(coryne28)
coryne28$SampleID <- row.names(coryne28)
coryne28 <- coryne28[,c("SampleID", "ammoniagenes_casei_stationis")]
summary(coryne28$ammoniagenes_casei_stationis)
(length(which(coryne28$ammoniagenes_casei_stationis!=0))/nrow(coryne28))

phy.coryne29 <- subset_taxa(phy.coryne, blast_species=="bovis")
ntaxa(phy.coryne29)
coryne29 <- data.frame(otu_table(phy.coryne29))
coryne29 <- data.frame(t(coryne29), check.names = FALSE)
coryne29$bovis <- rowSums(coryne29)
coryne29$SampleID <- row.names(coryne29)
coryne29 <- coryne29[,c("SampleID", "bovis")]
summary(coryne29$bovis)
(length(which(coryne29$bovis!=0))/nrow(coryne29))

phy.coryne30 <- subset_taxa(phy.coryne, blast_species=="suicordis")
ntaxa(phy.coryne30)
coryne30 <- data.frame(otu_table(phy.coryne30))
coryne30 <- data.frame(t(coryne30), check.names = FALSE)
coryne30$suicordis <- rowSums(coryne30)
coryne30$SampleID <- row.names(coryne30)
coryne30 <- coryne30[,c("SampleID", "suicordis")]
summary(coryne30$suicordis)
(length(which(coryne30$suicordis!=0))/nrow(coryne30))

phy.coryne31 <- subset_taxa(phy.coryne, blast_species=="tuscaniense")
ntaxa(phy.coryne31)
coryne31 <- data.frame(otu_table(phy.coryne31))
coryne31 <- data.frame(t(coryne31), check.names = FALSE)
coryne31$tuscaniense <- rowSums(coryne31)
coryne31$SampleID <- row.names(coryne31)
coryne31 <- coryne31[,c("SampleID", "tuscaniense")]
summary(coryne31$tuscaniense)
(length(which(coryne31$tuscaniense!=0))/nrow(coryne31))

phy.coryne32 <- subset_taxa(phy.coryne, blast_species=="confusum")
ntaxa(phy.coryne32)
coryne32 <- data.frame(otu_table(phy.coryne32))
coryne32 <- data.frame(t(coryne32), check.names = FALSE)
coryne32$confusum <- rowSums(coryne32)
coryne32$SampleID <- row.names(coryne32)
coryne32 <- coryne32[,c("SampleID", "confusum")]
summary(coryne32$confusum)
(length(which(coryne32$confusum!=0))/nrow(coryne32))

phy.coryne33 <- subset_taxa(phy.coryne, blast_species=="pseudotuberculosis_ulcerans")
ntaxa(phy.coryne33)
coryne33 <- data.frame(otu_table(phy.coryne33))
coryne33 <- data.frame(t(coryne33), check.names = FALSE)
coryne33$pseudotuberculosis_ulcerans <- rowSums(coryne33)
coryne33$SampleID <- row.names(coryne33)
coryne33 <- coryne33[,c("SampleID", "pseudotuberculosis_ulcerans")]
summary(coryne33$pseudotuberculosis_ulcerans)
(length(which(coryne33$pseudotuberculosis_ulcerans!=0))/nrow(coryne33))

phy.coryne34 <- subset_taxa(phy.coryne, blast_species=="timonense")
ntaxa(phy.coryne34)
coryne34 <- data.frame(otu_table(phy.coryne34))
coryne34 <- data.frame(t(coryne34), check.names = FALSE)
coryne34$timonense <- rowSums(coryne34)
coryne34$SampleID <- row.names(coryne34)
coryne34 <- coryne34[,c("SampleID", "timonense")]
summary(coryne34$timonense)
(length(which(coryne34$timonense!=0))/nrow(coryne34))

phy.coryne35 <- subset_taxa(phy.coryne, blast_species=="mycetoides_riegelii")
ntaxa(phy.coryne35)
coryne35 <- data.frame(otu_table(phy.coryne35))
coryne35 <- data.frame(t(coryne35), check.names = FALSE)
coryne35$mycetoides_riegelii <- rowSums(coryne35)
coryne35$SampleID <- row.names(coryne35)
coryne35 <- coryne35[,c("SampleID", "mycetoides_riegelii")]
summary(coryne35$mycetoides_riegelii)
(length(which(coryne35$mycetoides_riegelii!=0))/nrow(coryne35))

phy.coryne36 <- subset_taxa(phy.coryne, blast_species=="glaucum")
ntaxa(phy.coryne36)
coryne36 <- data.frame(otu_table(phy.coryne36))
coryne36 <- data.frame(t(coryne36), check.names = FALSE)
coryne36$glaucum <- rowSums(coryne36)
coryne36$SampleID <- row.names(coryne36)
coryne36 <- coryne36[,c("SampleID", "glaucum")]
summary(coryne36$glaucum)
(length(which(coryne36$glaucum!=0))/nrow(coryne36))

phy.coryne37 <- subset_taxa(phy.coryne, blast_species=="humireducens")
ntaxa(phy.coryne37)
coryne37 <- data.frame(otu_table(phy.coryne37))
coryne37 <- data.frame(t(coryne37), check.names = FALSE)
coryne37$humireducens <- rowSums(coryne37)
coryne37$SampleID <- row.names(coryne37)
coryne37 <- coryne37[,c("SampleID", "humireducens")]
summary(coryne37$humireducens)
(length(which(coryne37$humireducens!=0))/nrow(coryne37))

phy.coryne38 <- subset_taxa(phy.coryne, blast_species=="falsenii")
ntaxa(phy.coryne38)
coryne38 <- data.frame(otu_table(phy.coryne38))
coryne38 <- data.frame(t(coryne38), check.names = FALSE)
coryne38$falsenii <- rowSums(coryne38)
coryne38$SampleID <- row.names(coryne38)
coryne38 <- coryne38[,c("SampleID", "falsenii")]
summary(coryne38$falsenii)
(length(which(coryne38$falsenii!=0))/nrow(coryne38))

phy.coryne39 <- subset_taxa(phy.coryne, blast_species=="pyruviciproducens")
ntaxa(phy.coryne39)
coryne39 <- data.frame(otu_table(phy.coryne39))
coryne39 <- data.frame(t(coryne39), check.names = FALSE)
coryne39$pyruviciproducens <- rowSums(coryne39)
coryne39$SampleID <- row.names(coryne39)
coryne39 <- coryne39[,c("SampleID", "pyruviciproducens")]
summary(coryne39$pyruviciproducens)
(length(which(coryne39$pyruviciproducens!=0))/nrow(coryne39))

phy.coryne40 <- subset_taxa(phy.coryne, blast_species=="camporealensis_canis_epidermidicanis")
ntaxa(phy.coryne40)
coryne40 <- data.frame(otu_table(phy.coryne40))
coryne40 <- data.frame(t(coryne40), check.names = FALSE)
coryne40$camporealensis_canis_epidermidicanis <- rowSums(coryne40)
coryne40$SampleID <- row.names(coryne40)
coryne40 <- coryne40[,c("SampleID", "camporealensis_canis_epidermidicanis")]
summary(coryne40$camporealensis_canis_epidermidicanis)
(length(which(coryne40$camporealensis_canis_epidermidicanis!=0))/nrow(coryne40))

phy.coryne41 <- subset_taxa(phy.coryne, blast_species=="glucuronolyticum")
ntaxa(phy.coryne41)
coryne41 <- data.frame(otu_table(phy.coryne41))
coryne41 <- data.frame(t(coryne41), check.names = FALSE)
coryne41$glucuronolyticum <- rowSums(coryne41)
coryne41$SampleID <- row.names(coryne41)
coryne41 <- coryne41[,c("SampleID", "glucuronolyticum")]
summary(coryne41$glucuronolyticum)
(length(which(coryne41$glucuronolyticum!=0))/nrow(coryne41))

phy.coryne42 <- subset_taxa(phy.coryne, blast_species=="testudinoris")
ntaxa(phy.coryne42)
coryne42 <- data.frame(otu_table(phy.coryne42))
coryne42 <- data.frame(t(coryne42), check.names = FALSE)
coryne42$testudinoris <- rowSums(coryne42)
coryne42$SampleID <- row.names(coryne42)
coryne42 <- coryne42[,c("SampleID", "testudinoris")]
summary(coryne42$testudinoris)
(length(which(coryne42$testudinoris!=0))/nrow(coryne42))

phy.coryne43 <- subset_taxa(phy.coryne, blast_species=="glyciniphilum")
ntaxa(phy.coryne43)
coryne43 <- data.frame(otu_table(phy.coryne43))
coryne43 <- data.frame(t(coryne43), check.names = FALSE)
coryne43$glyciniphilum <- rowSums(coryne43)
coryne43$SampleID <- row.names(coryne43)
coryne43 <- coryne43[,c("SampleID", "glyciniphilum")]
summary(coryne43$glyciniphilum)
(length(which(coryne43$glyciniphilum!=0))/nrow(coryne43))

phy.coryne44 <- subset_taxa(phy.coryne, blast_species=="resistens_auriscanis")
ntaxa(phy.coryne44)
coryne44 <- data.frame(otu_table(phy.coryne44))
coryne44 <- data.frame(t(coryne44), check.names = FALSE)
coryne44$resistens_auriscanis <- rowSums(coryne44)
coryne44$SampleID <- row.names(coryne44)
coryne44 <- coryne44[,c("SampleID", "resistens_auriscanis")]
summary(coryne44$resistens_auriscanis)
(length(which(coryne44$resistens_auriscanis!=0))/nrow(coryne44))

phy.coryne_oth <- subset_taxa(phy.coryne, blast_species=="unknown")
ntaxa(phy.coryne_oth)
coryne_oth <- data.frame(otu_table(phy.coryne_oth))
coryne_oth <- data.frame(t(coryne_oth), check.names = FALSE)
coryne_oth$other <- rowSums(coryne_oth)
coryne_oth$SampleID <- row.names(coryne_oth)
coryne_oth <- coryne_oth[,c("SampleID", "other")]
summary(coryne_oth$other)
(length(which(coryne_oth$other!=0))/nrow(coryne_oth))

coryne_top <- merge(coryne_all, coryne1, by="SampleID")
coryne_top <- merge(coryne_top, coryne2, by="SampleID")
coryne_top <- merge(coryne_top, coryne3, by="SampleID")
coryne_top <- merge(coryne_top, coryne4, by="SampleID")
coryne_top <- merge(coryne_top, coryne5, by="SampleID")
coryne_top <- merge(coryne_top, coryne6, by="SampleID")
row.names(coryne_top) <- coryne_top[,1]
coryne_top <- coryne_top[,-1]
coryne_top <- tibble::rownames_to_column(coryne_top, "SampleID")
coryne_top$SampleID <- gsub( "X", "", coryne_top$SampleID)
coryne_top$SampleID <- gsub("[.]", "-", coryne_top$SampleID)

# Dolosigranulum pigrum

phy.dolo <- subset_taxa(phy.np.inf.relative_asv, genus=="Dolosigranulum")
ntaxa(phy.dolo)
dolo_all <- data.frame(otu_table(phy.dolo))
dolo_all <- data.frame(t(dolo_all), check.names = FALSE)
dolo_all$dolo_ALL <- rowSums(dolo_all)
dolo_all$SampleID2 <- row.names(dolo_all)
dolo_all <- dolo_all[,c("SampleID2", "dolo_ALL")]
dolo_all <- tibble::rownames_to_column(dolo_all, "SampleID")
dolo_all$SampleID <- gsub( "X", "", dolo_all$SampleID)
dolo_all$SampleID <- gsub("[.]", "-", dolo_all$SampleID)
dolo_all <- dolo_all[,-2]

biotypes_final <- merge(biotypes, coryne_top, by="SampleID")
biotypes_final <- merge(biotypes_final, dolo_all, by="SampleID")
remove(coryne1, coryne2, coryne3, coryne4, coryne5, coryne6, coryne7, coryne8, coryne9, coryne10, coryne11, coryne12, coryne13, 
       coryne14, coryne15, coryne16, coryne17, coryne18, coryne19, coryne20, coryne21, coryne22, coryne23, coryne24, coryne25,
       coryne26, coryne27, coryne28, coryne29, coryne30, coryne31, coryne32, coryne33, coryne34, coryne35, coryne36, coryne37,
       coryne38, coryne39, coryne40, coryne41, coryne42, coryne43, coryne44, coryne_oth)
remove(phy.coryne1, phy.coryne2, phy.coryne3, phy.coryne4, phy.coryne5, phy.coryne6, phy.coryne7, phy.coryne8, phy.coryne9, 
       phy.coryne10, phy.coryne11, phy.coryne12, phy.coryne13, 
       phy.coryne14, phy.coryne15, phy.coryne16, phy.coryne17, phy.coryne18, phy.coryne19, phy.coryne20, phy.coryne21, 
       phy.coryne22, phy.coryne23, phy.coryne24, phy.coryne25,
       phy.coryne26, phy.coryne27, phy.coryne28, phy.coryne29, phy.coryne30, phy.coryne31, phy.coryne32, phy.coryne33, 
       phy.coryne34, phy.coryne35, phy.coryne36, phy.coryne37,
       phy.coryne38, phy.coryne39, phy.coryne40, phy.coryne41, phy.coryne42, phy.coryne43, phy.coryne44, phy.coryne_oth)
remove(coryne_all, coryne_top, dolo_all, phy.np.inf.relative_asv, phy.coryne, phy.dolo)

write.csv(biotypes_final, "biotypes.csv")

# ********************************************************************************************************************** 
# ANALYSES OF MICROBIOME STABILITY
# ********************************************************************************************************************** 

lab_data <- read.csv("metadata_all.csv", stringsAsFactors = TRUE)
lab_data <- lab_data[,-c(1,2,69)]
lab_data <- subset(lab_data, Source=="INF" & SampleType=="NPS")
np_data <- biotypes_final[,c(1:5,68,70:81)]
combined <- merge(lab_data, np_data, by=c("study_id", "month", "Source", "SampleType"), all.x=TRUE)

combined$inf_sp_yn[combined$inf_sp>0] <- 1
combined$inf_sp_yn[combined$inf_sp==0] <- 0

combined$recent_uri_yn[combined$recent_uri=="Y"] <- 1
combined$recent_uri_yn[combined$recent_uri=="N"] <- 0

combined$pneumonia_yn[combined$pneumonia=="Y"] <- 1
combined$pneumonia_yn[combined$pneumonia=="N"] <- 0

combined$season2[combined$season=="Rainy"] <- "Summer"
combined$season2[combined$season=="Dry"] <- "Winter"

combined$pcv[combined$age_days<(combined$pcv1+15)] <- "0"  
combined$pcv[combined$age_days>=(combined$pcv1+15) & combined$age_days<(combined$pcv2+15)] <- "1"
combined$pcv[combined$age_days>=(combined$pcv2+15) & combined$age_days<(combined$pcv3+15)] <- "2"
combined$pcv[combined$age_days>=(combined$pcv3+15)] <- "3"
combined$pcv[combined$age_days<60 & is.na(combined$pcv)] <- "0"
combined$pcv <- as.numeric(combined$pcv)
table(combined$pcv, useNA="always")
combined <- combined[ , -which(names(combined) %in% c("pcv1","pcv2", "pcv3"))]

# Infant RV infections
combined$inf_rv_yn[combined$inf_rv=="NEG"] <- "N"
combined$inf_rv_yn[combined$inf_rv!="NEG" & !is.na(combined$inf_rv)] <- "Y"

# Maternal Sp colonization
combined$mom_sp_yn[combined$mom_sp>0] <- "Y"
combined$mom_sp_yn[combined$mom_sp==0] <- "N"

vars_to_keep <- c("SampleID", "study_id", "month", "Source", "SampleType", "age_days", "biotype")
stability <- combined[,(vars_to_keep)]
stability$biotype[is.na(stability$biotype)] <- "Missing"
stability <- suppressWarnings(slide(stability, "biotype", TimeVar="month", GroupVar="study_id", NewVar="biotype_lead", slideBy = 1))
stability <- subset(stability, biotype!="Missing" & !is.na(biotype_lead))
stability$stable[stability$biotype=="Biodiverse" & stability$biotype_lead=="Biodiverse"] <- "0"
stability$stable[stability$biotype=="Corynebacterium" & stability$biotype_lead=="Corynebacterium"] <- "0"
stability$stable[stability$biotype=="Dolosigranulum" & stability$biotype_lead=="Dolosigranulum"] <- "0"
stability$stable[stability$biotype=="Haemophilus" & stability$biotype_lead=="Haemophilus"] <- "0"
stability$stable[stability$biotype=="Moraxella" & stability$biotype_lead=="Moraxella"] <- "0"
stability$stable[stability$biotype=="Other" & stability$biotype_lead=="Other"] <- "0"
stability$stable[stability$biotype=="Staphylococcus" & stability$biotype_lead=="Staphylococcus"] <- "0"
stability$stable[stability$biotype=="Streptococcus" & stability$biotype_lead=="Streptococcus"] <- "0"
stability$stable[is.na(stability$stable)] <- "1"
stability$biotype <- as.factor(stability$biotype)
stability$stable <- as.factor(stability$stable)
table(stability$biotype)
table(stability$biotype, stability$stable)
prop.table(table(stability$biotype, stability$stable),1)
logit_biotype <- glm(stable ~ biotype + age_days, family=binomial(link='logit'), data=stability)
summary(logit_biotype)
exp(cbind(OR = coef(logit_biotype), confint(logit_biotype)))
remove(logit_biotype, stability, vars_to_keep)

# ********************************************************************************************************************** 
# ANALYSES OF STREP PNEUMO COLONIZATION
# ********************************************************************************************************************** 

sp_df <- subset(combined, study_id %in% infants)
table(sp_df$inf_sp_yn, useNA="always") # Missing Sp colonization data for 10 samples; all of these infants were colonized at >=1 other time points
sp_df <- subset(sp_df, !is.na(inf_sp_yn))
length(unique(sp_df$study_id))
sp_col_df <- subset(sp_df, inf_sp_yn=="1")
length(unique(sp_col_df$study_id))
length(unique(sp_col_df$study_id))/length(unique(sp_df$study_id))
sp_col_df <- sp_col_df %>% group_by(study_id) %>% filter(row_number()==1)
summary(sp_col_df$age_days)

# Create cumulative incidence curve for Strep pneumo colonization

sp_df <- sp_df[,c("study_id","age_days","inf_rv","inf_sp_yn")]

# Create dataset with last follow-up visit
DT <- sp_df
DT_exit <- DT[,-c(3,4)]
DT_exit <- as.data.table(DT_exit)[,.SD[which.max(age_days)], study_id]

# Create dataset with first Sp colonization event
DT_sp_pos <- DT
DT_sp_pos <- filter(DT_sp_pos, inf_sp_yn==1)
DT_sp_pos <- as.data.table(DT_sp_pos)[,.SD[which.min(age_days)], study_id]
DT_sp_pos$exit <- DT_sp_pos$age_days
DT_sp_pos <- DT_sp_pos[,-c(2,3)]

# Merge these datasets
DT_merge <- merge(DT_exit, DT_sp_pos, all.x=TRUE)
DT_merge$cause[is.na(DT_merge$exit)] <- 2

DT_merge$cause[DT_merge$exit>=0] <- 1
DT_merge$exit[is.na(DT_merge$exit)] <- DT_merge$age_days[is.na(DT_merge$exit)]
DT_merge <- DT_merge[,-2]
DT_merge$entry <- 0
DT_merge$group <- 0
DT_merge <- DT_merge[DT_merge$exit>0]
DT_merge <- subset(DT_merge, cause==1 | exit>=6)

require(etm)
DT_df <- as.data.frame.matrix(DT_merge)
cif.sp <- etmCIF(survival::Surv(entry, exit, cause != 0) ~ group, data = DT_df, etype = cause, failcode = 1)

strep <- expression(paste("Proportion of infants colonized with ", italic("S. pneumoniae")))

library(Hmisc)
png(file="Supplemental_Figure_1.png", 
    width = 6, height = 6, units = 'in', res = 600)
plot(cif.sp, xlab="Infant age (days)", ylab=strep, col="blue", lwd=2, ylim=c(0,1), xlim=c(0,400))
minor.tick(ny=2,nx=4)
dev.off()

remove(cif.sp, DT, DT_df, DT_exit, DT_merge, DT_sp_pos, strep, sp_col_df, sp_df)

# ********************************************************************************************************************** 
# Aim 2: Determine if nasopharyngeal microbiome diversity and composition are associated with subsequent 
# colonization by S. pneumoniae and other bacterial respiratory pathogens. 
# ********************************************************************************************************************** 

library(survival)
library(survminer)

# Create dataset for survival analyses
survival <- combined[,c("study_id", "month", "sex", "mat_hiv", "wood", "lbw", "breastmilk", "season2", "residence", "inf_abx_any",
                        "age_days", "pcv", "mom_sp_yn", "inf_rv_yn", "num_kids",
                        "inf_sp_yn", "recent_uri_yn", "pneumonia_yn", "biotype", "coryne_ALL", "accolens_macginleyi",
                        "pseudodiphtheriticum_propinquum", "tuberculostearicum", "dolo_ALL")]

# **************************
# S. pneumoniae colonization
# **************************

# Remove infants who were colonized with Sp at birth
nrow(survival)
table(survival$month, survival$inf_sp_yn, useNA="ifany")
drop_df <- subset(survival, month=="0" & inf_sp_yn=="1")
drop <- drop_df$study_id
survival_sp <- subset(survival, !(study_id %in% drop))
remove(drop_df, drop)
nrow(survival_sp)

# Lead NP microbiome data
survival_sp$study_id <- as.factor(survival_sp$study_id)
survival_sp <- survival_sp[order(survival_sp$study_id, survival_sp$month),]
survival_sp <- suppressWarnings(slide(survival_sp, "age_days", TimeVar="month", GroupVar="study_id", NewVar="age_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "biotype", TimeVar="month", GroupVar="study_id", NewVar="biotype_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "coryne_ALL", TimeVar="month", GroupVar="study_id", NewVar="coryne_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "accolens_macginleyi", TimeVar="month", 
                                      GroupVar="study_id", NewVar="accolens_macginleyi_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "pseudodiphtheriticum_propinquum", TimeVar="month", 
                                      GroupVar="study_id", NewVar="pseudodiphtheriticum_propinquum_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "tuberculostearicum", TimeVar="month", 
                                      GroupVar="study_id", NewVar="tuberculostearicum_lag", slideBy = -1))
survival_sp <- suppressWarnings(slide(survival_sp, "dolo_ALL", TimeVar="month", GroupVar="study_id", NewVar="dolo_lag", slideBy = -1))
names(survival_sp)[names(survival_sp) == "study_id"] <- "id"
names(survival_sp)[names(survival_sp) == "age_days"] <- "stop"
survival_sp$stop <- as.numeric(survival_sp$stop)
names(survival_sp)[names(survival_sp) == "age_lag"] <- "start"
names(survival_sp)[names(survival_sp) == "inf_sp_yn"] <- "colonized"
survival_sp <- subset(survival_sp, month!="0") # Do not predict outcome at month=0
nrow(survival_sp)

survival_sp_final <- group_by(survival_sp, id) %>% filter(cumsum(cumsum(colonized == 1))<2)
table(survival_sp_final$biotype_lag, useNA="always")
survival_sp_final$biotype_lag[is.na(survival_sp_final$biotype_lag)] <- "Missing"
table(survival_sp_final$biotype_lag, useNA="always")
survival_sp_final$pcv <- as.numeric(survival_sp_final$pcv)
survival_sp_final$pcv[survival_sp_final$start<=60 & is.na(survival_sp_final$pcv)] <- 0
survival_sp_final <- survival_sp_final[,c("id", "start", "stop", "sex", "mat_hiv", "wood", "lbw", "breastmilk", "season2", "residence", 
                                          "inf_abx_any", "pcv", "num_kids",
                                          "inf_rv_yn", "mom_sp_yn", "colonized",  
                                          "biotype_lag", "coryne_lag", "accolens_macginleyi_lag", "pseudodiphtheriticum_propinquum_lag", 
                                          "tuberculostearicum_lag", "dolo_lag")]

table(survival_sp_final$pcv, useNA="always")
table(survival_sp_final$inf_rv_yn, useNA="always")
summary(survival_sp_final$num_kids)

# BIOTYPE

table(survival_sp_final$biotype_lag, survival_sp_final$colonized)
prop.table(table(survival_sp_final$biotype_lag, survival_sp_final$colonized),1)
# Use of cluster function clusters the observations for the purpose of a robust variance
cox_biotype <- coxph(Surv(start, stop, colonized) ~ biotype_lag + 
                  sex + lbw + mat_hiv + residence + wood + num_kids + season2 + breastmilk + pcv + inf_abx_any + 
                  cluster(id), data=survival_sp_final)
summary(cox_biotype)
cox.zph(cox_biotype) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)

# CORYNEBACTERIUM RELATIVE ABUNDANCE IN QUARTILES

suppressWarnings(survival_sp_final$coryne_quartile <- NA)
survival_sp_final$coryne_quartile <- with(survival_sp_final, cut(coryne_lag, 
                                                 breaks=quantile(coryne_lag, probs=seq(0,1, by=0.25), na.rm=TRUE), 
                                                 include.lowest=TRUE))
levels(survival_sp_final$coryne_quartile)
table(survival_sp_final$coryne_quartile, survival_sp_final$colonized)
prop.table(table(survival_sp_final$coryne_quartile, survival_sp_final$colonized),1)
prop.table(table(survival_sp_final$coryne_quartile, survival_sp_final$colonized),1)
cox_coryneq <- coxph(Surv(start, stop, colonized) ~ coryne_quartile + 
                  sex + lbw + mat_hiv + residence + wood + num_kids + season2 + breastmilk + pcv + inf_abx_any + 
                  cluster(id), data=survival_sp_final)
summary(cox_coryneq)
cox.zph(cox_coryneq) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)

# SPECIFIC CORYNEBACTERIUM SPECIES

survival_sp_final$pseudodiphtheriticum_propinquum_lag10 <- (survival_sp_final$pseudodiphtheriticum_propinquum_lag)*10
survival_sp_final$accolens_macginleyi_lag10 <- (survival_sp_final$accolens_macginleyi_lag)*10
survival_sp_final$tuberculostearicum_lag10 <- (survival_sp_final$tuberculostearicum_lag)*10
cox_species <- coxph(Surv(start, stop, colonized) ~ pseudodiphtheriticum_propinquum_lag10 + accolens_macginleyi_lag10 + 
                       tuberculostearicum_lag10 +
                       sex + lbw + mat_hiv + residence + wood + num_kids + season2 + breastmilk + pcv + inf_abx_any + 
                  cluster(id), data=survival_sp_final)
summary(cox_species)
cox.zph(cox_species) # test of the proportional hazards assumption (P<0.05 is indicative of non-proportional hazards)

remove(lab_data, np_data, cox_biotype, cox_coryneq, cox_species, survival_sp, biotypes)

# ********************************************************************************************************************** 
# Aim 1: Evaluate associations between external factors and the composition of the nasopharyngeal microbiome.
# ********************************************************************************************************************** 

library(Maaslin2)
library(lme4)
options(warn=-1)

combined$sample[is.na(combined$SampleID)] <- "0"
combined$sample[!is.na(combined$SampleID)] <- "1"
combined$sample <- as.numeric(combined$sample)
maaslin <- group_by(combined, study_id) %>% mutate(visit = cumsum(sample))

# Create variable "interval" for time periods of interest
maaslin$interval[(maaslin$sample==0 & maaslin$visit==0) | (maaslin$sample==1 & maaslin$visit==1)] <- "1"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==1) | (maaslin$sample==1 & maaslin$visit==2)] <- "2"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==2) | (maaslin$sample==1 & maaslin$visit==3)] <- "3"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==3) | (maaslin$sample==1 & maaslin$visit==4)] <- "4"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==4) | (maaslin$sample==1 & maaslin$visit==5)] <- "5"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==5) | (maaslin$sample==1 & maaslin$visit==6)] <- "6"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==6) | (maaslin$sample==1 & maaslin$visit==7)] <- "7"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==7) | (maaslin$sample==1 & maaslin$visit==8)] <- "8"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==8) | (maaslin$sample==1 & maaslin$visit==9)] <- "9"
maaslin$interval[(maaslin$sample==0 & maaslin$visit==9) | (maaslin$sample==1 & maaslin$visit==10)] <- "10"

maaslin$breastmilk_num <- NA
maaslin$breastmilk_num[maaslin$breastmilk=="N"] <- 0
maaslin$breastmilk_num[maaslin$breastmilk=="Y"] <- 1
maaslin <- group_by(maaslin, study_id, interval) %>% mutate(bm_max = max(breastmilk_num))
maaslin$bm_max[maaslin$bm_max=="0"] <- "N"
maaslin$bm_max[maaslin$bm_max=="1"] <- "Y"
table(maaslin$breastmilk, maaslin$bm_max)

maaslin$season2_num <- NA
maaslin$season2_num[maaslin$season2=="Summer"] <- 0
maaslin$season2_num[maaslin$season2=="Winter"] <- 1
maaslin <- group_by(maaslin, study_id, interval) %>% mutate(season2_max = max(season2_num))
maaslin$season2_max[maaslin$season2_max=="0"] <- "Summer"
maaslin$season2_max[maaslin$season2_max=="1"] <- "Winter"
table(maaslin$season2, maaslin$season2_max)

maaslin$abx_num <- NA
maaslin$abx_num[maaslin$inf_abx_any=="N"] <- 0
maaslin$abx_num[maaslin$inf_abx_any=="Y"] <- 1
maaslin <- group_by(maaslin, study_id, interval) %>% mutate(abx_max = max(abx_num))
maaslin$abx_max[maaslin$abx_max=="0"] <- "N"
maaslin$abx_max[maaslin$abx_max=="1"] <- "Y"
table(maaslin$inf_abx_any, maaslin$abx_max)

maaslin$kids_cat[maaslin$num_kids<=1] <- "<=1"
maaslin$kids_cat[maaslin$num_kids>=2] <- ">=2"

maaslin <- maaslin[ , -which(names(maaslin) %in% c("breastmilk_num", "season2_num", "abx_num"))]
maaslin <- subset(maaslin, !is.na(SampleID))
maaslin_df <- data.frame(maaslin[,c("study_id", "month", "interval", "SampleID", "Source", "SampleType", "calendar_mo", 
                                "age_days", "sex", "lbw", "mat_hiv", "residence", "wood", "inf_rv_yn", "residence", "pcv", "kids_cat",
                                "bm_max", "season2_max", "abx_max")])
row.names(maaslin_df) <- maaslin_df$SampleID
remove(maaslin)

# Agglomerate ASVs on the genus level
phy.np.agglom.inf <- tax_glom(phy.np.normalized.inf, "genus")
nsamples(phy.np.agglom.inf)
ntaxa(phy.np.agglom.inf)
sample_data(phy.np.agglom.inf) <- maaslin_df

# Transform to relative abundances
ntaxa(phy.np.agglom.inf)
phy.maaslin <- transform_sample_counts(phy.np.agglom.inf, function(Abundance) Abundance/sum(Abundance))
head(sample_sums(phy.maaslin))
ntaxa(phy.maaslin)

input_data_df <- data.frame(otu_table(phy.maaslin))
input_data_df2 <- merge(taxtable, input_data_df, by="row.names", all.x=FALSE)
row.names(input_data_df2) <- input_data_df2$genus
input_data_df2 <- input_data_df2[,-c(1:9)]
input_metadata_df <- data.frame(sample_data(phy.maaslin)) 
row.names(input_metadata_df) <- gsub("-",".", rownames(input_metadata_df))

# SEX - MAASLIN2

input_metadata_df$sex = factor(input_metadata_df$sex,levels = c("sex", "F", "M"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/sex', 
  fixed_effects = c("sex"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
# No association between sex and the relative abundances of bacterial genera

# LBW - MAASLIN2

input_metadata_df$lbw = factor(input_metadata_df$lbw,levels = c("lbw", "N", "Y"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/lbw', 
  fixed_effects = c("lbw"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
# No association between LBW and the relative abundances of bacterial genera

# HIV EXPOSURE STATUS

input_metadata_df$mat_hiv = factor(input_metadata_df$mat_hiv,levels = c("mat_hiv", "N", "Y"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/mat_hiv', 
  fixed_effects = c("mat_hiv"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
# No association between HIV exposure status and the relative abundances of bacterial genera

# URBAN RESIDENCE

input_metadata_df$residence = factor(input_metadata_df$residence,levels = c("residence", "Rural", "Urban"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/residence', 
  fixed_effects = c("residence"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)

# HOUSEHOLD USE OF SOLID FUELS

input_metadata_df$wood = factor(input_metadata_df$wood,levels = c("wood", "N", "Y"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/wood', 
  fixed_effects = c("wood"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)

# NUMBER OF CHILD HOUSEHOLD MEMBERS

input_metadata_df$kids_cat = factor(input_metadata_df$kids_cat,levels = c("kids_cat", "<=1", ">=2"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/kids', 
  fixed_effects = c("kids_cat"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)

# SEASON

input_metadata_df$season2_max = factor(input_metadata_df$season2_max,levels = c("season2_max", "Summer", "Winter"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/season', 
  fixed_effects = c("season2_max"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
maaslin_season <- read.table(file = 'maaslin/season/significant_results.tsv', sep = '\t', header = TRUE)
maaslin_season$coef <- as.numeric(maaslin_season$coef)
maaslin_season$coef_pos[maaslin_season$coef>0] <- "Positive"
maaslin_season$coef_pos[maaslin_season$coef<0] <- "Negative"
maaslin_season$coef_pos <- as.factor(maaslin_season$coef_pos)
maaslin_season <- subset(maaslin_season, metadata=="season2_max")
maaslin_season <- ggplot(maaslin_season, aes(x=reorder(feature, -coef), y=coef, label=coef)) +
  geom_bar(stat='identity', aes(fill=coef_pos), width=.5) + 
  coord_flip() + scale_y_continuous(limits=c(-0.8,0.8), breaks=seq(-0.8,0.8,by=0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#003C67FF","#EFC000FF")) + ggtitle("Winter season") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5), axis.text.y=element_text(face = "italic")) + xlab("") + ylab("")
png(file="maaslin/maaslin_plot_season.png", width = 6, height = 4, units = 'in', res = 600)
plot(maaslin_season)
dev.off()

# BREASTFEEDING

input_metadata_df$bm_max = factor(input_metadata_df$bm_max,levels = c("bm_max", "N", "Y"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/breastfeeding', 
  fixed_effects = c("bm_max"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
maaslin_bm_max <- read.table(file = 'maaslin/breastfeeding/significant_results.tsv', sep = '\t', header = TRUE)
maaslin_bm_max$coef <- as.numeric(maaslin_bm_max$coef)
maaslin_bm_max$coef_pos[maaslin_bm_max$coef>0] <- "Positive"
maaslin_bm_max$coef_pos[maaslin_bm_max$coef<0] <- "Negative"
maaslin_bm_max$coef_pos <- as.factor(maaslin_bm_max$coef_pos)
maaslin_bm_max <- subset(maaslin_bm_max, metadata=="bm_max")
maaslin_bm_max <- ggplot(maaslin_bm_max, aes(x=reorder(feature, -coef), y=coef, label=coef)) +
  geom_bar(stat='identity', aes(fill=coef_pos), width=.5) + 
  coord_flip() + scale_y_continuous(limits=c(-0.8,0.8), breaks=seq(-0.8,0.8,by=0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#003C67FF","#EFC000FF")) + ggtitle("Infant breastfeeding") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5), axis.text.y=element_text(face = "italic")) + xlab("") + ylab("")
png(file="maaslin/maaslin_plot_bm_max.png", width = 6, height = 4, units = 'in', res = 600)
plot(maaslin_bm_max)
dev.off()

# PCV-13

fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/pcv', 
  fixed_effects = c("pcv", "age_days"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
maaslin_pcv <- read.table(file = 'maaslin/pcv/significant_results.tsv', sep = '\t', header = TRUE)
maaslin_pcv$coef <- as.numeric(maaslin_pcv$coef)
maaslin_pcv$coef_pos[maaslin_pcv$coef>0] <- "Positive"
maaslin_pcv$coef_pos[maaslin_pcv$coef<0] <- "Negative"
maaslin_pcv$coef_pos <- as.factor(maaslin_pcv$coef_pos)
maaslin_pcv <- subset(maaslin_pcv, metadata=="pcv")
maaslin_pcv <- ggplot(maaslin_pcv, aes(x=reorder(feature, -coef), y=coef, label=coef)) +
  geom_bar(stat='identity', aes(fill=coef_pos), width=.5) + 
  coord_flip() + scale_y_continuous(limits=c(-0.8,0.8), breaks=seq(-0.8,0.8,by=0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#003C67FF","#EFC000FF")) + ggtitle("Number of PCV-13 doses") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5), axis.text.y=element_text(face = "italic")) + xlab("") + ylab("")
png(file="maaslin/maaslin_plot_pcv.png", width = 6, height = 4, units = 'in', res = 600)
plot(maaslin_pcv)
dev.off()

# ANTIBIOTIC TREATMENT

input_metadata_df$abx_max = factor(input_metadata_df$abx_max,levels = c("abx_max", "N", "Y"))
fit_data <- Maaslin2(
  input_data_df2, input_metadata_df, 'maaslin/abx_max', 
  fixed_effects = c("abx_max"),
  random_effects = c("study_id"),
  min_prevalence = 0.1,
  max_significance = 0.2)
maaslin_abx_max <- read.table(file = 'maaslin/abx_max/significant_results.tsv', sep = '\t', header = TRUE)
maaslin_abx_max$coef <- as.numeric(maaslin_abx_max$coef)
maaslin_abx_max$coef_pos[maaslin_abx_max$coef>0] <- "Positive"
maaslin_abx_max$coef_pos[maaslin_abx_max$coef<0] <- "Negative"
maaslin_abx_max$coef_pos <- as.factor(maaslin_abx_max$coef_pos)
maaslin_abx_max <- subset(maaslin_abx_max, metadata=="abx_max")
maaslin_abx_max <- ggplot(maaslin_abx_max, aes(x=reorder(feature, -coef), y=coef, label=coef)) +
  geom_bar(stat='identity', aes(fill=coef_pos), width=.5) + 
  coord_flip() + scale_y_continuous(limits=c(-0.8,0.8), breaks=seq(-0.8,0.8,by=0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#003C67FF","#EFC000FF")) + ggtitle("Receipt of antibiotics") +
  theme(legend.position="none", plot.title=element_text(hjust=0.5), axis.text.y=element_text(face = "italic")) + xlab("") + ylab("")
png(file="maaslin/maaslin_plot_abx_max.png", width = 6, height = 4, units = 'in', res = 600)
plot(maaslin_abx_max)
dev.off()

png(file="R Plots/Figure_5.png", width = 14, height = 8, units = 'in', res = 600)
plot_grid(maaslin_abx_max, maaslin_pcv, maaslin_bm_max, maaslin_season, labels="auto", label_size=16, ncol=2) 
dev.off()

remove(maaslin_df, maaslin_abx_max, maaslin_pcv, maaslin_bm_max, maaslin_season, fit_data, input_data_df, input_data_df2, input_metadata_df)