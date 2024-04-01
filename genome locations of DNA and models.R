setwd("C:/study/RNA")
library(biomaRt)
library(minfi)
library(dplyr)
library(readxl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", force = TRUE)
BiocManager::install("minfi",force = TRUE)
BiocManager::install("plyranges")
BiocManager::install("qvalue")
BiocManager::install("sva")
BiocManager::install("GenomeInfoDb", force = TRUE)



annote450k <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annote450k$Islands_Name)
gene_names <- annote450k$UCSC_RefGene_Name

ph_hlhs <- readRDS('ph_hlhs.rds')
beta_hlhs <- readRDS('beta_hlhs.rds') #Methylation beta value
hlhs_data <- data.frame(read.table
                        (file = "35HLHSPatientMethylationAnalysisOutcomeTracking06192023.txt",
                          header = T, sep = "\t"))
hlhs_data <- hlhs_data[-30,]
rownames(hlhs_data) <- hlhs_data$sample_name
ph_hlhs <- ph_hlhs[rownames(hlhs_data),]
hlhs_data <- cbind(hlhs_data, ph_hlhs$Sex)

hlhs_data$meldxi[2] <- -6.54133
hlhs_data$meldxi[13] <- -1.34833

patient_gender <- ph_hlhs$Sex
numerical_gender <- rep(0,34)
for (i in 1:34)
{
  if (patient_gender[i] == "F")
  {
    numerical_gender[i] <- 1
  }
}

intercept <- rep(1,34)
design_matrix <- cbind(intercept, hlhs_data$arrhythmia,
                       hlhs_data$age_at_dna_collection_years,
                       numerical_gender)

bmi_design <- cbind(intercept, hlhs_data$bmi,
                       hlhs_data$age_at_dna_collection_years,
                       numerical_gender)

hp_design <- cbind(intercept, hlhs_data$heart_transplant,
                     hlhs_data$age_at_dna_collection_years,
                     numerical_gender)

meld_design <- cbind(intercept, hlhs_data$meldxi,
                   hlhs_data$age_at_dna_collection_years,
                   numerical_gender)
meld_design <- meld_design[-11,]
meld_beta <- meld_beta[,-9]
dim(meld_beta)


colnames(design_matrix) <- c("intercept", "arrhythmia", "age_at_collection", "gender")

devtools::install_github("chrismckennan/CorrConf/CorrConf")
devtools::install_github("chrismckennan/BCconf/BCconf")
library(parallel)

#With all cpgs in Hyperglycemia and MSX
k_bmi <- CorrConf::ChooseK(Y = bmi_beta, Cov = bmi_design, maxK = 20, B = NULL)
output_bmi <- BCconf::Correction(Y = bmi_beta, X = bmi_design, ind.cov = 2, r.confound = k_bmi$K.hat)
hist(output_bmi$p.values, main = "p-values for CpGs in Hyperglycemia and Metabolic Syndrome X")
min(output_bmi$p.values) #0.0001595732
1/dim(bmi_cpg)[1] #0.0001036699

#With all cpgs in Heart Failure
k_hp <- CorrConf::ChooseK(Y = hp_beta, Cov = hp_design, maxK = 20, B = NULL)
output_hp <- BCconf::Correction(Y = hp_beta, X = hp_design, ind.cov = 2, r.confound = k_hp$K.hat)
hist(output_hp$p.values, main = "p-values for CpGs in Heart Failure")
min(output_hp$p.values) #2.697794e-05
1/dim(hp_cpg)[1] #9.00171e-05

#With all cpgs in Liver Cirrosis and Kidney Disease
k_meld <- CorrConf::ChooseK(Y = meld_beta, Cov = meld_design, maxK = 20, B = NULL)
output_meld <- BCconf::Correction(Y = meld_beta, X = meld_design, ind.cov = 2, r.confound = k_meld$K.hat)
hist(output_meld$p.values, main = "p-values for CpGs in Liver Cirrosis and Kidney Disease")
min(output_meld$p.values) #6.359595e-06
1/dim(meld_cpg)[1] #0.0001090156 ?

#With all cpgs in Liver Cirrosis
k_lc <- CorrConf::ChooseK(Y = lc_beta, Cov = meld_design, maxK = 20, B = NULL)
output_lc <- BCconf::Correction(Y = lc_beta, X = meld_design, ind.cov = 2, r.confound = 3)
hist(output_lc$p.values, main = "p-values for CpGs in Liver Cirrosis")
min(output_lc$p.values) #6.359595e-06
1/dim(lc_cpg)[1] #0.0001090156 ?

#With cpgs corr > 0.5
k_bmi_1 <- CorrConf::ChooseK(Y = bmi_beta, Cov = bmi_design, maxK = 20, B = NULL)
output_bmi_1 <- BCconf::Correction(Y = bmi_beta, X = bmi_design, ind.cov = 2, r.confound = k_bmi_1$K.hat)
hist(output_bmi_1$p.values, main = "p-values for CpGs in Hyperglycemia and Metabolic Syndrome X")
min(output_bmi_1$p.values) #0.0005768913
1/dim(bmi_cpg)[1] #0.0007874016

#With cpgs corr > 0.5
k_hp_1 <- CorrConf::ChooseK(Y = hp_beta, Cov = hp_design, maxK = 20, B = NULL)
output_hp_1 <- BCconf::Correction(Y = hp_beta, X = hp_design, ind.cov = 2, r.confound = k_hp_1$K.hat)
hist(output_hp_1$p.values, main = "p-values for CpGs in Heart Failure")
min(output_hp_1$p.values) #0.0005284406
1/dim(hp_cpg)[1] #0.0006963788

#With cpgs corr > 0.5
k_meld_1 <- CorrConf::ChooseK(Y = meld_beta, Cov = meld_design, maxK = 20, B = NULL)
output_meld_1 <- BCconf::Correction(Y = meld_beta, X = meld_design, ind.cov = 2, r.confound = k_meld_1$K.hat)
hist(output_meld_1$p.values, main = "p-values for CpGs in Liver Cirrosis and Kidney Disease")
min(output_meld_1$p.values) #8.994944e-05
1/dim(meld_cpg)[1] #0.0008025682

FNME_beta <- "C:/study/DNAm/beta_hlhs.rds"
FNME_pdat <- "C:/study/DNAm/ph_hlhs.rds"
FNME_anno <- "C:/study/DNAm/35HLHSPatientMethylationAnalysisOutcomeTracking Ivy 06192023.xlsx"

beta <- readRDS(FNME_beta)
pdat <- readRDS(FNME_pdat)
pdat <- janitor::clean_names(pdat)

rownames(pdat) <- pdat$sample_name
se <- SummarizedExperiment( assays = list(beta = beta),
                            rowData = DataFrame(probe_id = rownames(beta)),
                            colData = pdat)

anno  <- readxl::read_xlsx(FNME_anno,
                           sheet = 'sheet_for_cecilia')
anno  <- janitor::clean_names(anno)

#Get genome locations of CpGs
cpg_list <- row.names(beta_hlhs)
head(cpg_list)
subset_annote450k <- annote450k[cpg_list,]
cpg_gene_name <- data.frame(subset_annote450k@listData[["UCSC_RefGene_Name"]])
row.names(cpg_gene_name) <- cpg_list
colnames(cpg_gene_name) <- c('Gene name')


disease_hit_gene <- read.table(file = "disease_hit.txt",
                            header = T, sep = "\t")
H <- read_xlsx("C0020456_disease_gda_summary.xlsx")$Hyperglycemia
LC <- read_xlsx("C0023893_disease_gda_summary.xlsx")$`Liver Cirrhosis, Experimental`
MSX <- read_xlsx("C0524620_disease_gda_summary.xlsx")$`Metabolic Syndrome X`
CHF <- read_xlsx("C0018802_disease_gda_summary.xlsx")$`Congestive heart failure`
CA <- read_xlsx("C0010054_disease_gda_summary.xlsx")$`Coronary Arteriosclerosis`
HF <- read_xlsx("C0018801_disease_gda_summary.xlsx")$`Heart failure`
KD <- read_xlsx("C0022658_disease_gda_summary.xlsx")$`Kidney Diseases`
disease_gene<- c(H,LC,MSX,CHF,CA,HF,KD)
bmi_group <- c(H,MSX)
HF_group <- c(CHF,CA,HF)
meld_group <- c(LC,KD)

hlhs_data$heart_transplant <- as.factor(hlhs_data$heart_transplant)

cpg_subset <- row.names(subset(data.corr, data.corr$Rho>0.5))

#Get a list of CpGs that are in the genes of diseases
disease_cpg <- subset(cpg_gene_name, cpg_gene_name$`Gene name` %in% disease_gene)
subset_beta <- beta_hlhs[row.names(disease_cpg)]

bmi_cpg <- subset(cpg_gene_name, cpg_gene_name$`Gene name` %in% bmi_group)
bmi_beta <- beta_hlhs[row.names(bmi_cpg),]
bmi_cpg <- subset(bmi_cpg, row.names(bmi_cpg) %in% cpg_subset)#1270 CpGs

lc_cpg <- subset(cpg_gene_name, cpg_gene_name$`Gene name` %in% LC)
lc_beta <- beta_hlhs[row.names(lc_cpg),]

hp_cpg <- subset(cpg_gene_name, cpg_gene_name$`Gene name` %in% HF_group)
hp_beta <- beta_hlhs[row.names(hp_cpg),]
hp_cpg <- subset(hp_cpg, row.names(hp_cpg) %in% cpg_subset)#1436 CpGs

meld_cpg <- subset(cpg_gene_name, cpg_gene_name$`Gene name` %in% meld_group)
meld_beta <- beta_hlhs[row.names(meld_cpg),]
meld_cpg <- subset(meld_cpg, row.names(meld_cpg) %in% cpg_subset)#1246 CpGs

cpg_list <- read.csv("cpgs in the disease genes.csv", header = TRUE)
beta_Subset <- beta_hlhs[cpg_list$X,]
