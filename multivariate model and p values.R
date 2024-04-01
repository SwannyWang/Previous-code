setwd("C:/study/DNAm")

library(LaplacesDemon)
library(ggplot2)
library(qvalue)
library(GenomicRanges)



#read files
ph_hlhs <- readRDS('ph_hlhs.rds')
beta_hlhs <- readRDS('beta_hlhs.rds') #Methylation beta value
hlhs_data <- data.frame(read.table
                        (file = "35HLHSPatientMethylationAnalysisOutcomeTracking06192023.txt",
                          header = T, sep = "\t"))

hlhs_data <- hlhs_data[-30,]

#Cell type proportions are the coefficients for linear combinations.
#coefficients for epigenetic clocks?

nk <- ph_hlhs$NK
cd4t <- ph_hlhs$CD4T
cd8t <- ph_hlhs$CD8T
bc <- ph_hlhs$Bcell
Mono <- ph_hlhs$Mono
neu <- ph_hlhs$Neu
ph_hlhs$DNAmAge
#plot x age y dnam age
plot(x = hlhs_data$age_at_dna_collection_years, ph_hlhs$DNAmAge, xlab = "real age", ylab = "DNAm age")
abline(a = 0, b = 1)

#derive acceleration
#take arrhythmia, regress onto dnam age + age + sex, dnam coefficient
#logistic regression
logic <- glm(hlhs_data$arrhythmia~ph_hlhs$DNAmAge+hlhs_data$age_at_dna_collection_years+numerical_sex, family = "binomial")
summary(logic)
#screening method to select cpgs
sum(hlhs_data$arrhythmia)

logic1 <- glm(hlhs_data$arrhythmia~ph_hlhs$DNAmGrimAge+hlhs_data$age_at_dna_collection_years+numerical_sex, family = "binomial")
summary(logic1)

rownames(hlhs_data) <- hlhs_data$sample_name
ph_hlhs <- ph_hlhs[rownames(hlhs_data),]
hlhs_data <- cbind(hlhs_data, ph_hlhs$Sex)

#convert beta values into m values: m = logit(beta)
#What to do with m values?
m_hlhs <- logit(beta_hlhs)
m_hlhs <- m_hlhs[, rownames(ph_hlhs)]
m_hlhs <- data.frame(m_hlhs)
min(beta_hlhs)
max(beta_hlhs)

dim(beta_hlhs)
#Histograms
nsample <- 34
nprobe <- dim(beta_hlhs)[1]

par(mar = c(1, 0.5, 1, 1))
par(mfrow = c(9,2))
for (i in 17:34) {
  hist(beta_hlhs[,i], main = "CpG")
}

#compute mean beta values
mean_beta <- rep(0,nprobe)
for (i in 1:nprobe)
{
  mean_beta[i] <- mean(beta_hlhs[i,])
}
#save functional CpGs
index <- {}
for (i in 1:nprobe)
{
  if (mean_beta[i] > 0.2 && mean_beta[i] < 0.8)
  {
    index <- append(index, i)
  }
}
functional <- beta_hlhs[index,]#225175 CpGs have mean values between 0.2 and 0.8
n_functional <- dim(functional)[1]
n_functional


patient_sex <- ph_hlhs$Sex
numerical_sex <- rep(0,34)
for (i in 1:34)
{
  if (patient_sex[i] == "F")
  {
    numerical_sex[i] <- 1
  }
}


#create design matrix
intercept <- rep(1,34)
design_matrix <- cbind(intercept, hlhs_data$arrhythmia,
                       hlhs_data$age_at_dna_collection_years,
                       numerical_sex)
colnames(design_matrix) <- c("intercept", "arrhythmia", "age_at_collection", "sex")

#SLR
devtools::install_github("chrismckennan/CorrConf/CorrConf")
devtools::install_github("chrismckennan/BCconf/BCconf")
#Using BCconf package
library(parallel)
k_value_all <- CorrConf::ChooseK(Y = beta_hlhs, Cov = design_matrix, maxK = 20, B = NULL)
output_all <- BCconf::Correction(Y = beta_hlhs, X = design_matrix, ind.cov = 2, r.confound = 3)
#histogram of p values
hist(output_all$p.values, main = "distribution of p-values")
#no signal
min(output_all$p.values)
#iid uniform expected value of min p vlaue : 1/700000
#very close
#intermediate
BCconf::Sig.Omega(Corr.Object = output_all)
#cell type proportions correlation with arrhythmia
#reduce dimension

#epigenetic age

#ppt, plots, relevant results, histogram of p values, min p value, no signal

#for functional cpgs
k_value_functional <- CorrConf::ChooseK(Y = functional, Cov = design_matrix, maxK = 20, B = NULL)
output_functional <- BCconf::Correction(Y = functional, X = design_matrix,
                                        ind.cov = 2, r.confound = 4)
#histogram of p values
hist(output_functional$p.values, main = "distribution of p-values for functional CpGs")
#no signal
min(output_functional$p.values)

#separate beta values into 353 horvath cpgs and the others
horvath_cpg <- read.csv("Horvath_cpgs.csv", header = T)
cpg_list <- horvath_cpg$CpGmarker

ind.horvath <- which(rownames(beta_hlhs)%in%cpg_list)
horvath_cpgs_list <- beta_hlhs[ind.horvath,]
#including only 315 cpgs
the_others <- beta_hlhs[-ind.horvath,]
dim(the_others)
#699297 others
length(output_all$p.values)
#get indices?
#expected p value = i/(n+1), n = 699612
n1 <- 315
n2 <- (699612-315)

real_p_values <- output_all$p.values[ind.horvath]
real_p_values <- -log10(real_p_values)
real_p_index <- order(real_p_values)

real_others <- output_all$p.values[-ind.horvath]
real_others_index <- order(real_others)
real_others <- -log10(real_others)

p_values <- real_p_index/(n1+1)
p_values <- -log10(p_values)
p_others <- -log10(real_others_index/(n2+1))
length(real_others)
df1 <- data.frame(cbind(sort(real_p_values), sort(p_values)))
df2 <- data.frame(cbind(sort(real_others), sort(p_others)))


ggplot()+
  geom_point(aes(x = df1$X2, y = df1$X1, col = 'Horvath')) +
  geom_point(aes(x = df2$X2, y = df2$X1, col = 'non-Horvath')) +
  labs(x ="-log10 expected p values", y = "-log10 real p values") +
  scale_color_manual(name="p values for Horvath and non-Horvath CpGs",
                     breaks=c('Horvath', 'non-Horvath'),
                     values=c('Horvath'='red', 'Horvath'='blue'))+
  geom_smooth(method = "lm",
               formula = y ~ x)


ggplot(data = df) + geom_point()

qqnorm(real_p_values)

#Use DNAm data to build epigenetic clocks?
#use existing epigenetic clocks

#Epic array
#353 cpgs #334 out of the 353
#Use Ho

#Restricted to intermediate cpgs mean between 0.2 and 0.8 Multiple testing Q values
#BH procedure Q value estimate pi_0
#package
#min Q value.
#Histogram of p values - uniform - no signal
#spike at 0 - do have signal

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
BiocManager::install("sva")

