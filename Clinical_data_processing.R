library(pheatmap)
library(fifer)
library(FactoMineR)

#I)Loading data 

#A)Cytokine + suppression data

Clinical_data = read.delim("Desktop/COVID_Italy/Clinical_data.txt",dec = ",",na.strings = "")
colnames(Clinical_data) = c("Date","Run","N_sample","Sex","Hospitalization","Severity","APACHE_clinical_score",
                            "CRP","Erythrocyte","Hemoglobin","Platelets","Leukocytes","Neutrophils",
                            "Lymphocytes","Monocytes","Eosinophils","Basophils","Ferritin","Transferrin",
                            "Trombin_clotting_time","Fibrinogen","D_Dimer","Creatinin","eGFR","AST",
                            "ALT","ALP","GGT","CK","LDH","Procalcitonin","pCO2",
                            "pO2","FiO2","F_ratio",
                            "Suppression_by_CD14","Suppression_by_CD14_SN","Suppression_by_CD66_LDN_SN","Suppression_by_CD66_NDN_SN",
                            "BNDF","EGF","EOTAXIN","FGF2","GM-CSF","CXCL1",
                            "HGF","IFNa","IFNg","IL1A","IL1B","IL1RA",
                            "IL2","IL4","IL5","IL6","IL9","IL21",
                            "IL22","IL23","IL27","IL31","TNFB","IL7",
                            "IL8","IL10","IL12p70","IL13","IL15","IL17A",
                            "IL18","CXCL10","LIF","CCL2","CCL3","CCL4",
                            "NGF_beta","PDGF_BB","PIGF-1","CCL5","SCF",
                            "CXCL12","TNFa","VEGF-A","VEGF-D")

x = as.character(Clinical_data$N_sample)
x = strsplit(x,split = " ",fixed = T)
x = unlist(lapply(x,FUN = function(x) {x[1]}))
x = paste("Patient_",x,sep = "")
rownames(Clinical_data) = x

y = as.character(Clinical_data$Severity)
y[y=="mild"] = "Mild"
y[y=="healthy donor"] = "Healthy control"
Clinical_data$Severity = factor(y)

Cytokine_data = Clinical_data[,40:ncol(Clinical_data)]
Suppression_data = Clinical_data[,36:39]
Cellcount_data = Clinical_data[,9:17]
#Removing hemoglobin ... : not a cell type
Cellcount_data = Cellcount_data[,-2]


#Removing not detected cytokines

Non_detected_cytokines = colSums(Cytokine_data==0,na.rm = T)
Cytokine_data = Cytokine_data[,Non_detected_cytokines<10]

#B)ICU data (Only for severe samples)
  
ICU_data = read.delim("Desktop/COVID_Italy/ICU_patient_data.txt",sep="\t",dec = ",",na.strings = "")
ICU_data_cleaned = data.frame(Age = ICU_data$age,
                              Sex = ICU_data$Sex,
                              BMI = ICU_data$BMI,
                              Outcome = ICU_data$Outcome.ICU,
                              SOFA_score = ICU_data$sequential.organ.failure.assessment.score..SOFA.score..at.sampling)

x = as.character(ICU_data$N.SAMPLE)
x = strsplit(x,split = " ",fixed = T)
x = unlist(lapply(x,FUN = function(x) {x[1]}))
x = paste("Patient_",x,sep = "")
rownames(ICU_data_cleaned) = x

#II)Looking for the most interesting variables to distinguish the different class of patients

#A)Difference between HC/Mild/Severes

#1)Computing the P-values

P_value_severity = function(x) {
  m = aov(x~Clinical_data$Severity)
  m = anova(m)
  m = m$`Pr(>F)`[1]
  return(m)
}

Pvalue_cytokine_severity = sapply(log(Cytokine_data+1),FUN=P_value_severity)
Pvalue_cytokine_severity = p.adjust(Pvalue_cytokine_severity,method = "bonf")
Pvalue_cytokine_severity = Pvalue_cytokine_severity[order(Pvalue_cytokine_severity)]
Pvalue_cytokine_severity = -log10(Pvalue_cytokine_severity)


Pvalue_cell_count_severity = sapply(log10(1+Cellcount_data),FUN=P_value_severity)
Pvalue_cell_count_severity = p.adjust(Pvalue_cell_count_severity,method = "bonf")
Pvalue_cell_count_severity = Pvalue_cell_count_severity[order(Pvalue_cell_count_severity)]
Pvalue_cell_count_severity = -log10(Pvalue_cell_count_severity)

pdf("Desktop/COVID_Italy/Figures/Significance_blood_feature.pdf",width = 5.5,height = 5)
par(las=1,bty="l")
barplot(Pvalue_cytokine_severity,ylim=c(0,4),names.arg = "",ylab="-Log10(P-value)",
        cex.lab=1.3,main="Cytokine p-value",cex.axis = 1.3)
abline(h=2,lwd=1.5,lty=2,col="red")
barplot(Pvalue_cell_count_severity,ylim=c(0,4),names.arg = "",ylab="-Log10(P-value)",
        cex.lab=1.3,main="Blood cell count p-value",cex.axis = 1.3)
abline(h=2,lwd=1.5,lty=2,col="red")
dev.off()

#2)Plot of the interesting cytokines/cell counts

pdf("Desktop/COVID_Italy/Figures/Cytokine_prism_plots.pdf",width = 4,height = 5.5,useDingbats = F)
par(las=1,bty="l")
prism.plots(`IL6`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
                bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
                col="black",cex=1.5,ylab="IL6 (pg/mL)",cex.lab=1.3,cex.axis=1.3)
prism.plots(`IL1RA`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
            bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
            col="black",cex=1.5,ylab="IL1RA (pg/mL)",cex.lab=1.3,cex.axis=1.3)
prism.plots(`VEGF-A`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
            bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
            col="black",cex=1.5,ylab="VEGF-A (pg/mL)",cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Desktop/COVID_Italy/Figures/Cell_count_prism_plots.pdf",width = 4,height = 5.5,useDingbats = F)
par(las=1,bty="l")
prism.plots(`Erythrocyte`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
            bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
            col="black",cex=1.5,ylab="Erythrocytes (10^12/L)",cex.lab=1.3,cex.axis=1.3)
prism.plots(`Neutrophils`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
            bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
            col="black",cex=1.5,ylab="Neutrophils (10^9/L)",cex.lab=1.3,cex.axis=1.3)
prism.plots(`Lymphocytes`~Severity,Clinical_data[order(Clinical_data$Severity),],yaxs="i",pch=21,
            bg=string.to.colors(Clinical_data$Severity[order(Clinical_data$Severity)],colors = c("skyblue2","orange","brown3")),
            col="black",cex=1.5,ylab="Lymphocytes (10^9/L)",cex.lab=1.3,cex.axis=1.3)
dev.off()


#B)Difference between survival/death in the ICU

P_value_survival = function(x) {
  m = aov(x~ICU_data_cleaned$Outcome)
  m = anova(m)
  m = m$`Pr(>F)`[1]
  return(m)
}

Pvalue_cytokine_outcome = sapply(sqrt(Cytokine_data[rownames(ICU_data_cleaned),]),FUN=P_value_survival)
Pvalue_cell_count_outcome = sapply(log10(1+Cellcount_data[rownames(ICU_data_cleaned),]),FUN=P_value_survival)
Pvalue_suppression_outcome = sapply(log10(1+Suppression_data[rownames(ICU_data_cleaned),]),FUN=P_value_survival)











