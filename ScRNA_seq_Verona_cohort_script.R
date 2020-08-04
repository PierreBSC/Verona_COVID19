#0) Loading packages and basic functions

library(RColorBrewer)
library(pagoda2)
library(uwot)
library(Matrix)
library(MASS)
library(igraph)
library(pheatmap)
library(FactoMineR)
library(leiden)
library(fifer)
library(Matrix)

rotate <- function(x) t(apply(x, 2, rev))

logit = function(x) {
  y = log(x/(1-x))
}

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}


color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

color_convertion_correlation=function(x=NULL) {
  f <- colorRamp(c("darkblue","cornflowerblue","white","yellow","red"))
  x=as.numeric(x)
  x = (x+1)/2
  x_color=f(x)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}


#I)Data loading and filtering

#A)Loading 

Data_blood_1 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_1/Sangue_1/", n.cores = 1, verbose = T)
Data_BAL_1 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_1/BAL1//", n.cores = 1, verbose = T)

Data_blood_2 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_2/Sangue_2/", n.cores = 1, verbose = T)
Data_BAL_2 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_2/BAL2//", n.cores = 1, verbose = T)

Data_blood_3 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_3/Sangue_3/", n.cores = 1, verbose = T)
Data_BAL_3 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_3/BAL3//", n.cores = 1, verbose = T)

Data_blood_4 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_4/Sangue_4/", n.cores = 1, verbose = T)
Data_BAL_4 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_4/BAL4//", n.cores = 1, verbose = T)

Data_blood_5 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_5/Sangue_5/", n.cores = 1, verbose = T)
Data_BAL_5 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_5/BAL5//", n.cores = 1, verbose = T)

Data_blood_6 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_6/Sangue_6/", n.cores = 1, verbose = T)
Data_BAL_6 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_6/BAL6//", n.cores = 1, verbose = T)

Data_blood_7 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_7/Sangue_7/", n.cores = 1, verbose = T)
Data_BAL_7 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_7/BAL7//", n.cores = 1, verbose = T)

Data_blood_8 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_8/Sangue_8/", n.cores = 1, verbose = T)
Data_BAL_8 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_8/BAL8///", n.cores = 1, verbose = T)

Data_blood_9 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_9/Sangue9//", n.cores = 1, verbose = T)
Data_BAL_9 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_9/BAL9//", n.cores = 1, verbose = T)

Data_blood_10 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_10/Sangue_10/", n.cores = 1, verbose = T)
Data_BAL_10 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_10/BAL10///", n.cores = 1, verbose = T)

Data_blood_11 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_11/Sangue_11/", n.cores = 1, verbose = T)
Data_BAL_11 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_11/BAL11///", n.cores = 1, verbose = T)

Data_blood_14 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_14/Sangue14//", n.cores = 1, verbose = T)
Data_BAL_14 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_14/BAL14/", n.cores = 1, verbose = T)

Data_blood_16 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_16/Sangue16///", n.cores = 1, verbose = T)
Data_BAL_16 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_16/BAL16/", n.cores = 1, verbose = T)

Data_blood_17 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_17/Sangue17/", n.cores = 1, verbose = T)
Data_BAL_17 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_17/BAL17/", n.cores = 1, verbose = T)

Data_blood_18 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_18/Sangue18/", n.cores = 1, verbose = T)
Data_BAL_18 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_18/BAL18//", n.cores = 1, verbose = T)

Data_blood_22 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_22/Sangue22//", n.cores = 1, verbose = T)
Data_BAL_22 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_22/BAL22//", n.cores = 1, verbose = T)

Data_blood_23 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_23/Sangue23//", n.cores = 1, verbose = T)
Data_BAL_23 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_23/BAL23///", n.cores = 1, verbose = T)

Data_blood_24 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_24/Sangue24//", n.cores = 1, verbose = T)
Data_BAL_24 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_24/BAL24//", n.cores = 1, verbose = T)

Data_blood_25 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_25/Sangue25/", n.cores = 1, verbose = T)
Data_BAL_25 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_25/BAL25/", n.cores = 1, verbose = T)

Data_blood_27 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_27/Sangue27/", n.cores = 1, verbose = T)
Data_BAL_27 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_27/BAL27/", n.cores = 1, verbose = T)

Data_blood_29 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_29/Sangue29/", n.cores = 1, verbose = T)
Data_BAL_29 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_29/BAL29/", n.cores = 1, verbose = T)

Data_blood_30 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_30//Sangue30//", n.cores = 1, verbose = T)
Data_BAL_30 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_30/BAL30/", n.cores = 1, verbose = T)


Data_blood_32 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_32//Sangue_32//", n.cores = 1, verbose = T)
Data_blood_33 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_33/Sangue_33//", n.cores = 1, verbose = T)
Data_blood_34 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_34/Sangue_34/", n.cores = 1, verbose = T)
Data_blood_35 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_35/Sangue_35/", n.cores = 1, verbose = T)
Data_blood_36 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_36/Sangue_36/", n.cores = 1, verbose = T)
Data_blood_37 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_37/Sangue_37//", n.cores = 1, verbose = T)
Data_blood_38 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_38/Sangue_38//", n.cores = 1, verbose = T)
Data_blood_39 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_39/Sangue_39//", n.cores = 1, verbose = T)
Data_blood_40 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_40/Sangue_40/", n.cores = 1, verbose = T)
Data_blood_41 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_41/Sangue41//", n.cores = 1, verbose = T)
Data_blood_42 = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_42/Sangue42//", n.cores = 1, verbose = T)


Data_merged = cbind(Data_blood_1,Data_blood_2,Data_blood_3,Data_blood_4,Data_blood_5,
                    Data_blood_6,Data_blood_7,Data_blood_8,Data_blood_9,Data_blood_10,
                    Data_blood_11,Data_blood_14,Data_blood_16,Data_blood_17,Data_blood_18,Data_blood_22,
                    Data_blood_23,Data_blood_24,Data_blood_25,Data_blood_27,Data_blood_29,
                    Data_blood_30,Data_blood_32,Data_blood_33,Data_blood_34,Data_blood_35,Data_blood_36,
                    Data_blood_37,Data_blood_38,Data_blood_39,Data_blood_40,Data_blood_41,Data_blood_42,
                    Data_BAL_1,Data_BAL_2,Data_BAL_3,Data_BAL_4,Data_BAL_5,
                    Data_BAL_6,Data_BAL_7,Data_BAL_8,Data_BAL_9,Data_BAL_10,
                    Data_BAL_11,Data_BAL_14,Data_BAL_16,Data_BAL_17,Data_BAL_18,Data_BAL_22,
                    Data_BAL_23,Data_BAL_24,Data_BAL_25,Data_BAL_27,Data_BAL_29,Data_BAL_30)

Merged_annotation = c(rep("Blood_1",ncol(Data_blood_1)),rep("Blood_2",ncol(Data_blood_2)),rep("Blood_3",ncol(Data_blood_3)),rep("Blood_4",ncol(Data_blood_4)),rep("Blood_5",ncol(Data_blood_5)),
                      rep("Blood_6",ncol(Data_blood_6)),rep("Blood_7",ncol(Data_blood_7)),rep("Blood_8",ncol(Data_blood_8)),rep("Blood_9",ncol(Data_blood_9)),rep("Blood_10",ncol(Data_blood_10)),
                      rep("Blood_11",ncol(Data_blood_11)),rep("Blood_14",ncol(Data_blood_14)),rep("Blood_16",ncol(Data_blood_16)),rep("Blood_17",ncol(Data_blood_17)),rep("Blood_18",ncol(Data_blood_18)),
                      rep("Blood_22",ncol(Data_blood_22)),rep("Blood_23",ncol(Data_blood_23)),rep("Blood_24",ncol(Data_blood_24)),rep("Blood_25",ncol(Data_blood_25)),rep("Blood_27",ncol(Data_blood_27)),
                      rep("Blood_29",ncol(Data_blood_29)),rep("Blood_30",ncol(Data_blood_30)),rep("Blood_32",ncol(Data_blood_32)),rep("Blood_33",ncol(Data_blood_33)),rep("Blood_34",ncol(Data_blood_34)),
                      rep("Blood_35",ncol(Data_blood_35)),rep("Blood_36",ncol(Data_blood_36)),rep("Blood_37",ncol(Data_blood_37)),rep("Blood_38",ncol(Data_blood_38)),rep("Blood_39",ncol(Data_blood_39)),
                      rep("Blood_40",ncol(Data_blood_40)),rep("Blood_41",ncol(Data_blood_41)),rep("Blood_42",ncol(Data_blood_42)),
                      rep("BAL_1",ncol(Data_BAL_1)),rep("BAL_2",ncol(Data_BAL_2)),rep("BAL_3",ncol(Data_BAL_3)),rep("BAL_4",ncol(Data_BAL_4)),rep("BAL_5",ncol(Data_BAL_5)),
                      rep("BAL_6",ncol(Data_BAL_6)),rep("BAL_7",ncol(Data_BAL_7)),rep("BAL_8",ncol(Data_BAL_8)),rep("BAL_9",ncol(Data_BAL_9)),rep("BAL_10",ncol(Data_BAL_10)),
                      rep("BAL_11",ncol(Data_BAL_11)),rep("BAL_14",ncol(Data_BAL_14)),rep("BAL_16",ncol(Data_BAL_16)),rep("BAL_17",ncol(Data_BAL_17)),rep("BAL_18",ncol(Data_BAL_18)),
                      rep("BAL_22",ncol(Data_BAL_22)),rep("BAL_23",ncol(Data_BAL_23)),rep("BAL_24",ncol(Data_BAL_24)),rep("BAL_25",ncol(Data_BAL_25)),rep("BAL_27",ncol(Data_BAL_27)), 
                      rep("BAL_29",ncol(Data_BAL_29)),rep("BAL_30",ncol(Data_BAL_30)))
X = strsplit(x = Merged_annotation,split = "_")
Tissue = unlist(lapply(X,FUN = function(x) {x[1]}))
Patient = unlist(lapply(X,FUN = function(x) {x[2]}))
Patient = paste("Patient_",Patient,sep = "")

Type_sample = c(rep("Severe",ncol(Data_blood_1)),rep("Severe",ncol(Data_blood_2)),rep("Severe",ncol(Data_blood_3)),rep("Severe",ncol(Data_blood_4)),rep("Severe",ncol(Data_blood_5)),
                rep("Severe",ncol(Data_blood_6)),rep("Severe",ncol(Data_blood_7)),rep("Severe",ncol(Data_blood_8)),rep("Severe",ncol(Data_blood_10)),rep("Severe",ncol(Data_blood_11)),
                rep("Mild",ncol(Data_blood_32)),rep("Mild",ncol(Data_blood_33)),rep("Mild",ncol(Data_blood_34)),rep("Mild",ncol(Data_blood_35)),rep("Mild",ncol(Data_blood_36)),
                rep("Mild",ncol(Data_blood_37)),rep("Control",ncol(Data_blood_38)),rep("Control",ncol(Data_blood_39)),rep("Control",ncol(Data_blood_40)),
                rep("Severe",ncol(Data_BAL_1)),rep("Severe",ncol(Data_BAL_2)),rep("Severe",ncol(Data_BAL_3)),rep("Severe",ncol(Data_BAL_4)),rep("Severe",ncol(Data_BAL_5)),
                rep("Severe",ncol(Data_BAL_6)),rep("Severe",ncol(Data_BAL_7)),rep("Severe",ncol(Data_BAL_8)),rep("Severe",ncol(Data_BAL_10)),rep("Severe",ncol(Data_BAL_11)))

Patient_annotation_table = read.delim("Desktop/COVID_Italy/Patient_annotation.txt")
rownames(Patient_annotation_table) = Patient_annotation_table$Patient

Patient_annotation_table$Patient = paste("Patient_",as.character(Patient_annotation_table$Patient),sep = "")
x = Patient_annotation_table
y = Patient_annotation_table
rownames(x) = paste("BAL_",rownames(x),sep = "")
rownames(y) = paste("Blood_",rownames(y),sep = "")
Patient_annotation_table_extended = rbind(x,y)
Patient_annotation_table_extended$Tissue = ifelse(grepl(rownames(Patient_annotation_table_extended),pattern = "Blood"),yes ="Blood",no="BAL")
Patient_annotation_table_extended$Type = as.character(Patient_annotation_table_extended$Type)
Patient_annotation_table_extended$Type[Patient_annotation_table_extended$Type=="home"] = "Healthy control"
Patient_annotation_table_extended$Type[Patient_annotation_table_extended$Type=="ICU"] = "Severe"
Patient_annotation_table_extended$Type[Patient_annotation_table_extended$Type=="Infective Disease"] = "Mild"


#B)Filtering

Lib_size = colSums(Data_merged)
Gene_size = rowSums(Data_merged)

pdf("Desktop/COVID_Italy/Figures/QC_plots.pdf",width = 6.5,height = 6.5)
par(las=1,bty="l")
hist(log10(1+Gene_size),40,main="Total gene UMIs distribution",xlab="Total UMIs (Log10)",
     xaxs="i",col="grey",yaxs='i')
abline(v=log10(50),lty=2,lwd=2,col="red")
hist(log10(1+Lib_size),40,main="Total cellular UMIs distribution",xlab="Total UMIs (Log10)",
     xaxs="i",col="grey",yaxs='i',xlim=c(2,6),ylim=c(0,8000))
abline(v=log10(500),lty=2,lwd=2,col="red")

boxplot(log10(Lib_size)~Merged_annotation,outline=F)

Number_gene_detected = Matrix::colSums(Data_merged>0)

MT_genes = rownames(Data_merged)
MT_genes = MT_genes[grepl(pattern = "MT-",MT_genes)]
MT_proportion = colSums(Data_merged[MT_genes,])/colSums(Data_merged)*100
hist(MT_proportion,40,main="",xlab="Proportion of MT genes (%)",
     xaxs="i",col="grey",yaxs='i')
abline(v=20,lty=2,lwd=2,col="red")
dev.off()

data_count = Data_merged[Gene_size >50 ,MT_proportion < 20 & Lib_size > 500]

plot(Lib_size,Number_gene_detected,log="xy",pch=21,bg=color_convertion(MT_proportion))
plot(Lib_size,Number_gene_detected,log="xy",pch=21,bg=string.to.colors(Tissue))

par(las=1)
boxplot(log10(Lib_size)~Merged_annotation,col=string.to.colors(string = grepl(pattern = "BAL",levels(factor(Merged_annotation))),c('firebrick2',"dodgerblue2")),
        outline=F,xlab="Sample",ylab="Total cell UMIs (Log10)",cex.lab=1.2,cex.axis=1.2)
boxplot(MT_proportion~(Merged_annotation),
        col=string.to.colors(string = grepl(pattern = "BAL",levels(factor(Merged_annotation))),c('firebrick2',"dodgerblue2")),
        outline=F,xlab="Sample",ylab="Proportion of MT genes (%)",cex.lab=1.2)

#C)Gene selection

Zero_proportion = rowSums(data_count==1)/ncol(data_count)
Zero_proportion = 1 -Zero_proportion
Mean_expression = rowMeans(data_count)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes = names(Zeros_excess_proportion[1:4000])

Transformed_zero_proportion = log(Zero_proportion/(1-Zero_proportion))
plot(log10(Mean_expression),Transformed_zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Logit_excess = lm(Transformed_zero_proportion~log10(Mean_expression))
abline(coef(Logit_excess),lwd=2,lty=2,col="red")

#II)Clinical properties of the cohort

#A)Global properties of the cohort

ICU_data = read.delim("Desktop/COVID_Italy/ICU_patient_data.txt",sep="\t",dec = ",",na.strings = "")
ICU_data$Mock_variable = " "

pdf("Desktop/COVID_Italy/Figures/Death_sex_ratio.pdf",width = 5.5,height = 8)
par(las=1,mfrow = c(1,2))
Death_ratio = table(ICU_data$Outcome.ICU)/sum(table(ICU_data$Outcome.ICU))
barplot(matrix(Death_ratio)*100,cex.axis = 1.3,col=c("grey","red3"),
        ylab = "Death ratio (%)",cex.lab=1.3)

Sex_ratio = table(ICU_data$Sex)/sum(table(ICU_data$Sex))
barplot(matrix(Sex_ratio)*100,cex.axis = 1.3,col=c("skyblue2","orange"),
        ylab = "Female/male ratio (%)",cex.lab=1.3)
dev.off()

Smoker_ratio = table(ICU_data$Smoker)/sum(table(ICU_data$Smoker))
par(las=1,mar=c(5,8,2,2))
barplot(matrix(Smoker_ratio)*100,xlim=c(0,3),cex.axis = 1.3,
        col=brewer.pal(3, "Spectral"), ylab = "Proportion of smokers (%)",cex.lab=1.3)
legend("right",fill = brewer.pal(3, "Spectral"),legend = c("No","Former","Yes"),bty='n',cex=1.3)

dev.off()

pdf("Desktop/COVID_Italy/Figures/Age_BMI.pdf",width = 6.5,height = 6.5)
par(las=1,mfrow = c(1,2))
prism.plots(BMI~Mock_variable,ICU_data,bty="l",cex.axis=1.3,
            pch=21,bg="orange",col="black",cex=2,cex.lab=1.3)
prism.plots(age~Mock_variable,ICU_data,bty="l",cex.axis=1.3,ylab="Age (year)",
            pch=21,bg="orange",col="black",cex=2,cex.lab=1.3)
dev.off()

prism.plots(age~Outcome.ICU,ICU_data,bty="l",cex.axis=1.3,ylab="Age (year)",
            pch=21,bg="orange",col="black",cex=2,cex.lab=1.3)


#B)Clinical properties/state of the patient 


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

#C)Cleaning the table so we can perfom joined analysis

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
y = ICU_data_cleaned[Patient_annotation_table_extended$Patient,]

z = Clinical_data[Patient_annotation_table_extended$Patient,]

Cytokine_data_blood_data = z
rownames(Cytokine_data_blood_data) = rownames(Patient_annotation_table_extended)

Patient_annotation_table_extended$BMI = y$BMI
Patient_annotation_table_extended$Age = y$Age
Patient_annotation_table_extended$Outcome = y$Outcome
Patient_annotation_table_extended$SOFA_score = y$SOFA_score
Patient_annotation_table_extended$APACHE_score = z$APACHE_clinical_score
Patient_annotation_table_extended$F_ratio = z$F_ratio

#III)Single-cell analysis per se 

#0) Filtering of the low quality cells

Patient_count = Patient[MT_proportion < 20 & Lib_size > 500]
Tissue_count = Tissue[MT_proportion < 20 & Lib_size > 500]
Merged_annotation_count = Merged_annotation[MT_proportion < 20 & Lib_size > 500]
Lib_size_count = Lib_size[MT_proportion < 20 & Lib_size > 500]

Severity_count = Patient_annotation_table_extended[Merged_annotation_count,"Type"]
names(Severity_count) = colnames(data_count)
Severity_count[is.na(Severity_count)] = Severe

names(Tissue_count) = colnames(data_count)
names(Merged_annotation_count) = colnames(data_count)
names(Lib_size_count) = colnames(data_count)

Number_cells = table(Merged_annotation_count)
Number_cells = data.frame(N = log10(as.numeric(Number_cells)), 
                          Annotation = paste(Patient_annotation_table_extended[names(Number_cells),"Tissue"],Patient_annotation_table_extended[names(Number_cells),"Type"]))
Number_cells = Number_cells[complete.cases(Number_cells),]
Number_cells = Number_cells[Number_cells$Annotation!="NA NA",]

pdf("Desktop/COVID_Italy/Figures/Number_cell_sequenced.pdf",width = 8,height = 6)
par(las=1,bty="l")
prism.plots(N~Annotation,Number_cells,pch=21,cex=1.4,bg="grey",col="black",ylim=c(2.5,4),
            yaxs="i",xaxs="i",ylab="Number of high-quality cells (Log10)",cex.lab=1.3)
dev.off()

#A)Pagoda2 analysis

Genes_to_remove = table(rownames(data_count))
Genes_to_remove = names(which(Genes_to_remove>1))
data_count = data_count[!rownames(data_count)%in%Genes_to_remove,]
colnames(data_count) = paste("Cell",1:ncol(data_count))

Selected_genes = intersect(Selected_genes,rownames(data_count))

r <- Pagoda2$new(data_count,log.scale=F)
r$adjustVariance(plot=T,gam.k=5)
r$calculatePcaReduction(nPcs=150,n.odgenes=0.5e3,odgenes = Selected_genes,)
r$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='PCA')

Leiden_clustering = leiden(r$graphs$PCA)
names(Leiden_clustering) = colnames(data_count)
r$clusters$PCA$Leiden = factor(Leiden_clustering)

r$getDifferentialGenes(type = "PCA",clusterType = "community",verbose = T,z.threshold = 3)
r$getDifferentialGenes(type = "PCA",clusterType = "Leiden",verbose = T,z.threshold = 3)

umap_plot = umap(r$reductions$PCA,n_neighbors = 30,spread = 3,
                 n_components = 2,metric = "cosine",verbose = T)

#Creating a nice color palette
N_cluster = length(unique(r$clusters$PCA$Leiden))
optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster)/N_cluster)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

#Various UMAP plot to control the quality of the data
plot(umap_plot,pch=21,bg=string.to.colors(r$clusters$PCA$Leiden,col=optimal_palette),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot,pch=21,bg=string.to.colors(Patient_count),main="Patient distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot,pch=21,bg=string.to.colors(Tissue_count,colors = c('firebrick2',"dodgerblue2")),main="Patient distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot,pch=21,bg=color_convertion(r$counts[,"IL6"]),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#B)Aggregation of the clusters (Figure S1E)

Mean_expression_cluster = aggregate(as.matrix(r$counts[,Selected_genes]),FUN = mean,by=list(r$clusters$PCA$Leiden))
rownames(Mean_expression_cluster) =Mean_expression_cluster$Group.1
Mean_expression_cluster = t(Mean_expression_cluster[,-1])
Meta_clustering = pheatmap(cor(Mean_expression_cluster,method = "spearman"),clustering_method = "ward.D2")
Meta_clustering = Meta_clustering$tree_col
Order_cluster = Meta_clustering$order

Cluster_distribution = table(Merged_annotation_count,r$clusters$PCA$Leiden)

#Low quality clusters (blood cells, dead cells, Platelets) need to be removed...
Cluster_to_remove = c(3,9,15,23,24) #Those corresponds to cells, dead cells, Platelets....
Cluster_distribution = Cluster_distribution[,-Cluster_to_remove]

Cluster_distribution_normalized = Cluster_distribution/rowSums(Cluster_distribution)
barplot(t(Cluster_distribution_normalized))
pheatmap(Cluster_distribution_normalized,clustering_method = "ward.D",
         annotation_row = data.frame(Type = Patient_annotation_table_extended[rownames(Cluster_distribution_normalized),"Type"],
                                     Tissue = Patient_annotation_table_extended[rownames(Cluster_distribution_normalized),"Tissue"],
                                     row.names = rownames(Cluster_distribution_normalized)))

#C)Isolating the main cellular compartments

#Lymphoid 

Lymphoid_cluster = c(21,10,11,18,12,14,25,16,19)
Lymphoid_cells = r$clusters$PCA$Leiden%in%Lymphoid_cluster

plot(umap_plot,pch=21,bg=string.to.colors(r$clusters$PCA$Leiden%in%Lymphoid_cluster),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#Neutrophil 

Neutrophil_cluster = c(20,15,4,1,7,2,8)
Neutrophil_cells = r$clusters$PCA$Leiden%in%Neutrophil_cluster

plot(umap_plot,pch=21,bg=string.to.colors(Neutrophil_cells),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#Monocytes/Macrophage 

Monocyte_macrophage_clusters = c(6,22,5,13)
Monocyte_macrophage_cells = r$clusters$PCA$Leiden%in%Monocyte_macrophage_clusters

plot(umap_plot,pch=21,bg=string.to.colors(Monocyte_macrophage_cells),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")

#Epithelial 

Epithelial_cluster = c(17,26)
Epithelial_cells = r$clusters$PCA$Leiden%in%Epithelial_cluster

Major_cell_type = rep("Unassigned",ncol(data_count))
Major_cell_type[Neutrophil_cells] = "Neutrophils"
Major_cell_type[Lymphoid_cells] = "Lymphoid cells"
Major_cell_type[Epithelial_cells] = "Epithelial cells"
Major_cell_type[Monocyte_macrophage_cells] = "Monocytes/Macrophages"


pdf("Desktop/COVID_Italy/Figures/Rough_clustering.pdf",width = 7,height = 7)
plot(umap_plot[Major_cell_type!="Unassigned",],pch=21,cex=0.5,lwd=0.5,
     bg=string.to.colors(Major_cell_type[Major_cell_type!="Unassigned"],colors = c("firebrick3","orange","dodgerblue3","green3")),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()

pdf("Desktop/COVID_Italy/Figures/Tissue_UMAP.pdf",width = 7,height = 7)
plot(umap_plot[Major_cell_type!="Unassigned" & Tissue_count=="Blood",],pch=21,cex=0.5,lwd=0.5,
     bg ="firebrick3",main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot[Major_cell_type!="Unassigned" & Tissue_count!="Blood",],pch=21,cex=0.5,lwd=0.5,
     bg ="dodgerblue",main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()

pdf("Desktop/COVID_Italy/Figures/Heatmap_metacluster.pdf",width = 6,height = 6)
Meta_clustering = pheatmap(cor(Mean_expression_cluster,method = "spearman")[Cluster_total_cells,Cluster_total_cells],
                           cluster_cols = F,cluster_rows = F,show_colnames = F,show_rownames = F,breaks = base::seq(0,1,length.out=100))
dev.off()

x = table(grepl(Merged_annotation_count,pattern = "BAL"),Major_cell_type)
x = x[,-5]
x = x/rowSums(x)
par(las=1,mar=c(6,6,4,4))
barplot(t(x),col=c("green3","orange","dodgerblue3","firebrick3"),names.arg = c("Blood","BAL"),
        cex.lab=1.7,ylab="Cell type distribution",cex.axis=1.3,cex.names=1.3)

#D)Generation of the heatmap (Figure 1C)

V = table(r$clusters$PCA$Leiden)
V = V/sum(V)*100
Cluster_total_cells = names(which(V>1))
Cluster_total_cells = Order_cluster[Order_cluster%in%Cluster_total_cells]
Cluster_total_cells = Cluster_total_cells[!Cluster_total_cells%in%Cluster_to_remove]

List_marker_genes = c()
for (k in Cluster_total_cells) {
  X = r$diffgenes$PCA$Leiden[[k]]
  X = X[X$highest,]
  X = X[X$M>1 & X$Z>10,]
  X = X[order(X$Z,decreasing = T),]
  List_marker_genes = c(List_marker_genes,rownames(X[1:5,]))
}

List_marker_genes = List_marker_genes[!is.na(List_marker_genes)]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NA")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MALAT1")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MT-")]
List_marker_genes = unique(List_marker_genes)

TPM_total= log2(1+t(data_count[List_marker_genes,])/colSums(data_count)*10^6)
Heatmap_data = c()

for (k in Cluster_total_cells) {
  X = TPM_total[r$clusters$PCA$Leiden==k,List_marker_genes]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Heatmap_data = cbind(Heatmap_data,X)
}


Heatmap_data=as.matrix(Heatmap_data)
lim_values=quantile(as.numeric((Heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Heatmap_data[Heatmap_data<lim_values[1]]=lim_values[1]
Heatmap_data[Heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r$clusters$PCA$Leiden)
png("Desktop/COVID_Italy/Figures//Heatmap_total.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(Heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_total_cells)])/sum(size_clusters[as.character(Cluster_total_cells)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Heatmap_data[,]))), 
     labels= ( rownames(Heatmap_data)[nrow(Heatmap_data):1]),
     las= 2, cex.axis=4.5, tick = F )
dev.off()

png("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Tissue_color_bar_total_2.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(400,40,4,4))
plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
x_position = seq(0,1,length.out = ncol(Heatmap_data))
abline(v=x_position,col=string.to.colors(Tissue_count[colnames(Heatmap_data)],colors = c("black","white")))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_total_cells)])/sum(size_clusters[as.character(Cluster_total_cells)]),lwd=10,lty=2)
dev.off()

png("Desktop/COVID_Italy/Figures/Severity_color_bar_total.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(400,40,4,4))
plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
x_position = seq(0,1,length.out = ncol(Heatmap_data))
abline(v=x_position,col=string.to.colors(Severity_count[colnames(Heatmap_data)],colors = c("brown3","skyblue2","orange")))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_total_cells)])/sum(size_clusters[as.character(Cluster_total_cells)]),lwd=10,lty=2)
dev.off()

N_cluster_high_quality = length(Cluster_total_cells)

optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster_high_quality)/N_cluster_high_quality)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

#Generation of Figure 1B
pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/UMAP_total_clusters.pdf",width = 12,height = 12)
plot(umap_plot[r$clusters$PCA$Leiden%in%Cluster_total_cells,],pch=21
     ,bg=string.to.colors(r$clusters$PCA$Leiden[r$clusters$PCA$Leiden%in%Cluster_total_cells],col=optimal_palette)
     ,main="Cluster distribution",xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()

umap_temp = data.frame(umap_plot[r$clusters$PCA$Leiden%in%Cluster_total_cells,])
umap_temp$tissue = Tissue_count[r$clusters$PCA$Leiden%in%Cluster_total_cells]

X = string.to.colors(r$clusters$PCA$Leiden[r$clusters$PCA$Leiden%in%Cluster_total_cells],optimal_palette)
names(X) = r$clusters$PCA$Leiden[r$clusters$PCA$Leiden%in%Cluster_total_cells]
X = X[as.character(Cluster_total_cells)]

Bar_value = size_clusters[Cluster_total_cells]
Bar_value = Bar_value/sum(Bar_value)

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures//Barplot_cluster_total_cell.pdf",width = 8,height = 2)
barplot(matrix(Bar_value),horiz = T,col=X,xaxt="n")
dev.off()

#E)General analysis for the main cell types

#1)For the blood : composition (Figure 1E)

Population_distribution_Blood = table(Merged_annotation_count[Tissue_count=="Blood"],Major_cell_type[Tissue_count=="Blood"])
Population_distribution_Blood = Population_distribution_Blood[,c(-1,-5)]
Population_distribution_Blood_normalized = Population_distribution_Blood/rowSums(Population_distribution_Blood)

X = Patient_annotation_table_extended[rownames(Cluster_distribution_Blood),"Type"]
X[is.na(X)] = "Severe"
Severity_factor = X

Proportion_change = function(x) {
  m = aov(x~Severity_factor)
  m = anova(m)
  return(m$`Pr(>F)`[1])
}
x = apply(Population_distribution_Blood_normalized,MARGIN = 2,FUN =Proportion_change )

U = data.frame(Lymphoid = Population_distribution_Blood_normalized[,1]*100,
               Mon_Mac = Population_distribution_Blood_normalized[,2]*100,
               Neutrophils = Population_distribution_Blood_normalized[,3]*100,Severity = Severity_factor )
U = U[order(U$Severity),]

pdf("Desktop/COVID_Italy/Figures/Proportion_total_blood.pdf",width = 4,height = 6.5,useDingbats = F)
par(las=1,bty="l")
prism.plots(Neutrophils~Severity,U,pch=21,cex=2,col="black",ylim=c(0,100),
            bg=string.to.colors(U$Severity,colors = c("skyblue2","orange","brown3")))
prism.plots(Lymphoid~Severity,U,pch=21,cex=2,col="black",ylim=c(0,100),
            bg=string.to.colors(U$Severity,colors = c("skyblue2","orange","brown3")))
prism.plots(Mon_Mac~Severity,U,pch=21,cex=2,col="black",ylim=c(0,100),
            bg=string.to.colors(U$Severity,colors = c("skyblue2","orange","brown3")))
dev.off()

#2)For the BAL : production of cytokine by the inflammatory macrophage compartment (Figure 1F)

Proportion_Inflamatory_macrophages = table(Merged_annotation_count[Tissue_count=="BAL"],r$clusters$PCA$Leiden[Tissue_count=="BAL"]==5)
Proportion_Inflamatory_macrophages = Proportion_Inflamatory_macrophages/rowSums(Proportion_Inflamatory_macrophages)
Proportion_Inflamatory_macrophages = Proportion_Inflamatory_macrophages[,2]

U = data.frame(Proportion = Proportion_Inflamatory_macrophages*100,Outcome = Patient_annotation_table_extended[names(Proportion_Inflamatory_macrophages),"Outcome"] )
U = U[order(U$Outcome),]
U = U[complete.cases(U),]

pdf("Desktop/COVID_Italy/Figures/Outcome_macrophages.pdf",width = 3.5,height = 5.5,useDingbats = F)
par(las=1,bty="l")
prism.plots(Proportion~Outcome,U,cex=2,pch=21,bg=string.to.colors(U$Outcome,colors = c("green3","red3")),col="black",
            ylab="Proportion of macrophages in BAL")
dev.off()

TPM_cytokine= log2(1+t(data_count[c("SPP1","IL1B","IL6","IL1RN","CCL3"),])/colSums(data_count)*10^6)
Mean_cytokine_expression = aggregate(as.matrix(TPM_cytokine[r$clusters$PCA$Leiden%in%Cluster_total_cells,]),
                                     FUN = mean,by=list(r$clusters$PCA$Leiden[r$clusters$PCA$Leiden%in%Cluster_total_cells]))
rownames(Mean_cytokine_expression) = as.character(Mean_cytokine_expression$Group.1)
Mean_cytokine_expression = Mean_cytokine_expression[as.character(Cluster_total_cells),]
Mean_cytokine_expression = Mean_cytokine_expression[complete.cases(Mean_cytokine_expression),]
barplot(Mean_cytokine_expression$IL6,names.arg = Mean_cytokine_expression$Group.1)
Mean_cytokine_expression = Mean_cytokine_expression[order(Mean_cytokine_expression$SPP1,decreasing = F),]

pdf("Desktop/COVID_Italy/Figures/SPP1_cluster_expression.pdf",width = 4,height = 6)
par(las=1,bty="l")
barplot(Mean_cytokine_expression$SPP1,names.arg = Mean_cytokine_expression$Group.1,xlim=c(0,6),
        xlab="SPP1 expression (log2(1+TPM))",horiz = T)
dev.off()

#F) CA analysis of total blood (Figure 4G)

Cluster_distribution_Blood = table(Merged_annotation_count[Tissue_count=="Blood"],r$clusters$PCA$Leiden[Tissue_count=="Blood"])
Cluster_distribution_Blood = Cluster_distribution_Blood[,-Cluster_to_remove]
Cluster_distribution_Blood = Cluster_distribution_Blood[,!(colSums(Cluster_distribution_Blood)==0)]
Cluster_distribution_Blood_normalized = Cluster_distribution_Blood/rowSums(Cluster_distribution_Blood)

CA_total_Blood = FactoMineR::CA(Cluster_distribution_Blood[])

Type_samples = rownames(Cluster_distribution_Blood)
Type_samples[Type_samples=="Blood_8"]="23"
Type_samples[Type_samples=="Blood_4"]="24"
Type_samples[Type_samples=="Blood_25"]="24"
Type_samples[!Type_samples%in%c("23","24")] =21

pdf("Desktop/COVID_Italy/Figures/CA_global_virus.pdf",width = 6,height = 6)
par(las=1,bty="l")
plot(CA_total_Blood$row$coord[,1:2],pch=as.numeric(Type_samples),bg="firebrick3",cex=2)
dev.off()

Correlation_virus_status = apply(Cluster_distribution_Blood,MARGIN = 2,
                                 FUN = function(x) {cor(CA_total_Blood$row$coord[,1],x,method = "spearman")})

pdf("//home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Scree_plot_CA_blood_total.pdf",width = 8,height = 6)
par(las=1,bty="l")
barplot(CA_total_Blood$eig[1:10,2],ylim=c(0,40),names.arg = c(""),
        cex.lab=1.3,cex.axis = 1.3,ylab="Proportion of variance explained (%)",
        xlab="Ranked components")
dev.off()

Virus_associated_cluster = names(which(Correlation_virus_status==max(Correlation_virus_status,na.rm = T)))

#Generation of the heatmap of ISGs vs resting neutrophils (Figure 4H)
List_marker_genes = c()
for (k in c(Virus_associated_cluster)) {
  X = r$diffgenes$PCA$Leiden[[k]]
  X = X[X$highest,]
  X = X[X$M>1 & X$Z>10,]
  X = X[order(X$M,decreasing = T),]
  List_marker_genes = c(List_marker_genes,rownames(X[1:30,]))
}

List_marker_genes = List_marker_genes[!is.na(List_marker_genes)]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NA")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MALAT1")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MT-")]
List_marker_genes = unique(List_marker_genes)

TPM_total= log2(1+t(data_count[List_marker_genes,])/Matrix::colSums(data_count)*10^6)
Heatmap_data = c()

for (k in c(1,Cluster_virus)) {
  X = TPM_total[r$clusters$PCA$Leiden==k,List_marker_genes]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Heatmap_data = cbind(Heatmap_data,X)
}


Heatmap_data=as.matrix(Heatmap_data)
lim_values=quantile(as.numeric((Heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Heatmap_data[Heatmap_data<lim_values[1]]=lim_values[1]
Heatmap_data[Heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r$clusters$PCA$Leiden)
png("Desktop/COVID_Italy/Figures//Heatmap_neutro_virus.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(Heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_virus)])/sum(size_clusters[as.character(Cluster_virus)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Heatmap_data[,]))), 
     labels= ( rownames(Heatmap_data)[nrow(Heatmap_data):1]),
     las= 2, cex.axis=4.5, tick = F )
dev.off()



#IV)Lymphoid analysis (Figure 3)

#A)Pagoda2 analysis 

data_lymphoid = data_count[,Lymphoid_cells]
Zero_proportion = rowSums(data_lymphoid==0)/ncol(data_lymphoid)
Mean_expression = rowMeans(data_lymphoid)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2,subset = Mean_expression>0)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes_lymphoid = names(Zeros_excess_proportion[1:3000])


r_lymphoid <- Pagoda2$new(data_lymphoid,log.scale=F)
r_lymphoid$adjustVariance(plot=T,gam.k=10)
r_lymphoid$calculatePcaReduction(nPcs=70,odgenes = Selected_genes_lymphoid)
r_lymphoid$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r_lymphoid$getKnnClusters(method=multilevel.community,type='PCA')

Leiden_clustering = leiden(r_lymphoid$graphs$PCA)
names(Leiden_clustering) = colnames(data_lymphoid)
r_lymphoid$clusters$PCA$Leiden = factor(Leiden_clustering)

r_lymphoid$getDifferentialGenes(type = "PCA",upregulated.only = T,clusterType = "community",verbose = T,z.threshold = 3)
r_lymphoid$getDifferentialGenes(type = "PCA",upregulated.only = T,clusterType = "Leiden",verbose = T,z.threshold = 3)

umap_plot_lymphoid = umap(r_lymphoid$reductions$PCA,n_neighbors = 30,
                          n_components = 2,metric = "cosine",verbose = T)

#Figure 3A
pdf("Desktop/COVID_Italy/Figures/UMAP_lymphoid_cluster.pdf",width = 6.5,height = 6.5)
plot(umap_plot_lymphoid,pch=21,bg=string.to.colors(r_lymphoid$clusters$PCA$Leiden==1),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",cex=0.8)
plot(umap_plot_lymphoid,pch=21,bg=string.to.colors(Tissue_count[Lymphoid_cells],colors = c("firebrick3","dodgerblue")),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",cex=0.8)

dev.off()

pdf("Desktop/COVID_Italy/Figures/UMAP_lymphoid_gene_expression.pdf",width = 6,height = 6)
plot(umap_plot_lymphoid,pch=21,bg=color_convertion(r_lymphoid$counts[,"IFNG"]),main="IFNG",
     lwd=1,cex=1,xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_plot_lymphoid,pch=21,bg=color_convertion(r_lymphoid$counts[,"ZNF683"]),main="ZNF683",
     lwd=1,cex=1,xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()

#B)Heatmaps generation 

#1)Removing low quality clusters

Mean_expression_cluster_lymphoid = aggregate(as.matrix(r_lymphoid$counts[,Selected_genes_lymphoid]),FUN = mean,by=list(r_lymphoid$clusters$PCA$Leiden))
rownames(Mean_expression_cluster_lymphoid) = Mean_expression_cluster_lymphoid$Group.1
Mean_expression_cluster_lymphoid = t((Mean_expression_cluster_lymphoid)[,-1])
X = pheatmap(cor((Mean_expression_cluster_lymphoid),method = "spearman"),clustering_method = "ward.D2")
Order_cluster_lympho = X$tree_col$order

#Remove non lymphoid cluster : RBCs + basophils + cancer cells + platelet + myleoid/lymphoid doublets
Cluster_to_remove_lymphoid = c(9,18,17,13,20,8)

#2) Heatmap of the lymphoid compartment in the BAL (Figure 3D)

V = table(Tissue_count[Lymphoid_cells],r_lymphoid$clusters$PCA$Leiden)
V = V/rowSums(V)
V = round(100*V["BAL",])
Cluster_BAL_lympho = names(which(V>5))
Cluster_BAL_lympho = Cluster_BAL_lympho[!Cluster_BAL_lympho%in%Cluster_to_remove_lymphoid]
Cluster_BAL_lympho = Order_cluster_lympho[Order_cluster_lympho%in%Cluster_BAL_lympho]

List_marker_genes = c()
for (k in Cluster_BAL_lympho) {
  X = r_lymphoid$diffgenes$PCA$Leiden[[k]]
  X = X[X$highest,]
  X = X[X$M>1 & X$Z>10,]
  X = X[order(X$M,decreasing = T),]
  List_marker_genes = c(List_marker_genes,rownames(X[1:10,]))
}

List_marker_genes = List_marker_genes[!is.na(List_marker_genes)]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NA")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MALAT1")]

TPM_lymphoid = log2(1+t(data_lymphoid)/colSums(data_lymphoid)*10^6)
Heatmap_data = c()

for (k in Cluster_BAL_lympho) {
  X = TPM_lymphoid[r_lymphoid$clusters$PCA$Leiden==k,List_marker_genes]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Heatmap_data = cbind(Heatmap_data,X)
}


Heatmap_data=as.matrix(Heatmap_data)
lim_values=quantile(as.numeric((Heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Heatmap_data[Heatmap_data<lim_values[1]]=lim_values[1]
Heatmap_data[Heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r_lymphoid$clusters$PCA$Leiden)
png("Desktop/COVID_Italy/Figures//Heatmap_BAL_lymphoid_2.png", width=8000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(Heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_BAL_lympho)])/sum(size_clusters[as.character(Cluster_BAL_lympho)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Heatmap_data[,]))), 
     labels= ( rownames(Heatmap_data)[nrow(Heatmap_data):1]),
     las= 2, cex.axis=4.5, tick = F )
dev.off()

#3) Manually picked up genes based heatmap (Figure S3A)

List_marker_genes_reduced = c("IL7R","NKG7","TRBC1","TRDC","CD8A","CD4","GZMA","GZMB","IFNG","PDCD1","CTLA4","MKI67","CD74","CD79A","TCF4","IGKC")
Reduced_heatmap_data = c()

Order_cluster_lympho_bis = Order_cluster_lympho[!Order_cluster_lympho%in%Cluster_to_remove_lymphoid]
for (k in Order_cluster_lympho_bis) {
  X = TPM_lymphoid[r_lymphoid$clusters$PCA$Leiden==k,List_marker_genes_reduced]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Reduced_heatmap_data = cbind(Reduced_heatmap_data,X)
}

Reduced_heatmap_data=as.matrix(Reduced_heatmap_data)
lim_values=quantile(as.numeric((Reduced_heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Reduced_heatmap_data[Reduced_heatmap_data<lim_values[1]]=lim_values[1]
Reduced_heatmap_data[Reduced_heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r_lymphoid$clusters$PCA$Leiden)
png("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures//Heatmap_BAL_lymphoid_reduced.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(330,40,4,4))
image(rotate(Reduced_heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Order_cluster_lympho_bis)])/sum(size_clusters[as.character(Order_cluster_lympho_bis)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Reduced_heatmap_data[,]))), 
     labels= ( rownames(Reduced_heatmap_data)[nrow(Reduced_heatmap_data):1]),
     las= 2, cex.axis=4.5, tick = F )
dev.off()

#C) CA analysis (Figure 3B, C, E, F and S3B to E)

#1)Study of cell stype distribution

Distribution_lymphoid = table(Merged_annotation_count[Lymphoid_cells],r_lymphoid$clusters$PCA$Leiden)
Distribution_lymphoid = Distribution_lymphoid[,-Cluster_to_remove_lymphoid]
Distribution_lymphoid_normalized = Distribution_lymphoid/rowSums(Distribution_lymphoid)

pheatmap(Distribution_lymphoid_normalized,clustering_method = "ward.D2",
         annotation_row = Patient_annotation_table_extended[rownames(Distribution_lymphoid),2:4])

#2)CA analysis per se

CA_lymphoid = CA(Distribution_lymphoid,ncp = 5,graph = F)

pdf("//home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Scree_plot_CA_lymphoid.pdf",width = 8,height = 6)
par(las=1,bty="l")
barplot(CA_lymphoid$eig[1:10,2],ylim=c(0,40),names.arg = c(""),
        cex.lab=1.3,cex.axis = 1.3,ylab="Proportion of variance explained (%)",
        xlab="Ranked components")
dev.off()

plot(CA_lymphoid,axes = c(1,2))
pdf("Desktop/COVID_Italy/Figures/CA_lymphoid.pdf",width = 6,height = 6)
par(las=1,bty="l")
plot(CA_lymphoid$row$coord[,c(1,2)],pch=21,xlab="CA dimension 1",ylab="CA dimension 2",cex=1.5,
     bg=string.to.colors(grepl(rownames(Distribution_lymphoid),pattern = "BAL"),colors = c("dodgerblue3","firebrick3")))

dev.off()

plot(CA_lymphoid$row$coord[,1:2],pch=21,xlab="CA dimension 2",ylab="CA dimension 3",cex=1.5,
     bg=string.to.colors(Patient_annotation_table_extended[rownames(Distribution_lymphoid),"Type"]))
plot(CA_lymphoid$row$coord[,1:2],pch=21,xlab="CA dimension 2",ylab="CA dimension 3",cex=1.5,
     bg=color_convertion(Patient_annotation_table_extended[rownames(Distribution_lymphoid),"SOFA_score"]))


#C)Link between majors axes of CA and clinical features (Figure 3C)

#1) CA dimension 1 : association with SOFA score

CA_lymphoid_CA_1 = CA_lymphoid$row$coord[,1]
CA_lymphoid_CA_1 = CA_lymphoid_CA_1[grepl(names(CA_lymphoid_CA_1),pattern = "BAL")]

Correlation_CA_1_lymphoid_clinical = apply(Patient_annotation_table_extended[names(CA_lymphoid_CA_1),sapply(Patient_annotation_table_extended, is.numeric)],
                                           FUN = function(x) {cor(x,CA_lymphoid_CA_1,use = "complete")},MARGIN = 2)
Correlation_CA_1_lymphoid_clinical = Correlation_CA_1_lymphoid_clinical[order(Correlation_CA_1_lymphoid_clinical,decreasing = F)]

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA_lymphoid_1_correlation.pdf",width = 8,height = 6)
par(las=1,mar=c(6,8,4,4))
barplot(Correlation_CA_1_lymphoid_clinical,horiz = T,xlim=c(-1,1)
        ,xlab="Correlation with CA 1 (BAL)",cex.lab=1.3,col="grey70")
dev.off()
Correlation_CA_1_lymphoid_cytokine = apply(Cytokine_data_blood_data[names(CA_lymphoid_CA_1),sapply(Cytokine_data_blood_data, is.numeric)],
                                           FUN = function(x) {cor(sqrt(x),CA_lymphoid_CA_1,use = "complete")},MARGIN = 2)
Correlation_CA_1_lymphoid_cytokine = Correlation_CA_1_lymphoid_cytokine[order(Correlation_CA_1_lymphoid_cytokine)]
par(las=1,mar=c(6,14,4,4))
barplot(Correlation_CA_1_lymphoid_cytokine,horiz = T,xlim=c(-1,1)
        ,xlab="Correlation with CA 2 (BAL)",cex.lab=1.3,col="dodgerblue3")

m = lm(Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"]~CA_lymphoid_CA_1)

pdf("Desktop/COVID_Italy/Figures/SOFA_CA_lymphoid.pdf",width = 6,height = 6)
par(las=1,bty="l")
plot(CA_lymphoid_CA_1,Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"],pch=21,bg="dodgerblue3",
     ylab="SOFA score",xlab="CA dimension 1",cex=1.5,cex.lab=1.3,cex.axis=1.3)
abline(coef(m),lwd=2,lty=2,col="red")
dev.off()

pdf("Desktop/COVID_Italy/Figures/SOFA_regression.pdf",width = 6,height = 6)
par(las=1,bty="l")
plot(y=Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"],x=Distribution_lymphoid_normalized[names(CA_lymphoid_CA_1),"6"],
     ylab="SOFA score",xlab="Proportion of cytotoxic T cells in BAL",pch=21,bg="dodgerblue3",cex=1.5)
m = lm(Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"]~Distribution_lymphoid_normalized[names(CA_lymphoid_CA_1),"6"])
abline(coef(m),lwd=2,lty=2,col="red")

plot(y=Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"],x=Distribution_lymphoid_normalized[names(CA_lymphoid_CA_1),"1"],
     ylab="SOFA score",xlab="Proportion of Naive T cells in BAL",pch=21,bg="dodgerblue3",cex=1.5)
m = lm(Patient_annotation_table_extended[names(CA_lymphoid_CA_1),"SOFA_score"]~Distribution_lymphoid_normalized[names(CA_lymphoid_CA_1),"1"])
abline(coef(m),lwd=2,lty=2,col="red")
dev.off()

Correlation_CA_1_lymphoid_population = apply(Distribution_lymphoid_normalized[names(CA_lymphoid_CA_1),],
                                             FUN = function(x) {cor(x,CA_lymphoid_CA_1,method = "pearson")},MARGIN = 2)
Correlation_CA_1_lymphoid_population = Correlation_CA_1_lymphoid_population[as.character(Cluster_BAL_lympho)]

size_clusters_CA_1 = size_clusters[as.character(Cluster_BAL_lympho)]
pdf("Desktop/COVID_Italy/Figures/Color_bar_BAL_lymphocytes.pdf",width = 8,height = 4)
par(las=1)
barplot(matrix(size_clusters_CA_1),horiz = T,xaxt="n",col=color_convertion_correlation(Correlation_CA_1_lymphoid_population))
dev.off()

#B) CA dimension 2 : association with blood activation (Figure )

#1)CA analysis by itself
CA_lymphoid_CA_2 = CA_lymphoid$row$coord[,2]
CA_lymphoid_CA_2 = CA_lymphoid_CA_2[grepl(names(CA_lymphoid_CA_2),pattern = "Blood")]

Correlation_CA_2_lymphoid_clinical = apply(Patient_annotation_table_extended[names(CA_lymphoid_CA_2),sapply(Patient_annotation_table_extended, is.numeric)],
                                           FUN = function(x) {cor(x,CA_lymphoid_CA_2,use = "complete")},MARGIN = 2)
Correlation_CA_2_lymphoid_clinical = Correlation_CA_2_lymphoid_clinical[order(Correlation_CA_2_lymphoid_clinical)]
par(las=1,mar=c(6,8,4,4))
barplot(Correlation_CA_2_lymphoid_clinical,horiz = T,xlim=c(-1,1)
        ,xlab="Correlation with CA 2 (blood)",cex.lab=1.3,col="firebrick3")


Correlation_CA_2_lymphoid_cytokine = apply(Cytokine_data_blood_data[names(CA_lymphoid_CA_2),sapply(Cytokine_data_blood_data, is.numeric)],
                                           FUN = function(x) {cor(sqrt(x),CA_lymphoid_CA_2,use = "complete")},MARGIN = 2)
Correlation_CA_2_lymphoid_cytokine = Correlation_CA_2_lymphoid_cytokine[order(Correlation_CA_2_lymphoid_cytokine,decreasing = T)]
pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA_lymphoid_2_correlation.pdf",width = 8,height = 6)
par(las=1,mar=c(6,8,4,4))
barplot(Correlation_CA_2_lymphoid_clinical,horiz = T,xlim=c(-1,1)
        ,xlab="Correlation with CA 2 (Blood)",cex.lab=1.3,col="grey70")
dev.off()

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA2_lymphoid_correlation_features.pdf",width = 8,height = 8,useDingbats = F)
par(las=1,bty="l")
plot(Correlation_CA_2_lymphoid_cytokine,cex=1.5,pch=21,bg="firebrick3",xlab="Rank of the variable",
     ylim=c(-0.8,0.8),xaxs="i",yaxs="i",xlim=c(-1,80),ylab="Correlation with CA dimension 2",cex.lab=1.3,cex.axis=1.2)
dev.off()

pdf("Desktop/COVID_Italy/Figures/CA3_lymphoid_associated_features.pdf",width = 6,height = 6)
par(las=1,bty="l")
plot(CA_lymphoid_CA_2,log10(Cytokine_data_blood_data[names(CA_lymphoid_CA_2),"Neutrophils"]),pch=21,bg="firebrick3",
     ylab="Neutrophils (Log10)",xlab="CA dimension 3",cex=1.5)
plot(CA_lymphoid_CA_2,sqrt(Cytokine_data_blood_data[names(CA_lymphoid_CA_2),"IL6"]),pch=21,bg="firebrick3",
     ylab="IL1RA (Log10)",xlab="CA dimension 3",cex=1.5)
plot(CA_lymphoid_CA_2,log10(Cytokine_data_blood_data[names(CA_lymphoid_CA_2),"PDGF_BB"]),pch=21,bg="firebrick3",
     ylab="PDGF BB",xlab="CA dimension 3",cex=1.5)
abline(coef(m),lwd=2,lty=2,col="red")
dev.off()


U = data.frame(Score = CA_lymphoid_CA_2,
               Activated_NK = Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),"7"],
               Resting_NK = Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),"3"],
               Activated_gd = Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),"4"],
               Resting_gd = Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),"2"],
               
               Type = Patient_annotation_table_extended[names(CA_lymphoid_CA_2),"Type"])
U$Type[is.na(U$Type)] = "Severe"
U = U[order(U$Type),]

prism.plots(Activated_NK~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),cex=2,col="black",
            ylab="CA dimension 2",cex.lab=1.4,cex.axis=1.3)


pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA_dimension_2_blood.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
prism.plots(Score~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),cex=2,col="black",
            ylab="CA dimension 2",cex.lab=1.4,cex.axis=1.3)
plotSigBars(Score~Type,U)
dev.off()

#Proportion of activated/resting NK/gd-T cells (Figure 3E, 3F and S3E)
pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA2_lymphoid_associated_population.pdf",width = 9,height = 6,useDingbats = F)
par(las=1,bty="l",mfrow=c(1,2),mar=c(3,5,3,1))
prism.plots(Resting_NK~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=1.5,col="black",
            ylab="Proportion of resting NK",cex.lab=1.2,cex.axis=1.3)
plotSigBars(Resting_NK~Type,U)
prism.plots(Activated_NK~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=1.5,col="black",
            ylab="Proportion of activated NK",cex.lab=1.2,cex.axis=1.3)
plotSigBars(Activated_NK~Type,U)
prism.plots(Activated_gd~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=1.5,col="black",
            ylab="Proportion of resting gd T-cells",cex.lab=1.2,cex.axis=1.3)
plotSigBars(Activated_gd~Type,U)
prism.plots(Resting_gd~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=1.5,col="black",
            ylab="Proportion of activated gd T-cells",cex.lab=1.2,cex.axis=1.3)
plotSigBars(Resting_gd~Type,U)

dev.off()


#2)Heatmap of the blood lymphoid compartment (Figure 3G)

Correlation_CA_2_lymphoid_population = apply(Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),],
                                             FUN = function(x) {cor(x,CA_lymphoid_CA_2)},MARGIN = 2)
Correlation_CA_2_lymphoid_population = Correlation_CA_2_lymphoid_population[!is.na(Correlation_CA_2_lymphoid_population)]
Correlation_CA_2_lymphoid_population = Correlation_CA_2_lymphoid_population[abs(Correlation_CA_2_lymphoid_population)>0.5]

Cluster_Blood_lympho = names(Correlation_CA_2_lymphoid_population)
Cluster_Blood_lympho = Order_cluster_lympho[Order_cluster_lympho%in%Cluster_Blood_lympho]

List_marker_genes = c()
for (k in Cluster_Blood_lympho) {
  X = r_lymphoid$diffgenes$PCA$Leiden[[k]]
  X = X[X$highest,]
  X = X[X$M>1 & X$Z>10,]
  X = X[order(X$M,decreasing = T),]
  List_marker_genes = c(List_marker_genes,rownames(X[1:20,]))
}

List_marker_genes = List_marker_genes[!is.na(List_marker_genes)]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NA")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MALAT1")]

Heatmap_data = c()

for (k in Cluster_Blood_lympho) {
  X = TPM_lymphoid[r_lymphoid$clusters$PCA$Leiden==k,List_marker_genes]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Heatmap_data = cbind(Heatmap_data,X)
}


Heatmap_data=as.matrix(Heatmap_data)
lim_values=quantile(as.numeric((Heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Heatmap_data[Heatmap_data<lim_values[1]]=lim_values[1]
Heatmap_data[Heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r_lymphoid$clusters$PCA$Leiden)
png("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures//Heatmap_Blood_lymphoid.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(Heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Cluster_Blood_lympho)])/sum(size_clusters[as.character(Cluster_Blood_lympho)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Heatmap_data[,]))), 
     labels= ( rownames(Heatmap_data)[nrow(Heatmap_data):1]),
     las= 2, cex.axis=6, tick = F )
dev.off()

Correlation_CA_2_lymphoid_population = apply(Distribution_lymphoid_normalized[names(CA_lymphoid_CA_2),],
                                             FUN = function(x) {cor(x,CA_lymphoid_CA_2,method = "pearson")},MARGIN = 2)
Correlation_CA_2_lymphoid_population = Correlation_CA_2_lymphoid_population[as.character(Cluster_Blood_lympho)]

size_clusters_CA_2 = size_clusters[as.character(Cluster_Blood_lympho)]
pdf("Desktop/COVID_Italy/Figures/Color_bar_BAL_blood.pdf",width = 8,height = 4)
par(las=1)
barplot(matrix(size_clusters_CA_2),horiz = T,xaxt="n",col=color_convertion_correlation(Correlation_CA_2_lymphoid_population))
dev.off()


#C) Correlation between CA1 and CA2 ? (not shown in the paper)

CA_lymphoid_CA_1_renamed = CA_lymphoid_CA_1
names(CA_lymphoid_CA_1_renamed) = substr(names(CA_lymphoid_CA_1_renamed),start = 5,stop = 6)

CA_lymphoid_CA_2_renamed = CA_lymphoid_CA_2
names(CA_lymphoid_CA_2_renamed) = substr(names(CA_lymphoid_CA_2_renamed),start = 7,stop = 8)

par(las=1,bty="l")
plot(CA_lymphoid_CA_1_renamed[intersect(names(CA_lymphoid_CA_1_renamed),names(CA_lymphoid_CA_2_renamed))],
     CA_lymphoid_CA_2_renamed[intersect(names(CA_lymphoid_CA_1_renamed),names(CA_lymphoid_CA_2_renamed))],
     pch=21,bg="dodgerblue3",cex=1.5,xlab="CA 1",ylab="CA 2")
m = lm(CA_lymphoid_CA_2_renamed[intersect(names(CA_lymphoid_CA_1_renamed),names(CA_lymphoid_CA_2_renamed))]~CA_lymphoid_CA_1_renamed[intersect(names(CA_lymphoid_CA_1_renamed),names(CA_lymphoid_CA_2_renamed))])
abline(coef(m),lwd=2,lty=2,col="red")


#V)Analysing the blood neutrophil compartment (Figure 2 and S2)

#A)Single-cell blood neutrophil analysis 

Blood_neutrophil_cells = Neutrophil_cells & Tissue_count=="Blood"

data_neutrophil_Blood = data_count[,Blood_neutrophil_cells]

Zero_proportion = Matrix::rowSums(data_neutrophil_Blood==0)/ncol(data_neutrophil_Blood)
Mean_expression = Matrix::rowMeans(data_neutrophil_Blood)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2,subset = Mean_expression>0)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes_neutrophils_Blood = names(Zeros_excess_proportion[1:700])

r_neutrophil_Blood<- Pagoda2$new(data_neutrophil_Blood,log.scale=F)
Selected_genes_neutrophils_Blood = intersect(Selected_genes_neutrophils_Blood,colnames(r_neutrophil_Blood$counts))
r_neutrophil_Blood$adjustVariance(plot=T,gam.k=5)
r_neutrophil_Blood$calculatePcaReduction(nPcs=50,odgenes = Selected_genes_neutrophils_Blood)
r_neutrophil_Blood$makeKnnGraph(k=30,type='PCA',distance = "cosine",)
r_neutrophil_Blood$getKnnClusters(method=multilevel.community,type='PCA')

r_neutrophil_Blood$getDifferentialGenes(type = "PCA",upregulated.only = T,clusterType = "community",verbose = T,z.threshold = 3)
save(r_neutrophil_Blood,file = "/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Pagoda2_blood_neutro.RData")

#B)Cluster analysis + UMAP plot

Mean_expression_cluster_Blood_neutrophil = aggregate(as.matrix(r_neutrophil_Blood$counts[,Selected_genes_neutrophils_Blood]),
                                                     FUN = mean,by=list(r_neutrophil_Blood$clusters$PCA$community))
rownames(Mean_expression_cluster_Blood_neutrophil) = Mean_expression_cluster_Blood_neutrophil$Group.1
Mean_expression_cluster_Blood_neutrophil = t((Mean_expression_cluster_Blood_neutrophil)[,-1])
X = pheatmap(cor((Mean_expression_cluster_Blood_neutrophil),method = "spearman"),clustering_method = "ward.D2")
Order_cluster_neutro_Blood = X$tree_col$order

#Remove non lymphoid cluster : RBCs + basophils + cancer cells + platelet + myleoid/lymphoid doublets
Cluster_to_remove_neutrophil_Blood = c(14,11)
Distribution_neutrophil_Blood = table(Merged_annotation_count[Blood_neutrophil_cells],r_neutrophil_Blood$clusters$PCA$community)
Distribution_neutrophil_Blood = Distribution_neutrophil_Blood[,-Cluster_to_remove_neutrophil_Blood]
Distribution_neutrophil_Blood_normalized = Distribution_neutrophil_Blood/rowSums(Distribution_neutrophil_Blood)

N_cluster = length(unique(r_neutrophil_Blood$clusters$PCA$community))
N_cluster = N_cluster-2
optimal_palette =colorRamp(brewer.pal(11, "Spectral"))
optimal_palette = optimal_palette((1:N_cluster)/N_cluster)
optimal_palette = optimal_palette / 255
optimal_palette = rgb(optimal_palette)

#Generation of Figure 2A

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Umap_blood_neutrophils.pdf",width = 8,height = 8)
plot(umap_plot_neutrophils_Blood,pch=21,
     bg=string.to.colors(r_neutrophil_Blood$clusters$PCA$community[!r_neutrophil_Blood$clusters$PCA$community%in%Cluster_to_remove_neutrophil_Blood],colors = optimal_palette),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2",cex=0.7)
dev.off()

#C) CA analysis (Figure 2C and S2A/B)

pheatmap(Distribution_neutrophil_Blood_normalized,clustering_method = "ward.D2",
         annotation_row = Patient_annotation_table_extended[rownames(Distribution_neutrophil_Blood),2:4])
CA_neutrophil_Blood = CA(Distribution_neutrophil_Blood,ncp = 5,graph = F)
barplot(CA_neutrophil_Blood$eig[1:10,2],xlab="Ranked component",ylab="Proportion of variance explained (%)",ylim=c(0,30),cex.lab=1.4)
plot(CA_neutrophil_Blood,axes=c(2,3))

X = Patient_annotation_table_extended[rownames(Distribution_neutrophil_Blood_normalized),"Type"]
X[is.na(X)] = "Severe"

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA_blood_neutrophil.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
barplot(CA_neutrophil_Blood$eig[1:10,2],xlab="Ranked component",ylab="Proportion of variance explained (%)",ylim=c(0,30),cex.lab=1.4)
plot(CA_neutrophil_Blood$row$coord[,c(1,2)],cex=2,pch=21,bg=string.to.colors(X,c("brown3","orange","skyblue2")))
plot(CA_neutrophil_Blood$row$coord[,c(2,3)],cex=2,pch=21,bg=string.to.colors(X,c("brown3","orange","skyblue2")))
dev.off()

#D)Correlation of CA axis with clinical parameters (Figure S2D)

#CA 2 : spread the patient in a homogenous way...

CA_neutrophil_Blood_2 = CA_neutrophil_Blood$row$coord[,2]

Correlation_CA_neutrophil_Blood_2_clinical = apply(Patient_annotation_table_extended[names(CA_neutrophil_Blood_2),sapply(Patient_annotation_table_extended, is.numeric)],
                                                   FUN = function(x) {cor(x,CA_neutrophil_Blood_2,use = "complete")},MARGIN = 2)
Correlation_CA_neutrophil_Blood_2_cytokine= apply(Cytokine_data_blood_data[names(CA_neutrophil_Blood_2),sapply(Cytokine_data_blood_data, is.numeric)],
                                                  FUN = function(x) {cor(sqrt(x),CA_neutrophil_Blood_2,use = "complete")},MARGIN = 2)
Correlation_CA_neutrophil_Blood_2_cytokine = Correlation_CA_neutrophil_Blood_2_cytokine[order(Correlation_CA_neutrophil_Blood_2_cytokine,decreasing = T)]

Correlation_CA_neutrophil_Blood_2_population = apply(Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),],
                                                     FUN = function(x) {cor((x),CA_neutrophil_Blood_2,use = "complete",method = "spearman")},MARGIN = 2)

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/CA2_blood_neutrophil_correlation_features.pdf",width = 7,height = 7,useDingbats = F)
par(las=1,bty="l")
plot(Correlation_CA_neutrophil_Blood_2_cytokine,pch=21,bg="red",ylab="Correlation with CA2",xlab="Rank of the variable")
plot(Correlation_CA_neutrophil_Blood_2_cytokine,pch=21,bg="red",ylab="Correlation with CA2",xlab="Rank of the variable")
text(Correlation_CA_neutrophil_Blood_2_cytokine,labels = names(Correlation_CA_neutrophil_Blood_2_cytokine))
dev.off()

#E)Plotting proportion of specific neutrophil population (Figure 2D and S2C)

U = data.frame(Type = X,
               Outcome = Patient_annotation_table_extended[names(CA_neutrophil_Blood_2),"Outcome"],
               Proportion_LDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"13"],
               Proportion_resting_NDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"12"],
               Proportion_ISG_NDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"2"],
               Proportion_Protease_inhibitor_NDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"1"],
               Proportion_CD177_NDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"3"],
               Proportion_Infla_NDN = Distribution_neutrophil_Blood_normalized[names(CA_neutrophil_Blood_2),"10"])
U = U[order(U$Type),]

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Blood_neutrophil_prism.pdf",width = 4,height = 7,useDingbats = F)
par(las=1,bty="l")
prism.plots(Proportion_LDN~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=2,col="black",yaxs="i",ylab="Proportion of LDN",ylim=c(0,0.3))
kruskal.test(Proportion_LDN~Type,U)
prism.plots(Proportion_resting_NDN~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=2,col="black",ylim=c(0,1),yaxs="i",ylab="Proportion of resting NDN")
kruskal.test(Proportion_resting_NDN~Type,U)
prism.plots(Proportion_ISG_NDN~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=2,col="black",ylim=c(0,1),yaxs="i",ylab="Proportion of ISG+ NDN")
kruskal.test(Proportion_ISG_NDN~Type,U)
prism.plots(Proportion_Protease_inhibitor_NDN~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=2,col="black",ylim=c(0,0.2),yaxs="i",ylab="Proportion of Protease inhibitor NDN")
kruskal.test(Proportion_Protease_inhibitor_NDN~Type,U)
prism.plots(Proportion_CD177_NDN~Type,U,pch=21,bg=string.to.colors(U$Type,colors = c("skyblue2","orange","brown3")),
            cex=2,col="black",ylim=c(0,0.2),yaxs="i",ylab="Proportion of CD177+ NDN")
kruskal.test(Proportion_CD177_NDN~Type,U)

dev.off()


#E)Drawing of the heatmap (Figure 2A)

size_clusters = table(r_neutrophil_Blood$clusters$PCA$community)
size_clusters = size_clusters/sum(size_clusters)

Selected_clusters_neutrophils_Blood = Order_cluster_neutro_Blood[!Order_cluster_neutro_Blood%in%Cluster_to_remove_neutrophil_Blood]
Selected_clusters_neutrophils_Blood = Selected_clusters_neutrophils_Blood[Selected_clusters_neutrophils_Blood%in%names(which(size_clusters>0.01))]


List_marker_genes = c()
for (k in Selected_clusters_neutrophils_Blood) {
  X = r_neutrophil_Blood$diffgenes$PCA$community[[k]]
  X = X[X$highest,]
  X = X[X$M>1 & X$Z>10,]
  X = X[order(X$M,decreasing = T),]
  List_marker_genes = c(List_marker_genes,rownames(X[1:10,]))
}

List_marker_genes = List_marker_genes[!is.na(List_marker_genes)]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NA")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MALAT1")]
List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "NEAT1")]

List_marker_genes = List_marker_genes[!grepl(List_marker_genes,pattern = "MT-")]
List_marker_genes = unique(List_marker_genes)

TPM_total= log2(1+t(data_neutrophil_Blood[List_marker_genes,])/Matrix::colSums(data_neutrophil_Blood)*10^6)
Heatmap_data = c()

for (k in Selected_clusters_neutrophils_Blood) {
  X = TPM_total[r_neutrophil_Blood$clusters$PCA$community==k,List_marker_genes]
  X = t(X)
  X = X[,sample(1:ncol(X),size = ncol(X),replace = F)]
  Heatmap_data = cbind(Heatmap_data,X)
}


Heatmap_data=as.matrix(Heatmap_data)
lim_values=quantile(as.numeric((Heatmap_data)),c(0.0,0.995)) ##Trimming extreme values

Heatmap_data[Heatmap_data<lim_values[1]]=lim_values[1]
Heatmap_data[Heatmap_data>lim_values[2]]=lim_values[2]

size_clusters = table(r_neutrophil_Blood$clusters$PCA$community)
png("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures//Heatmap_Blood_neutrophils.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(Heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(Selected_clusters_neutrophils_Blood)])/sum(size_clusters[as.character(Selected_clusters_neutrophils_Blood)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((Heatmap_data[,]))), 
     labels= ( rownames(Heatmap_data)[nrow(Heatmap_data):1]),
     las= 2, cex.axis=4.5, tick = F )
dev.off()

##Creating the color bar for below the heatmap
List_cluster = r_neutrophil_Blood$clusters$PCA$community[!r_neutrophil_Blood$clusters$PCA$community%in%Cluster_to_remove_neutrophil_Blood]
X = string.to.colors(r_neutrophil_Blood$clusters$PCA$community[!r_neutrophil_Blood$clusters$PCA$community%in%Cluster_to_remove_neutrophil_Blood],optimal_palette)
names(X) = r_neutrophil_Blood$clusters$PCA$community[!r_neutrophil_Blood$clusters$PCA$community%in%Cluster_to_remove_neutrophil_Blood]
X = X[as.character(Selected_clusters_neutrophils_Blood)]

Bar_value = size_clusters[Selected_clusters_neutrophils_Blood]
Bar_value = Bar_value/sum(Bar_value)

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures//Barplot_cluster_blood_neutro.pdf",width = 8,height = 2)
barplot(matrix(Bar_value),horiz = T,col=X,xaxt="n")
dev.off()


#VI)Viral landscape analysis (Figure 4)

#A)Global VT analysis (Figure 4A)

VT_analysis = read.delim("Desktop/COVID_Italy/VT_analysis.txt",row.names = 1)
VT_analysis = VT_analysis[complete.cases(VT_analysis),]
VT_analysis = VT_analysis[order(VT_analysis$SARS.COV.2.reads),]

pdf("Desktop/COVID_Italy/Figures/Viral_load_BAL.pdf",width = 7,height =   7)
par(las=1,mar=c(5,6,4,4))
barplot(log10(1+VT_analysis$SARS.COV.2.reads),horiz = T,names.arg = rownames(VT_analysis),
        xlim=c(0,6),xlab="Number of SARS-CoV-2 reads detected (log10)",
        cex.lab=1.3,cex.axis = 1.3,col="dodgerblue3")
dev.off()

#B)Coverage analysis of SARS-CoV-2 in BAL8 samples (Figure 4B)

SARS_COV_2_patient_8= readGAlignments("/hdd/VT_analysis_Italy/VT_output/Merged_BAL8_R2/Viral_BAM_files/refseq|NC_045512|29903nt|Severe.bam",param = ScanBamParam(what =scanBamWhat()))
Coverage_SARS_COV_2_pos = coverage(SARS_COV_2_patient_8[strand(SARS_COV_2_patient_8)=="+",])
Coverage_SARS_COV_2_pos = Coverage_SARS_COV_2_pos$'refseq|NC_045512|29903nt|Severe'
Coverage_SARS_COV_2_neg = coverage(SARS_COV_2_patient_8[strand(SARS_COV_2_patient_8)=="-",])
Coverage_SARS_COV_2_neg = Coverage_SARS_COV_2_neg$"refseq|NC_045512|29903nt|Severe"

range_coverage = max(c(max(Coverage_SARS_COV_2_pos),max(Coverage_SARS_COV_2_neg)))
range_coverage = c(-range_coverage,range_coverage)*1.1


pdf("Desktop/COVID_Italy/Figures/Coverage_SARS_CoV_2.pdf",width = 10,height = 5.5)
par(bty="n",las=1)
plot(Coverage_SARS_COV_2_pos,ylim= range_coverage,yaxt="n",cex.axis=1.2,
     ylab="Coverage (10K reads)",xlab="Genome position",cex.lab=1.3,type="l")
points(-Coverage_SARS_COV_2_neg,type="l")
polygon(Coverage_SARS_COV_2_pos,col="firebrick3")
polygon(Coverage_SARS_COV_2_neg,col="firebrick3")
axis(side = 2,at = seq(from = round(range_coverage[1]/10000)*10000,to = round(range_coverage[2]/10000)*10000,length.out = 5),
     labels = abs(seq(from = round(range_coverage[1]/10000),to = round(range_coverage[2]/10000),length.out = 5)))
dev.off()

#C)Single-cell analysis of SARS-CoV-2 in BAL8 samples (Figure S4B and 4F)

Data_BAL_8_Corona = read.10x.matrices("//home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Data/Patient_8/BAL8_Corona/", n.cores = 1, verbose = T)
Corona_load_8 = colSums(Data_BAL_8_Corona[!grepl(rownames(Data_BAL_8_Corona),pattern = "ENSG"),])
Host_load_8 = colSums(Data_BAL_8_Corona[grepl(rownames(Data_BAL_8_Corona),pattern = "ENSG"),])

pdf("/home/data/Pierre/Untitled/home/pbost/Desktop/COVID_Italy/Figures/Cell_SARS_CoV_2_load.pdf",width = 7,height = 7)
par(las=1,bty="l")
plot(Host_load_8,Corona_load_8+1,log="xy",pch=21,bg="dodgerblue3",ylim=c(1,5000),xaxs="i",yaxs="i",
     xlab="Total host UMIs",ylab="Total SARS-CoV-2 UMIs",cex=1.4,cex.lab=1.3)
dev.off()

#Converting ENSEMBL gene ID to Gene symbol

x = read.delim("Correspondance_ENSEMBL.txt")
Correspondance_ENSEMBL = as.character(x$hgnc_symbol)
names(Correspondance_ENSEMBL) = x$ensembl_gene_id
Correspondance_ENSEMBL = Correspondance_ENSEMBL[Correspondance_ENSEMBL!=""]

rownames(Data_BAL_8_Corona) = Correspondance_ENSEMBL[rownames(Data_BAL_8_Corona)]
Data_BAL_8_Corona = Data_BAL_8_Corona[!is.na(rownames(Data_BAL_8_Corona)),]
Data_BAL_8_Corona = aggregate(as.matrix(Data_BAL_8_Corona),FUN = sum,by=list(rownames(Data_BAL_8_Corona)))
rownames(Data_BAL_8_Corona) = Data_BAL_8_Corona$Group.1
Data_BAL_8_Corona = as(as.matrix(Data_BAL_8_Corona[,-1]),"dgCMatrix")

Gene_size_8 = rowSums(Data_BAL_8_Corona)
hist(log10(Gene_size_8+1),100)
Data_BAL_8_Corona_count = Data_BAL_8_Corona[Gene_size_8>10,]

Zero_proportion = rowSums(Data_BAL_8_Corona_count==0)/ncol(Data_BAL_8_Corona_count)
Mean_expression = rowMeans(Data_BAL_8_Corona_count)

par(las=1)
plot(log10(Mean_expression),Zero_proportion,pch=21,bg="orange",xaxs='i',yaxs='i')
Zeros_excess_proportion = loess(Zero_proportion~log10(Mean_expression),degree = 2,subset = Mean_expression>0)
Zeros_excess_proportion = Zeros_excess_proportion$residuals
Zeros_excess_proportion = Zeros_excess_proportion[order(Zeros_excess_proportion,decreasing = T)]
Selected_genes_patient_8 = names(Zeros_excess_proportion[1:1000])


r_patient_8 <- Pagoda2$new(Data_BAL_8_Corona_count,log.scale=F)
r_patient_8$adjustVariance(plot=T,gam.k=10)
r_patient_8$calculatePcaReduction(nPcs=50,odgenes = Selected_genes_patient_8)
r_patient_8$makeKnnGraph(k=30,type='PCA',distance = "cosine")
r_patient_8$getKnnClusters(method=multilevel.community,type='PCA')

r_patient_8$getDifferentialGenes(type = "PCA",upregulated.only = T,clusterType = "community",verbose = T,z.threshold = 3)

umap_patient_8 = umap(r_patient_8$reductions$PCA,n_neighbors = 30,spread = 3,
                      n_components = 2,metric = "cosine",verbose = T)

plot(umap_patient_8,pch=21,bg=string.to.colors(r_patient_8$clusters$PCA$community==6),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(umap_patient_8,pch=21,bg=color_convertion(r_patient_8$counts[,"S100A8"]),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")


plot(umap_patient_8,pch=21,bg=color_convertion(log10(1+Corona_load_8)),main="Cluster distribution",
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
Viral_load = aggregate(Corona_load_8,by=list(r_patient_8$clusters$PCA$community),FUN = mean)
Viral_load = Viral_load$x
names(Viral_load) = 1:length(Viral_load)
Viral_load = Viral_load[order(Viral_load,decreasing = F)]
pdf("Desktop/COVID_Italy/Figures/Viral_load_Patient_8.pdf",width = 5,height = 5)
par(las=1)
barplot(Viral_load,horiz = T,xlim=c(0,30),xlab="Number of SARS-CoV-2 UMIs")
dev.off()

#D)Coverage analysis of HSV-1 in BAL4 sample (Figure 4D)

HSV1_patient_4= readGAlignments("/hdd/VT_analysis_Italy/VT_output/Merged_BAL4_R2/Viral_BAM_files/refseq|NC_001806|152222nt|Human.bam",param = ScanBamParam(what =scanBamWhat()))
Coverage_HSV_1_pos = coverage(HSV1_patient_4[strand(HSV1_patient_4)=="+",])
Coverage_HSV_1_pos = Coverage_HSV_1_pos$'refseq|NC_001806|152222nt'
Coverage_HSV_1_neg = coverage(HSV1_patient_4[strand(HSV1_patient_4)=="-",])
Coverage_HSV_1_neg = Coverage_HSV_1_neg$'refseq|NC_001806|152222nt'

range_coverage = max(c(max(Coverage_HSV_1_neg),max(Coverage_HSV_1_pos)))
range_coverage = c(-range_coverage,range_coverage)*1.1

pdf("Desktop/COVID_Italy/Figures/Coverage_HSV_1.pdf",width = 10,height = 5.5)
par(bty="n",las=1)
plot(Coverage_HSV_1_pos,ylim= range_coverage,yaxt="n",cex.axis=1.2,
     ylab="Coverage (1K reads)",xlab="Genome position",cex.lab=1.3,type="l")
points(-as.numeric(Coverage_HSV_1_neg),type="l")
polygon(Coverage_HSV_1_pos,col="firebrick3")
polygon(-as.numeric(Coverage_HSV_1_neg),col="firebrick3")
axis(side = 2,at = seq(from = round(range_coverage[1]/1000)*1000,to = round(range_coverage[2]/1000)*1000,length.out = 5),
     labels = abs(seq(from = round(range_coverage[1]/1000),to = round(range_coverage[2]/1000),length.out = 5)))
dev.off()


#E)Single-cell analysis of HSV-1 in BAL4 sample


Data_BAL_4_HSV = read.10x.matrices("Desktop/COVID_Italy/Data/Patient_4/BAL4_HSV_1//", n.cores = 1, verbose = T)
HVS_load_4 = colSums(Data_BAL_4_HSV[!grepl(rownames(Data_BAL_4_HSV),pattern = "ENSG"),])
Host_load_4 = colSums(Data_BAL_4_HSV[grepl(rownames(Data_BAL_4_HSV),pattern = "ENSG"),])

pdf("Desktop/COVID_Italy/Figures/Cell_HSV_1_load.pdf",width = 7,height = 7)
par(las=1,bty="l")
plot(Host_load_4,HVS_load_4+1,log="xy",pch=21,bg="dodgerblue3",ylim=c(1,1000),xaxs="i",yaxs="i",
     xlab="Total host UMIs",ylab="Total HSV-1 UMIs",cex=1.4,cex.lab=1.3)
dev.off()

Gene_size_4 = rowSums(Data_BAL_4_HSV[grepl(rownames(Data_BAL_4_HSV),pattern = "ENSG"),])
hist(log10(1+Gene_size_4),n=50)

y = rowSums(Data_BAL_4_HSV[!grepl(rownames(Data_BAL_4_HSV),pattern = "ENSG"),])


#F)Antibody titer across patients and association with viral-status (Figure 4E)

Ig_concentration = read.delim("Desktop/COVID_Italy/Immuno_globulin_data.txt",dec=",")
rownames(Ig_concentration) = rownames(Clinical_data)
Ig_concentration$Severity = Clinical_data$Severity

Type_samples = rownames(Ig_concentration[order(Ig_concentration$Severity),])
Type_samples[Type_samples=="Patient_8"]="23"
Type_samples[Type_samples=="Patient_4"]="24"
Type_samples[Type_samples=="Patient_25"]="24"
Type_samples[!Type_samples%in%c("23","24")] =21

pdf("Desktop/COVID_Italy/Figures/IgG_RBD.pdf",width = 5,height = 5,useDingbats = F)
par(las=1,bty="l")
prism.plots(IgG_RBD~Severity,Ig_concentration[order(Ig_concentration$Severity),],yaxs="i",col="black",cex=1.5,
            ylab = "IgG RBD (OD)",pch = as.numeric(Type_samples),
            bg=string.to.colors(Ig_concentration$Severity[order(Ig_concentration$Severity)],colors = c("skyblue2","orange","brown3")))
dev.off()         


#VI)Genetic analysis (Figure 3G and H)

#A)Transcriptomic based approach 

#1)Wich cellular compartment express the genetic hit ?

GWAS_hit = c("SLC6A20","LZTFL1","FYCO1","CXCR6","XCR1","CCR9")
GWAS_hit = intersect(GWAS_hit,rownames(data_count))
GWAS_hit_expression = log2(1+t(data_count[GWAS_hit,])/colSums(data_count[,])*10^6)
GWAS_hit_mean_expression = aggregate(as.matrix((GWAS_hit_expression)),
                                     FUN = mean,by =list(paste(Major_cell_type,Tissue_count)))
rownames(GWAS_hit_mean_expression) = GWAS_hit_mean_expression$Group.1
#Removing the 'blood epithelial' and 'unsassigned' types
GWAS_hit_mean_expression = GWAS_hit_mean_expression[c(-1,-9,-10),-1]

pdf("Desktop/COVID_Italy/Figures/GWAS_hit_expression.pdf",width = 6.5,height = 4.5)
par(las=1,mar=c(4,6,3,1))
for (k in GWAS_hit) {
  barplot(GWAS_hit_mean_expression[,k],ylim=c(0,max(GWAS_hit_mean_expression[,k]*1.2)),
          main=k,ylab="Mean gene expression \n (Log2(1+TPM))",cex.lab=1,
          density = c(50,rep(c(50,10),5)),
          col = c("green3",rep(c("orange","dodgerblue3","firebrick3"),each=2)))
  
}
dev.off()

#2)Among the BAL lymphocyte : which cells express CXCR6 ?

library(ggplot2)
U=data.frame(Expression=TPM_lymphoid[r_lymphoid$clusters$PCA$Leiden%in%Cluster_BAL_lympho,"CXCR6"],
             Condition=factor(r_lymphoid$clusters$PCA$Leiden[r_lymphoid$clusters$PCA$Leiden%in%Cluster_BAL_lympho],levels = Cluster_BAL_lympho))

pdf("Desktop/COVID_Italy/Figures/CXCR6_expression_lung_lymphocytes.pdf",width = 7,height = 5.5)
ggplot(U, aes(x=Condition, y=Expression))  +  geom_violin(trim=TRUE,scale = "width",fill='grey', color="black",bw=1.5)+ 
  scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(U[,"Expression"]))) +
  theme_classic() + ggtitle("CXCR6 expression") + theme(plot.title = element_text(size=15),axis.text=element_text(size = 15),axis.title = element_text(size = 12)) + scale_x_discrete(labels=1:length(Cluster_BAL_lympho),name=" ")  
dev.off()
GWAS_hit_mean_expression_BAL_lymphocyte = aggregate(GWAS_hit_expression[Lymphoid_cells,"CXCR6"], by = list(r_lymphoid$clusters$PCA$Leiden),FUN = mean)
GWAS_hit_mean_expression_BAL_lymphocyte = GWAS_hit_mean_expression_BAL_lymphocyte[Cluster_BAL_lympho,2]

#3)Among the BAL lymphocyte : which gene correlates with CXCR6 expression ?

CXCR6_gene_correlation = apply(r_lymphoid$counts[Tissue_count[Lymphoid_cells]=="BAL",Selected_genes_lymphoid],
                               MARGIN = 2,FUN = function(x) {cor(log2(1+x),log2(1+r_lymphoid$counts[Tissue_count[Lymphoid_cells]=="BAL","CXCR6"]))})
CXCR6_gene_correlation =  CXCR6_gene_correlation[order(CXCR6_gene_correlation,decreasing = T)]
