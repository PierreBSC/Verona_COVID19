library(fifer)
library(pheatmap)


#I)Prism-plot for immuno-suppression

Immuno_supression_data = read.delim("Desktop/Immune_suppression_COVID19/Immune_suppression.txt",dec = ",")


#A)Across cell type

Rearranged_data = data.frame(Immuno_suppression = c(Immuno_supression_data$SUPPRESSION.BY.CD14.,
                                                    Immuno_supression_data$SUPPRESSION.BY.SN.CD14.,
                                                    Immuno_supression_data$SUPPRESSION.BY.SN.CD66b..LDN,
                                                    Immuno_supression_data$SUPPRESSION.BY.SN.CD66b..NDN),
               Cell_type = c(rep("CD14",nrow(Immuno_supression_data)),rep("CD14_SN",nrow(Immuno_supression_data)),rep("CD66_LDN_SN",nrow(Immuno_supression_data)),rep("CD66_NDN_SN",nrow(Immuno_supression_data))))


pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/IS_by_cell_type.pdf",width = 7.5,height = 7.5,useDingbats = F)
par(las=1)
prism.plots(Immuno_suppression~Cell_type,Rearranged_data,bty="l",
            pch=21,bg=string.to.colors(Rearranged_data$Cell_type,c("skyblue2","orange","brown3","grey")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression (%)",cex.lab=1.3)
dev.off()

m = aov(Immuno_suppression~Cell_type,Rearranged_data)
anova(m)

#B)Across conditions/survival

pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/IS_by_type.pdf",width = 5,height = 7.5,useDingbats = F)
par(las=1)
Immuno_supression_data = Immuno_supression_data[order(Immuno_supression_data$Type),]
prism.plots(SUPPRESSION.BY.CD14.~Type,Immuno_supression_data,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data$Type,c("skyblue2","orange","brown3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression (%)",cex.lab=1.3)
plotSigBars(SUPPRESSION.BY.CD14.~Type,Immuno_supression_data)
prism.plots(SUPPRESSION.BY.SN.CD14.~Type,Immuno_supression_data,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data$Type,c("skyblue2","orange","brown3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by CD14 SN(%)",cex.lab=1.3)
prism.plots(SUPPRESSION.BY.SN.CD66b..NDN~Type,Immuno_supression_data,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data$Type,c("skyblue2","orange","brown3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by NDN SN(%)",cex.lab=1.3)
prism.plots(SUPPRESSION.BY.SN.CD66b..LDN~Type,Immuno_supression_data,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data$Type,c("skyblue2","orange","brown3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by NDN SN(%)",cex.lab=1.3)

dev.off()

Immuno_supression_data_filtered = Immuno_supression_data[Immuno_supression_data$Type=="Severe",]
Immuno_supression_data_filtered = Immuno_supression_data_filtered[order(Immuno_supression_data_filtered$Outcome),]

pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/IS_by_outcome.pdf",width = 4,height = 7.5,useDingbats = F)
par(las=1)
prism.plots(SUPPRESSION.BY.CD14.~Outcome,Immuno_supression_data_filtered,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data_filtered$Outcome,c("green3","red3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by monocytes (%)",cex.lab=1.3)
plotSigBars(SUPPRESSION.BY.CD14.~Outcome,Immuno_supression_data_filtered)
prism.plots(SUPPRESSION.BY.SN.CD14.~Outcome,Immuno_supression_data_filtered,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data_filtered$Outcome,c("green3","red3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by monocytes supernatant (%)",cex.lab=1.3)
plotSigBars(SUPPRESSION.BY.SN.CD14.~Outcome,Immuno_supression_data_filtered)
prism.plots(SUPPRESSION.BY.SN.CD66b..LDN~Outcome,Immuno_supression_data_filtered,bty="l",
            pch=21,bg=string.to.colors(Immuno_supression_data_filtered$Outcome,c("green3","red3")),col="black",cex.axis=1.3,
            ylim=c(0,100),xaxs="i",yaxs="i",cex=2,ylab="T-cell suppression by LDN supernatant (%)",cex.lab=1.3)
plotSigBars(SUPPRESSION.BY.SN.CD66b..LDN~Outcome,Immuno_supression_data_filtered)

dev.off()


#II)Study of the possible mechanisms 


#A)FACS staining

FACS_data = read.delim("Desktop/Immune_suppression_COVID19/FACS_CD14_monocytes.txt",dec = ",")
Outcome_point = as.character(FACS_data$Outcome.ICU)
Outcome_point = ifelse(Outcome_point=="alive",yes = 21,no = 24)

#1)ARG1 : best predictor -> Fit by a Hill function

par(las=1,bty="l")
plot(log(FACS_data$MFI.ARG1.monocytes.CD14.+1),FACS_data$SUPPRESSION.BY.CD14.,pch=Outcome_point,
     bg=string.to.colors(FACS_data$Severity.of.Disease,c("brown3","orange","skyblue")),
     xlim=c(4,10),ylim=c(0,100),xaxs="i",yaxs="i",cex=2,
     xlab="ARG1 MFI (log)",ylab="Suppression by monocytes (%)")

x = log(FACS_data$MFI.ARG1.monocytes.CD14.+1)
y = FACS_data$SUPPRESSION.BY.CD14.
kept_sample = !is.na(x) & !is.na(y)
x = x[kept_sample]
y = y[kept_sample]
cor(x,y,method = "spearman")

nonlin_mod = nls(y~E_min + E/(1+(K/x)^n),start=list(E_min=10,E=100,n=30,K=7)) 
x = seq(4,10,length.out = 100)
y = coef(nonlin_mod)[1] + coef(nonlin_mod)[2] / (1+(coef(nonlin_mod)[4]/x)^coef(nonlin_mod)[3])


pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/ARG1_suppression.pdf",width = 8,height = 6,useDingbats = F)
par(las=1,bty="l")
plot(log(FACS_data$MFI.ARG1.monocytes.CD14.+1),FACS_data$SUPPRESSION.BY.CD14.,pch=Outcome_point,
     bg=string.to.colors(FACS_data$Severity.of.Disease,c("brown3","orange","skyblue")),
     xlim=c(4,10),ylim=c(0,100),xaxs="i",yaxs="i",cex=2,cex.lab=1.3,cex.axis=1.3,
     xlab="Monocytes ARG1 MFI (log)",ylab="Suppression by monocytes (%)")
lines(x,y,lwd=2,lty=2)
dev.off()


#2)HLA-DR : no clear relation but nice clustering

pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/HLA_DR_suppression.pdf",width = 8,height = 6,useDingbats = F)
par(las=1,bty="l")
plot(log(FACS_data$MFI.HLA_DR.on.CD14..monocytes+1),FACS_data$SUPPRESSION.BY.CD14.,pch=Outcome_point,
     bg=string.to.colors(FACS_data$Severity.of.Disease,c("brown3","orange","skyblue")),
     xlim=c(1,6),ylim=c(0,100),xaxs="i",yaxs="i",cex=2,
     xlab="HLA-DR MFI (log)",ylab="Suppression by monocytes (%)")
dev.off()


#3)PDL1 : weak association


pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/PDL1_suppression.pdf",width = 8,height = 6,useDingbats = F)
par(las=1,bty="l")
plot(log(FACS_data$MFI.PD.L1.monocytes.CD14.),FACS_data$SUPPRESSION.BY.CD14.,pch=Outcome_point,
     bg=string.to.colors(FACS_data$Severity.of.Disease,c("brown3","orange","skyblue")),
     ylim=c(0,100),xaxs="i",yaxs="i",cex=2,xlim=c(4,10),
     xlab="PD-L1 MFI (log)",ylab="Suppression by monocytes (%)")
dev.off()

x = log(FACS_data$MFI.PD.L1.monocytes.CD14.+1)
y = FACS_data$SUPPRESSION.BY.CD14.
kept_sample = !is.na(x) & !is.na(y)
x = x[kept_sample]
y = y[kept_sample]
cor(x,y,method = "spearman")


#B)Cytokine production

Cytokine_monocyte_production = read.delim("Desktop/Immune_suppression_COVID19/Cytokine_monocytes_CD14.txt",dec = ",")
Cytokine_monocyte_production[,-c(1:4)] = log(1+Cytokine_monocyte_production[,-c(1:4)] )
Cytokine_monocyte_production = Cytokine_monocyte_production[complete.cases(Cytokine_monocyte_production),]

pheatmap(cor(Cytokine_monocyte_production[,-c(1:4)]),clustering_method = "ward.D2")

Outcome_point = as.character(Cytokine_monocyte_production$Outcome.ICU)
Outcome_point = ifelse(Outcome_point=="alive",yes = 21,no = 24)

Correlation_suppression = apply(Cytokine_monocyte_production[,-c(1:4)],MARGIN = 2,
                                FUN = function(x) {cor(x,Cytokine_monocyte_production$SUPPRESSION.BY.SN.CD14.,method = "spearman")})
Correlation_suppression = Correlation_suppression[order(Correlation_suppression,decreasing = T)]

pdf("Desktop/Clinical_COVID_paper/Figure_2_v2/Sub_panel/Cytokine_correlation.pdf",width = 8,height = 6,useDingbats = F)
par(las=1,bty="l")
plot(Correlation_suppression,pch=21,bg="red",ylim=c(-1,1),yaxs="i",cex=2,
     ylab="Spearman's correlation with Immuno-suppresssion",xlab="Rank of the variable")
dev.off()
