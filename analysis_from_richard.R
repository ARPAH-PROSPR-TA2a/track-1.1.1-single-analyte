# lmvm_crVSal_covsexagebmirin
# Initial: February 21, 2023
# Latest:  February 21, 2023
# Richard Friedman
# based on
# based upon age_diet_lmvm_stage1_sex.R
# To analyze dan Belsky's mRNA data
# Gives 
# 1. CR 12 months compared to AL 12 months
#    relative to CR 0 months compared to AL 0 months.
# 2. CR 24 months compared to AL 24 months
#    relative to CR 0 months compared to AL 0 months.
#    using a single model for all 135 sammples, 
#    and taking into account the following
#    covariates:
#    1. Sex (b.female or a.male) (categorical fixed effect)
#    2. Agebl (continous fixed effect)
#    3. B_mbmi (continuous fixed effect)
#    4. Rin (contiuous fixed effect)
#    The duplicate correlation method is used, using subject
#    as a blocking effect.
######

# Input:
# data/comparisons.txt list of comparisons
# data/targets.3.wcov.txt target file with 2 outliers removed and above covariates included.
# data/age_diet_mrna_counts.3.txt countfile with 2 outliers removed.

library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(xlsx)


comparisons<-read.table("data/comparisons.txt",header=T,as.is=T)
no_genes_fdrle05<-rep(0,nrow(comparisons))
no_genes_ple001<-rep(0,nrow(comparisons))
fdr_at_ple001<-rep(1,nrow(comparisons))
no_genes<-cbind(no_genes_fdrle05,no_genes_ple001,fdr_at_ple001)
colnames(no_genes)<-c("no_genes_fdrle05","no_genes_ple001","fdr_at_ple001")
rownames(no_genes)<-comparisons$comparison


targets<-readTargets("data/targets.3.wcov.txt")


TreatTime<-factor(targets$TreatTime)
Sex<-factor(targets$Sex,levels=c("a.male","b.female"))

# Input counts 

dir.create("results/")
dir.create("results/global/")
counts<-read.table("data/age_diet_mrna_counts.3.txt",header=T,sep="\t",row.names=1,as.is=T)


y<-DGEList(counts=counts, genes=rownames(counts))
isexpr<-rowSums(y$counts>=10) >= 3
hasannot<-rowSums(is.na(y$genes))==0
y<-y[isexpr & hasannot,]
y$samples$lib.size<- colSums(y$counts)
y$samples$lib.size
y<-calcNormFactors(y)
des<-model.matrix(~0+TreatTime+ Sex+Agebl+Bmbmi+RIN,targets)

write.xlsx(des,"results/global/design.matrix.xlsx")

png("results/global/voomplot.png")
v<-voomWithQualityWeights(y,design=des,normalize.method="none",plot=TRUE)
dev.off()

cor <- duplicateCorrelation(v, des, block=targets$Subject)
sink("results/global/correlation.txt")
cat("within subject correlation =", cor$consensus.correlation)
sink()

fit <- lmFit(object=v, design=des, 
             block=targets$Subject, correlation=cor$consensus.correlation)

contrast.matrix <- makeContrasts(cr12moVSal12moVScr0moVSal0mo="TreatTimecr12mo-TreatTimeal12mo-TreatTimecr0mo+TreatTimeal0mo",
                                 cr24moVSal24moVScr0moVSal0mo="TreatTimecr24mo-TreatTimeal24mo-TreatTimecr0mo+TreatTimeal0mo",
                                 levels=des)


write.xlsx(contrast.matrix,"results/global/contrast.matrix.xlsx")

fit<-contrasts.fit(fit,contrast.matrix)
fit<-eBayes(fit)


####################################All comparisons #################################
for(i in 1:nrow(comparisons)){
  
  
  comp<-topTable(fit,coef=i,number=Inf,adjust.method="BH",sort.by="P")
  ids <- comp$genes
  
  annot<- select(org.Hs.eg.db, keys=ids, columns=c("SYMBOL","ENTREZID","GENENAME"), keytype="ENTREZID")
  
  comp.annot<-merge(comp,annot,by.x="genes",by.y="ENTREZID",all.x=T)
  
  
  comp.annot<-cbind(comp.annot[,1],comp.annot[,8:9],comp.annot[,2:7])
  colnames(comp.annot)[1]<-"ENTREZID"
  colnames(comp.annot)[8]<-"fdr"
  colnames(comp.annot)[4]<-"log2FC"
  
  comp.annot<-comp.annot[order(comp.annot$P.Value),]
  
  comp.annot.pein<-subset(comp.annot,select=c(SYMBOL,log2FC,fdr))
  comp.annot.pein.p <-subset(comp.annot,select=c(SYMBOL,log2FC,P.Value))
  
  dir<-paste0("results/",comparisons[i,1])
  dir.create(dir)
  
  #
  write.csv(comp.annot,paste0(dir,"/",comparisons[i,1],".csv"),row.names=F)
  png(paste0(dir,"/",comparisons[i,1],".his.png"))
  hist(comp$P.Value, breaks=seq(0,1,0.05),col="black",border="white")
  dev.off()
  
  comp.annot.pein<-comp.annot.pein[!duplicated(comp.annot.pein$SYMBOL),]
  comp.annot.pein<-comp.annot.pein[!is.na(comp.annot.pein$SYMBOL),]
  write.table(comp.annot.pein,paste0(dir,"/",comparisons[i,1],".pein.txt"),row.names=F,sep="\t",quote=F)
  
  comp.annot.pein.p<-comp.annot.pein.p[!duplicated(comp.annot.pein.p$SYMBOL),]
  comp.annot.pein.p<-comp.annot.pein.p[!is.na(comp.annot.pein.p$SYMBOL),]
  write.table(comp.annot.pein.p,paste0(dir,"/",comparisons[i,1],".pein.p.txt"),row.names=F,sep="\t",quote=F)
  
  comp.annot.fdrl05<-comp.annot[comp.annot$fdr<=0.05,]
  if(nrow(comp.annot.fdrl05)>0)
  {
    write.csv(comp.annot.fdrl05,paste0(dir,"/",comparisons[i,1],".fdrl05.csv"),row.names=F)
    no_genes[i,1]<-nrow(comp.annot.fdrl05)
  }
  
  comp.annot.pl001<-comp.annot[comp.annot$P.Value <=0.001,]
  if(nrow(comp.annot.pl001)>0)
  {
    
    write.csv(comp.annot.pl001,paste0(dir,"/",comparisons[i,1],"pl001.csv"),row.names=FALSE)
    no_genes[i,2]<-length(webstaltinp.pl001.sym)
    no_genes[i,3]<-comp.annot[nrow(comp.annot.pl001),8]
    write.table(webstaltinp.pl001.sym,paste0(dir,"/",comparisons[i,1],".wbgsltinp.pl001.sym.txt"),row.names=F,col.names=F,sep="\t",quote=F)
  }
  
  
  de<-subset(comp.annot,select=c(log2FC,P.Value))
  de<-na.omit(de)
  
  de$col<-8
  de$col[de$P.Value<0.05 & de$log2FC >0]<- "red"
  de$col[de$P.Value<0.05 &  de$log2FC  < 0]<- "blue"
  
  
  de$P.Value <- -log(de$P.Value,10)
  
  pdf(paste0(dir,"/",comparisons[i,1],".volcano.pdf"))
  par (mar=c(5,5,3,1))
  plot(de$log2FC,de$P.Value,col=de$col,ylim=c(0,6),xlim=c(-3,3),xlab="log2 (Fold Change)",ylab="-log10(Pvalue)",pch=16,cex=1.0)
  dev.off()
  
  
  rm(ids)
  rm(annot)
  rm(comp)
  rm(comp.annot)
  rm(webstaltinp.pl001.sym)
  rm(comp.annot.pl001)
  rm(comp.annot.fdrl05)
  rm(webstaltinp.fdrl05.abslfcg.6.sym)
}

write.xlsx(no_genes,"results/global/number.genes.xlsx")

sink("results/global/sessioninfo.txt")
sessionInfo()
sink()

rm(list=ls())
