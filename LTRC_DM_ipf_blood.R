library(gdsfmt)
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)


library(tidyverse)
library(data.table)
library(sjPlot)
library(janitor)
library(tableone)
library(hablar)
library(survival)
library("RColorBrewer")
library(VennDiagram)
library(ggsignif)
library(survminer)
library(mitch)
# library(EnsDb.Hsapiens.v79)
library(readr)

library(edgeR)
library(limma)
library(sva)
library(gplots)
library(BiocManager)
library(fgsea)
library(DESeq2)
library(GSVA)
# bigsheet <- read_csv("../../../proj/regeps/regep00/studies/COPDGene/analyses/rebdh/ltrc/cleaning/ltrcSamplePhenoIdMap_20200818.csv")
meta <- read_csv("/proj/regeps/regep00/studies/LTRC/data/phenotype/data/freezes/20210310/ltrcSamplePhenoIdMap_20210310.csv")

bigsheet<-read.table('/proj/regeps/regep00/studies/LTRC/data/phenotype/data/freezes/20210310/ltrcLongTopMedHarm_20210310.csv',sep=',',header=TRUE)
# 
# 
# meta3<-meta %>% mutate(newid = str_extract(newid, "[^/]+$"))
# 
big<-bigsheet[c('patid','race','gender','age')]
BIG<-drop_na(big)
# 
meta4<-meta[c('patid','topmedId','saphId.cdnm','Project')]
meta4$patid<-as.numeric(meta4$patid)
meta5<-meta4[meta4['Project']=='Methylation.blood',]


jj<-merge(meta5,BIG,by='patid')
jj<-unique(jj)

write.table(jj,"~/analyses/LTRC/pheno_meta.txt", sep="\t")

pheno_meta <- read.csv("~/analyses/LTRC/pheno_meta.txt", sep="\t")


lung<-readRDS('~/../resiq/meffil/QC/LTRC.blood/minfi/results/ltrc.blood.beta.rcp')
# 
# lung<-readRDS('~/../resiq/meffil/QC/LTRC.lung/minfi/results/ltrc.lung.beta.rcp')

# analysis.final.modCopd <- fread(paste0('/proj/edith/regeps/regep00/studies/COPDGene/analyses/reagh/LTRC_ILD/modCopd_pathConservRna.pheno_022221.csv'),colClasses = list(character=1))


analysis.final.modCopd <- fread(paste0('/proj/edith/regeps/regep00/studies/COPDGene/analyses/reagh/LTRC_ILD/ipfRna.pheno_022221.csv'),colClasses = list(character=1))


RNApheno.modCopd <- analysis.final.modCopd %>% mutate(modCopd.control=as.factor(ipf.clinpath),
                                                      gender=as.factor(gender),
                                                      race=as.factor(race),
                                                      smoking_current=as.factor((smoking_current)))

table(RNApheno.modCopd$modCopd.control, exclude =F)

RNApheno.modCopd$patid<-as.numeric(RNApheno.modCopd$patid)
# pheno_meta$patid<-pheno_meta$ALIAS

RNApheno.modCopd<-merge(pheno_meta,RNApheno.modCopd,on=patid)
date<-date()
date<-str_replace_all(string=date, pattern=" ", repl="")
## Read in and prepare count data

## Clean and prep counts file 
### Limit samples using rna bam ids from phenotype file
rnaBamId.modCopd <- as.character(RNApheno.modCopd$topmedId)
ww<-str_split(colnames(lung),'_')#[[1:length(colnames(lung))]][1]
qq<-data.frame(ww)
qq<-data.frame(t(qq))
colnames(lung)<-qq$X1
col.num <- which(colnames(lung) %in% RNApheno.modCopd$topmedId)
counts.clean <- lung[,col.num]



### Make DGEList object 
# counts.modCopd <- DGEList(counts.clean)

# dim(counts.modCopd)

col.num <- which(RNApheno.modCopd$topmedId %in% colnames(counts.clean))
meta <- RNApheno.modCopd[col.num,]
meta$bmi<-NULL
meta$ht_cm<-NULL
meta$wt_kg<-NULL
meta<-unique(meta)
## Create design matrices

jj<-colnames(design.modCopd)
jj<-jj[2:6]

AA<-t(combn(jj,5))
AA<-rbind(AA,cbind(t(combn(jj,4)),0))
AA<-rbind(AA,cbind(t(combn(jj,3)),0,0))
AA<-rbind(AA,cbind(t(combn(jj,2)),0,0,0))
AA<-rbind(AA,cbind(t(combn(jj,1)),0,0,0,0))
AA<-data.frame(AA)
AA$Z<-AA %>% unite("Z", X1:X2:X3:X4:X5, na.rm = TRUE, remove = TRUE,sep='+')

for(i in 1:nrow(AA)) {
    row <- AA[i,]
    assign("model",row$Z)
    # do stuff with row

design.modCopd <- model.matrix(as.formula(paste("~", paste(model),sep='')), data=meta)

# nullmod <- model.matrix(~ipf.clinpath+age+race+gender+smoking_packyears, data=meta)


### linear modeling 
fit.modCopd <- lmFit(counts.clean, design.modCopd)#sv)


fit.modCopd <- eBayes(fit.modCopd)

## Examine results
results.modCopd<-topTable(fit.modCopd, coef=2, number=Inf) ## make sure this isnt (intercept) ## hardcode this
results.modCopd <- results.modCopd %>% add_rownames(var = 'ensemble') ## threshold of meaningful difference
results.modCopd$ensemble <- gsub("\\..*","",results.modCopd$ensemble) ## 
results0.01FDR.modCopd <- subset(results.modCopd, results.modCopd$adj.P.Val < 0.1,)

write.csv(results.modCopd,paste("~/analyses/LTRC/blood_methyl_IPF_",date,"_",paste(row$Z),".csv",sep=''), row.names = F)
# write.csv(results0.01FDR.modCopd,"~/analyses/LTRC/blood_methyl_results0.1FDR.IPF_032921_svX.csv", row.names = F)
}
