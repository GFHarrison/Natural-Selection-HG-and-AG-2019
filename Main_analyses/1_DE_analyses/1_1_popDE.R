####################################################################
### 1. Load dependencies & declare/load local functions ############
####################################################################

library(preprocessCore)
library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(reshape2)
library(cobs)

reproduce=TRUE
iterations=1000
name=paste0("1_1_popDE")

pretty_up_cols=function(cols){
    
    cols=refactorize(cols)
    
    cols=mean_center_clean(cols,"CD14")
    cols=mean_center_clean(cols,"CD4")
    cols=mean_center_clean(cols,"CD20")
    cols=mean_center_clean(cols,"CD3")
    cols=mean_center_clean(cols,"CD8")
    cols=mean_center_clean(cols,"CD56")
    cols=mean_center_clean(cols,"fraction_assigned")
    
    return(cols)
}
refactorize=function(cols_local){
    
    cols_local$Genotyping_ID=factor(cols_local$Genotyping_ID,levels=unique(as.character(cols_local$Genotyping_ID)))
    cols_local$Condition=factor(cols_local$Condition)
    cols_local$Flowcell=factor(cols_local$Flowcell)
    cols_local$Sex=factor(cols_local$Sex)
    
    return(cols_local)
}
mean_center_clean=function(cols_local,column){
    index=which(colnames(cols_local)==column)
    
    cols_local[,index]=cols_local[,index]-mean(cols_local[,index])
    return(cols_local)
}

current=getwd()
setwd("Github_codes/common_functions")
source("permFDR.R")

####################################################################
### 2. Load input files & initialize stochasticity seed ############
####################################################################

setwd(current)
cols = read.table(paste0("Inputs/1_DE_analyses/RNAseq_metadata.txt"))
reads=read.table("Inputs/1_DE_analyses/RNAseq_reads_matrix.txt")

####################################################################
### 3. Format reads matrices & metadata ############################
####################################################################

cols=cols[which(cols$PopDE_set==1),]
rownames(cols)=cols$Sample
cols=cols[order(rownames(cols)),]
reads=reads[,which(colnames(reads) %in% rownames(cols))]
reads=reads[,order(colnames(reads))]

length(which(rownames(cols)!=colnames(reads)))
# 0

### Get condition-specific reads matrices
reads_CTL=reads[,which(cols$Condition=="CTL")]
reads_GARD=reads[,which(cols$Condition=="GARD")]
reads_LPS=reads[,which(cols$Condition=="LPS")]

### Get condition-specific metadata tables
cols_CTL=cols[which(cols$Condition=="CTL"),]
cols_GARD=cols[which(cols$Condition=="GARD"),]
cols_LPS=cols[which(cols$Condition=="LPS"),]
### Drop absent levels & Mean-center technical covariates in condition-specific metadata
cols_CTL=pretty_up_cols(cols_CTL)
cols_LPS=pretty_up_cols(cols_LPS)
cols_GARD=pretty_up_cols(cols_GARD)

####################################################################
#### 4. Model DE within each condition (real dta) ##################
####################################################################

DE_WC=function(reads,cols){
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    design=model.matrix(~Sex+fraction_assigned+Admixture+CD14+CD20+CD4,data=cols)
    v <- voom(dge,design,plot=FALSE)
    v_combat = ComBat(dat=as.matrix(v$E), batch=cols$Flowcell, mod=design, par.prior=TRUE)
    v$E=v_combat
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    return(list(fit,v))
}

DE_CTL=DE_WC(reads=reads_CTL,cols=cols_CTL)
fit_CTL=DE_CTL[[1]]
v_CTL=DE_CTL[[2]]

DE_LPS=DE_WC(reads=reads_LPS,cols=cols_LPS)
fit_LPS=DE_LPS[[1]]
v_LPS=DE_LPS[[2]]

DE_GARD=DE_WC(reads=reads_GARD,cols=cols_GARD)
fit_GARD=DE_GARD[[1]]
v_GARD=DE_GARD[[2]]

###################################################################################
#### 5. Do permutation tests to define null model for PopDE #######################
###################################################################################

if(reproduce)
{
    ### to reproduce results reported in Harrison et al. 2019, load these tables.
    shuffled_pvals_adm_CTL=read.table( paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/popDE_CTL.txt"))
    shuffled_pvals_adm_GARD=read.table(paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/popDE_GARD.txt"))
    shuffled_pvals_adm_LPS=read.table( paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/popDE_LPS.txt"))
}else{
    
    cols_CTL_random=cols_CTL
    cols_GARD_random=cols_GARD
    cols_LPS_random=cols_LPS
    
    for(iter in 1:iterations)
    {
        if(iter%%100==0)print(iter)
        cols_GARD_random$Admixture=sample(cols_GARD_random$Admixture)
        cols_LPS_random$Admixture=sample(cols_LPS_random$Admixture)
        cols_CTL_random$Admixture=sample(cols_CTL_random$Admixture)
        
        DE_CTL_rand=DE_WC(reads=reads_CTL,cols=cols_CTL_random)
        DE_LPS_rand=DE_WC(reads=reads_LPS,cols=cols_LPS_random)
        DE_GARD_rand=DE_WC(reads=reads_GARD,cols=cols_GARD_random)
        
        permuted_adm_GARD=topTable(DE_GARD_rand[[1]], coef="Admixture", sort.by="none", adjust="BH",n=nrow(DE_GARD_rand[[1]]$coefficients))
        permuted_adm_CTL=topTable(DE_CTL_rand[[1]], coef="Admixture", sort.by="none", adjust="BH",n=nrow(DE_CTL_rand[[1]]$coefficients))
        permuted_adm_LPS=topTable(DE_LPS_rand[[1]], coef="Admixture", sort.by="none", adjust="BH",n=nrow(DE_LPS_rand[[1]]$coefficients))
        
        if(iter==1)
        {
            shuffled_pvals_adm_GARD <-data.frame(x=permuted_adm_GARD$P.Value)
            shuffled_pvals_adm_LPS <-data.frame(x=permuted_adm_LPS$P.Value)
            shuffled_pvals_adm_CTL <-data.frame(x=permuted_adm_CTL$P.Value)
            
            rownames(shuffled_pvals_adm_GARD)=rownames(permuted_adm_GARD)
            rownames(shuffled_pvals_adm_LPS)=rownames(permuted_adm_LPS)
            rownames(shuffled_pvals_adm_CTL)=rownames(permuted_adm_CTL)
            
        } else {
            shuffled_pvals_adm_GARD <- cbind(shuffled_pvals_adm_GARD,x=permuted_adm_GARD$P.Value)
            shuffled_pvals_adm_LPS <- cbind(shuffled_pvals_adm_LPS,x=permuted_adm_LPS$P.Value)
            shuffled_pvals_adm_CTL <- cbind(shuffled_pvals_adm_CTL,x=permuted_adm_CTL$P.Value)
            
        }
    }
    system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/permuted_p_values"))
    #Writing permutation tests p-values tables
    write.table(shuffled_pvals_adm_GARD,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/popDE_GARD.txt"))
    write.table(shuffled_pvals_adm_LPS,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/popDE_LPS.txt"))
    write.table(shuffled_pvals_adm_CTL,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/popDE_CTL.txt"))
}


###################################################################################
#### 6. Estimate FDRs using a permutation-tests-based null ########################
###################################################################################

system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/results"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/CTL"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/LPS"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/GARD"))

results_GARD=topTable(fit_GARD,coef="Admixture",number=nrow(fit_GARD$coefficients))[,1:4]
results_LPS=topTable(fit_LPS,coef="Admixture",number=nrow(fit_LPS$coefficients))[,1:4]
results_CTL=topTable(fit_CTL,coef="Admixture",number=nrow(fit_CTL$coefficients))[,1:4]

fdrs_GARD=permFDR(full_data=results_GARD,full_column_id="P.Value",perm_data=shuffled_pvals_adm_GARD,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD"),significance_threshold=0.05)

fdrs_LPS=permFDR(full_data=results_LPS,full_column_id="P.Value",perm_data=shuffled_pvals_adm_LPS,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS"),significance_threshold=0.05)

fdrs_CTL=permFDR(full_data=results_CTL,full_column_id="P.Value",perm_data=shuffled_pvals_adm_CTL,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/CTL"),significance_threshold=0.05)

###################################################################################
#### 7. Write results #############################################################
###################################################################################

## Format results tables.
results_GARD=fdrs_GARD$fdrs[,c(1,3,4,5)]
results_LPS=fdrs_LPS$fdrs[,c(1,3,4,5)]
results_CTL=fdrs_CTL$fdrs[,c(1,3,4,5)]

colnames(results_GARD)=c("PopDE_log2FC","t","P_value","Fdr")
colnames(results_LPS)=c("PopDE_log2FC","t","P_value","Fdr")
colnames(results_CTL)=c("PopDE_log2FC","t","P_value","Fdr")

## Write results tables
write.table(results_GARD,paste0("Outputs/1_DE_analyses/",name,"/results/popDE_GARD.txt"))
write.table(results_LPS,paste0("Outputs/1_DE_analyses/",name,"/results/popDE_LPS.txt"))
write.table(results_CTL,paste0("Outputs/1_DE_analyses/",name,"/results/popDE_CTL.txt"))

#Writing hits tables
write.table(fdrs_GARD$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/popDE_hits.txt"))
write.table(fdrs_LPS$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/popDE_hits.txt"))
write.table(fdrs_CTL$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/CTL/popDE_hits.txt"))

#Writing pi_o (estimated fraction of true null hypotheses)
write.table(fdrs_GARD$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/popDE_pi_o.txt"))
write.table(fdrs_LPS$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/popDE_pi_o.txt"))
write.table(fdrs_CTL$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/CTL/popDE_pi_o.txt"))


