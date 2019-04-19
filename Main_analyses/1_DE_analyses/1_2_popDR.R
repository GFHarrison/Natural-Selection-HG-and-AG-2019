####################################################################
### 1. Load dependencies & declare/load local functions ############
####################################################################

library(preprocessCore)
library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(qvalue)

reproduce=TRUE
iterations=1000
name=paste0("1_2_popDR")

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
#### 4. Regress out technical-covariate efects within condition ####
####################################################################

## First: estimate tech. covariate effects within condition.

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

## Second: regress them out.

correct_exp=function(v,fit){
    corrected_exp=v$E
    indexes=c(2,3,5,6,7)

    for(i in indexes)
    corrected_exp <- corrected_exp - fit$coefficients[,i]%*%t(fit$design[,i])
    return(corrected_exp)
}

corrected_exp_CTL=correct_exp(v_CTL,fit_CTL)
corrected_exp_LPS=correct_exp(v_LPS,fit_LPS)
corrected_exp_GARD=correct_exp(v_GARD,fit_GARD)

######################################################################################
#### 5. Build Fold-change matrices & weights & metadata tables from paired samples ###
######################################################################################

get_FC=function(exp_stim,exp_ref,weights_stim,weights_ref,cols_stim,cols_ref,cond_stim,cond_ref){
    
    FC=exp_stim
    weights_FC=FC
    
    cont=0
    
    cols_all=rbind(cols_ref,cols_stim)
    exp_all=cbind(exp_ref,exp_stim)
    weights_all=cbind(weights_ref,weights_stim)

    cols_all$Genotyping_ID=factor(cols_all$Genotyping_ID)
    cols_all$Condition=factor(cols_all$Condition)

    for(i in 1:length(levels(cols_all$Genotyping_ID)))
    {
        indexes=which(cols_all$Genotyping_ID==levels(cols_all$Genotyping_ID)[i])
        ref_index=indexes[which(cols_all$Condition[indexes]==cond_ref)]
        stim_index=indexes[which(cols_all$Condition[indexes]==cond_stim)]
        
        if(length(ref_index)==1 & length(stim_index)==1)
        {
            cont=cont+1
            
            FC[,cont]=exp_all[,stim_index]-exp_all[,ref_index]
            weights_FC[,cont]=1/(1/weights_all[,stim_index]+1/weights_all[,ref_index])
            
            colnames(FC)[cont]=levels(cols_all$Genotyping_ID)[i]
            colnames(weights_FC)[cont]=levels(cols_all$Genotyping_ID)[i]
        }
    }
    
    FC=FC[,1:cont]
    weights_FC=weights_FC[,1:cont]

    
    return(list(FC,weights_FC))
}

FC_LPS=get_FC(exp_stim=corrected_exp_LPS,exp_ref=corrected_exp_CTL,weights_stim=v_LPS$weights,weights_ref=v_CTL$weights,cols_stim=cols_LPS,cols_ref=cols_CTL,cond_stim="LPS",cond_ref="CTL")

weights_LPS=FC_LPS[[2]]
FC_LPS=FC_LPS[[1]]

FC_GARD=get_FC(exp_stim=corrected_exp_GARD,exp_ref=corrected_exp_CTL,weights_stim=v_GARD$weights,weights_ref=v_CTL$weights,cols_stim=cols_GARD,cols_ref=cols_CTL,cond_stim="GARD",cond_ref="CTL")

weights_GARD=FC_GARD[[2]]
FC_GARD=FC_GARD[[1]]

cols_GARD_FC=cols[which(cols$Condition=="GARD"),]
cols_LPS_FC=cols[which(cols$Condition=="LPS"),]

cols_LPS_FC=cols_LPS_FC[which(cols_LPS_FC$Genotyping_ID %in% colnames(FC_LPS)),]
cols_GARD_FC=cols_GARD_FC[which(cols_GARD_FC$Genotyping_ID %in% colnames(FC_GARD)),]

rownames(cols_LPS_FC)=cols_LPS_FC$Genotyping_ID
rownames(cols_GARD_FC)=cols_GARD_FC$Genotyping_ID

length(which(cols_LPS_FC$Genotyping_ID!=colnames(FC_LPS)))
length(which(cols_LPS_FC$Genotyping_ID!=colnames(weights_LPS)))

length(which(cols_GARD_FC$Genotyping_ID!=colnames(FC_GARD)))
length(which(cols_GARD_FC$Genotyping_ID!=colnames(weights_GARD)))


##################################################################
#### 6. Model PopDR effects of admixture on responses to ligands #
##################################################################

DR=function(cols,exp,weights,corrected)
{
    design=model.matrix(~Admixture,data=cols)
    fit <-lmFit(exp,weights=weights,design)
    fit <- eBayes(fit)
    return(fit)
}

DR_LPS=DR(cols_LPS_FC,FC_LPS,weights_LPS)
DR_GARD=DR(cols_GARD_FC,FC_GARD,weights_GARD)

###################################################################################
#### 7. Randomize admixture in permutation tests to define null model for PopDR ###
###################################################################################

if(reproduce)
{
    ### to reproduce results reported in Harrison et al. 2019, load these tables.
    shuffled_pvals_adm_GARD=read.table(paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/popDR_GARD.txt"))
    shuffled_pvals_adm_LPS=read.table( paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/popDR_LPS.txt"))
}else{
    
    cols_GARD_FC_random=cols_GARD_FC
    cols_LPS_FC_random=cols_LPS_FC
    for(iter in 1:iterations)
    {
        print(iter)
        cols_GARD_FC_random$Admixture=sample(cols_GARD_FC_random$Admixture)
        cols_LPS_FC_random$Admixture=sample(cols_LPS_FC_random$Admixture)
        
        GARD_rand=DR(cols=cols_GARD_FC_random,exp=FC_GARD,weights=weights_GARD)
        LPS_rand=DR(cols=cols_LPS_FC_random,exp=FC_LPS,weights=weights_LPS)
        
        permuted_adm_GARD=topTable(GARD_rand, coef="Admixture", sort.by="none", adjust="BH",n=nrow(GARD_rand$coefficients))
        permuted_adm_LPS=topTable(LPS_rand, coef="Admixture", sort.by="none", adjust="BH",n=nrow(LPS_rand$coefficients))
        
        if(iter==1)
        {
            shuffled_pvals_adm_GARD <-data.frame(x=permuted_adm_GARD$P.Value)
            shuffled_pvals_adm_LPS <-data.frame(x=permuted_adm_LPS$P.Value)
            
            rownames(shuffled_pvals_adm_GARD)=rownames(permuted_adm_GARD)
            rownames(shuffled_pvals_adm_LPS)=rownames(permuted_adm_LPS)
        } else {
            shuffled_pvals_adm_GARD <- cbind(shuffled_pvals_adm_GARD,x=permuted_adm_GARD$P.Value)
            shuffled_pvals_adm_LPS <- cbind(shuffled_pvals_adm_LPS,x=permuted_adm_LPS$P.Value)
            
        }
    }
    
    system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/permuted_p_values"))
    #Writing permutation tests p-values tables
    write.table(shuffled_pvals_adm_GARD,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/popDR_GARD.txt"))
    write.table(shuffled_pvals_adm_LPS,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/popDR_LPS.txt"))
}

###################################################################################
#### 8. Estimate FDRs using a permutation-tests-based null ########################
###################################################################################

system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/results"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/LPS"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/GARD"))

results_GARD=topTable(DR_GARD,coef="Admixture",number=nrow(DR_GARD$coefficients))[,1:4]
results_LPS=topTable(DR_LPS,coef="Admixture",number=nrow(DR_LPS$coefficients))[,1:4]

fdrs_GARD=permFDR(full_data=results_GARD,full_column_id="P.Value",perm_data=shuffled_pvals_adm_GARD,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD"),significance_threshold=0.1)

fdrs_LPS=permFDR(full_data=results_LPS,full_column_id="P.Value",perm_data=shuffled_pvals_adm_LPS,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS"),significance_threshold=0.1)

###################################################################################
#### 9. Write results #############################################################
###################################################################################

## Format results tables.
results_GARD=fdrs_GARD$fdrs[,c(1,3,4,5)]
results_LPS=fdrs_LPS$fdrs[,c(1,3,4,5)]
colnames(results_GARD)=c("PopDR_effect_log2FC","t","P_value","Fdr")
colnames(results_LPS)=c("PopDR_effect_log2FC","t","P_value","Fdr")

## Write results tables
write.table(results_GARD,paste0("Outputs/1_DE_analyses/",name,"/results/popDR_GARD.txt"))
write.table(results_LPS,paste0("Outputs/1_DE_analyses/",name,"/results/popDR_LPS.txt"))

#Writing hits tables
write.table(fdrs_GARD$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/popDR_hits.txt"))
write.table(fdrs_LPS$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/popDR_hits.txt"))

#Writing pi_o (estimated fraction of true null hypotheses)
write.table(fdrs_GARD$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/popDR_pi_o.txt"))
write.table(fdrs_LPS$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/popDR_pi_o.txt"))
