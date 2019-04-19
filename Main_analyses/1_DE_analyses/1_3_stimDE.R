####################################################################
### 1. Load dependencies & declare/load local functions ############
####################################################################

library(preprocessCore)
library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)

reproduce=TRUE
iterations=1000
name=paste0("1_3_stimDE")


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
### 2. Load input files ############################################
####################################################################

setwd(current)
cols = read.table(paste0("Inputs/1_DE_analyses/RNAseq_metadata.txt"))
reads = read.table("Inputs/1_DE_analyses/RNAseq_reads_matrix.txt")

####################################################################
### 3. Format reads matrices & metadata ############################
####################################################################

cols_LPS=cols[which(cols$PopDE_set==1 & cols$Condition %in% c("CTL","LPS")),]
cols_GARD=cols[which(cols$GARD_stim_set==1 & cols$Condition %in% c("CTL","GARD")),]

reads_LPS=reads[,which(colnames(reads) %in% rownames(cols_LPS))]
reads_GARD=reads[,which(colnames(reads) %in% rownames(cols_GARD))]

cols_LPS=pretty_up_cols(cols_LPS)
cols_GARD=pretty_up_cols(cols_GARD)

cols_LPS=cols_LPS[order(rownames(cols_LPS)),]
cols_GARD=cols_GARD[order(rownames(cols_GARD)),]
reads_GARD=reads_GARD[,order(colnames(reads_GARD))]
reads_LPS=reads_LPS[,order(colnames(reads_LPS))]

####################################################################
### 4. Run stimDE for each ligand (real data) ######################
####################################################################

stim_DE=function(reads,cols){
    
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    design=model.matrix(~Condition+Sex+fraction_assigned+Admixture+CD14+CD20+CD4,data=cols)
    v <- voom(dge,design,plot=FALSE)
    v_combat = ComBat(dat=as.matrix(v$E), batch=cols$Flowcell, mod=design, par.prior=TRUE)
    v$E=v_combat
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    return(list(fit,v))
}

DE_LPS=stim_DE(reads=reads_LPS,cols=cols_LPS)
fit_LPS=DE_LPS[[1]]
v_LPS=DE_LPS[[2]]

DE_GARD=stim_DE(reads=reads_GARD,cols=cols_GARD)
fit_GARD=DE_GARD[[1]]
v_GARD=DE_GARD[[2]]

if(reproduce)
{
    ### to reproduce results reported in Harrison et al. 2019, load these tables.
    shuffled_pvals_GARD=read.table(paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/stimDE_GARD.txt"))
    shuffled_pvals_LPS=read.table( paste0("Inputs/1_DE_analyses/permuted_pvalues/",name,"/stimDE_LPS.txt"))

}else{
    
    permute=function(cols){
        for(i in levels(cols$Individual))
        {
            set=which(cols$Individual %in% i)
            if(length(set)>1)
            {
                cols$Condition[set]=sample(cols$Condition[set])
            }else{
                if(runif(1)>0.5){
                    cols$Condition[set]=levels(cols$Condition)[1]}else{
                        cols$Condition[set]=levels(cols$Condition)[2]}
            }
        }
        return(cols)
    }
    
    cols_GARD_random=cols_GARD
    cols_LPS_random=cols_LPS
    
    
    for(iter in 1:iterations)
    {
        print(iter)
        cols_GARD_random=permute(cols_GARD_random)
        cols_LPS_random=permute(cols_LPS_random)
        
        DE_LPS_rand=stim_DE(reads=reads_LPS,cols=cols_LPS_random)
        DE_GARD_rand=stim_DE(reads=reads_GARD,cols=cols_GARD_random)
        
        permuted_LPS=topTable(DE_LPS_rand[[1]], coef="ConditionLPS", sort.by="none", adjust="BH",n=nrow(DE_LPS_rand[[1]]$coefficients))
        permuted_GARD=topTable(DE_GARD_rand[[1]], coef="ConditionGARD", sort.by="none", adjust="BH",n=nrow(DE_GARD_rand[[1]]$coefficients))
        
        if(iter==1)
        {
            shuffled_pvals_GARD <-data.frame(x=permuted_GARD$P.Value)
            shuffled_pvals_LPS <-data.frame(x=permuted_LPS$P.Value)
            
            rownames(shuffled_pvals_GARD)=rownames(permuted_GARD)
            rownames(shuffled_pvals_LPS)=rownames(permuted_LPS)
            
        } else {
            shuffled_pvals_GARD <- cbind(shuffled_pvals_GARD,x=permuted_GARD$P.Value)
            shuffled_pvals_LPS <- cbind(shuffled_pvals_LPS,x=permuted_LPS$P.Value)
        }
    }
    
    system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/permuted_p_values"))
    
    #Writing permuted p-value tables

    write.table(shuffled_pvals_GARD,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/stimDE_GARD.txt"))
    write.table(shuffled_pvals_LPS,paste0("Outputs/1_DE_analyses/",name,"/permuted_p_values/stimDE_LPS.txt"))


}

system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/results"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/LPS"))
system(paste0("mkdir -p Outputs/1_DE_analyses/",name,"/summary_stats/GARD"))

results_GARD=topTable(fit_GARD,coef="ConditionGARD",number=nrow(fit_GARD$coefficients))[,1:4]
results_LPS=topTable(fit_LPS,coef="ConditionLPS",number=nrow(fit_LPS$coefficients))[,1:4]


fdrs_GARD=permFDR(full_data=results_GARD,full_column_id="P.Value",perm_data=shuffled_pvals_GARD,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD"),significance_threshold=0.05)

fdrs_LPS=permFDR(full_data=results_LPS,full_column_id="P.Value",perm_data=shuffled_pvals_LPS,perm_column_ids="all",output_name=paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS"),significance_threshold=0.05)


###################################################################################
#### 7. Write results #############################################################
###################################################################################

## Format results tables.
results_GARD=fdrs_GARD$fdrs[,c(1,3,4,5)]
results_LPS=fdrs_LPS$fdrs[,c(1,3,4,5)]

colnames(results_GARD)=c("StimDE_log2FC","t","P_value","Fdr")
colnames(results_LPS)=c("StimDE_log2FC","t","P_value","Fdr")

## Write results tables
write.table(results_GARD,paste0("Outputs/1_DE_analyses/",name,"/results/stimDE_GARD.txt"))
write.table(results_LPS,paste0("Outputs/1_DE_analyses/",name,"/results/stimDE_LPS.txt"))

#Writing hits tables
write.table(fdrs_GARD$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/stimDE_hits.txt"))
write.table(fdrs_LPS$hits,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/stimDE_hits.txt"))

#Writing pi_o (estimated fraction of true null hypotheses)
write.table(fdrs_GARD$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/GARD/stimDE_pi_o.txt"))
write.table(fdrs_LPS$pi_o,paste0("Outputs/1_DE_analyses/",name,"/summary_stats/LPS/stimDE_pi_o.txt"))
