####################################################################
### 1. Load dependencies ###########################################
####################################################################

library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(relaimpo)
options(width=10000,max.print=1000000000)

get_pvalues_and_qvalues=function(tab){
    
    tab$Delta_relimp_p=empPvals(tab$Delta_relimp_true,tab$Delta_relimp_null)
    tab$Delta_relimp_fdr=qvalue(tab$Delta_relimp_p)$qvalue
    
    return(tab)
}

system(paste0("mkdir -p Outputs/4_Delta_PVE/"))

####################################################################
### 2. Get DeltaPVE stats within each condition ####################
####################################################################

for(cond in c("LPS","GARD","CTL"))
{
    
    #######################################################################
    ### 2.1. Load & Prepare input files format ############################
    #######################################################################
    
    
    metadata = read.table(paste0("Inputs/1_DE_analyses/RNAseq_metadata.txt"))
    reads=read.table("Inputs/1_DE_analyses/RNAseq_reads_matrix.txt")
    CTL_true=read.table("Inputs/4_Delta_PVE/genotype_files/true/CTL.txt")
    LPS_true=read.table("Inputs/4_Delta_PVE/genotype_files/true/LPS.txt")
    GARD_true=read.table("Inputs/4_Delta_PVE/genotype_files/true/GARD.txt")
    CTL_perm=read.table("Inputs/4_Delta_PVE/genotype_files/permuted/CTL.txt")
    LPS_perm=read.table("Inputs/4_Delta_PVE/genotype_files/permuted/LPS.txt")
    GARD_perm=read.table("Inputs/4_Delta_PVE/genotype_files/permuted/GARD.txt")
    
    genotypes_ok=get(paste0(cond,"_true"))
    genotypes_null=get(paste0(cond,"_perm"))
    
    ### Metadata & reads treatments pre-analyses:
    
    metadata=metadata[which(metadata$EQTL_set==1 & metadata$Condition==cond & metadata$Genotyping_ID %in% colnames(genotypes_ok)),]
    metadata$Condition=factor(metadata$Condition)
    metadata$Flowcell=factor(metadata$Flowcell)
    reads=reads[,which(colnames(reads) %in% rownames(metadata))]
    reads=reads[,order(metadata$Genotyping_ID)]
    reads=reads[which(rownames(reads) %in% genotypes_ok$Gene),]
    reads=reads[order(rownames(reads)),]
    
    metadata=metadata[order(metadata$Genotyping_ID),]
    metadata$Genotyping_ID=factor(metadata$Genotyping_ID,levels=unique(as.character(metadata$Genotyping_ID)))
    rownames(metadata)=metadata$Genotyping_ID
    colnames(reads)=metadata$Genotyping_ID
    
    model=read.table(paste0("Inputs/4_Delta_PVE/DE_stats/popDE_",cond,".txt"))
    
    model=model[which(rownames(model)%in% rownames(reads)),]
    model=model[order(rownames(model)),]
    
    #genotypes_ok=read.table(paste0("clean_input_tables/",cond,"/genotypes_true.txt"))
    #genotypes_null=read.table(paste0("clean_input_tables/",cond,"/genotypes_permuted.txt"))
    #reads=read.table(paste0("clean_input_tables/",cond,"/reads_ok.txt"))
    #metadata=read.table(paste0("clean_input_tables/",cond,"/metadata_ok.txt"))
    
    model$relimp=0
    model$relimp_genotype_true=0
    model$relimp_genotype_null=0
    
    model$Delta_relimp_true=0
    model$Delta_relimp_null=0
    
    rownames(genotypes_ok)=genotypes_ok$Gene
    rownames(genotypes_null)=genotypes_null$Gene
    rownames(metadata)=metadata$Genotyping_ID
    
    genotypes_ok=genotypes_ok[,which(colnames(genotypes_ok) %in% rownames(metadata))]
    genotypes_null=genotypes_null[,which(colnames(genotypes_null) %in% rownames(metadata))]
    
    metadata$CD4=metadata$CD4-mean(metadata$CD4)
    metadata$CD14=metadata$CD14-mean(metadata$CD14)
    metadata$CD20=metadata$CD20-mean(metadata$CD20)
    metadata$fraction_assigned=metadata$fraction_assigned-mean(metadata$fraction_assigned)
    
    length(which(rownames(model)!=rownames(reads)))
    length(which(rownames(model)!=rownames(genotypes_ok)))
    length(which(rownames(model)!=rownames(genotypes_null)))
    
    length(which(colnames(reads)!=colnames(genotypes_ok)))
    length(which(colnames(reads)!=colnames(genotypes_null)))
    length(which(colnames(reads)!=colnames(metadata$Genotyping_ID)))
    
    ################################################################
    ### 2.2. Get PVE by admixture in basic design ##################
    ################################################################

    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    design = model.matrix(~Sex+fraction_assigned+CD14+CD4+CD20+Admixture,data=metadata)
    v=voom(dge,design,plot=FALSE)
    v_combat = ComBat(dat=as.matrix(v$E), batch=metadata$Flowcell, mod=design, par.prior=TRUE)
    v$E=v_combat
    
    ############################################################################
    ### 2.3 Recompute PVE after adding for each gene top true of random EQTLs ##
    ############################################################################
    
    hits_DE=which(model$Fdr<0.2)


    ## ~ 3-4 s per gene.
    for(iter in hits_DE)
    {
        print(iter)
        
        ###Build an extended metadata matrix containing the genotypes of the top snps (full and perm)
        
        metadata_iter=cbind(metadata,SNP_true=t(genotypes_ok[iter,]),SNP_null=t(genotypes_null[iter,]))
        colnames(metadata_iter)[(length(metadata_iter)-1):length(metadata_iter)]=c("SNP_true","SNP_null")
        
        ### Remove samples having missing entries for some of the genotypes
        set_missing_data_true=which(metadata_iter$SNP_true==-9)
        set_missing_data_null=which(metadata_iter$SNP_null==-9)
        
        if(length(set_missing_data_true)>0){
            metadata_iter_true=metadata_iter[-set_missing_data_true,]}else{metadata_iter_true=metadata_iter}
        
        if(length(set_missing_data_null)>0){
            metadata_iter_null=metadata_iter[-set_missing_data_null,]}else{metadata_iter_null=metadata_iter}
        
        set_full=which(colnames(reads) %in% rownames(metadata_iter_true))
        set_null=which(colnames(reads) %in% rownames(metadata_iter_null))
        
        reads_true=reads[,set_full]
        reads_null=reads[,set_null]
        
        metadata_iter_true$Flowcell=factor(metadata_iter_true$Flowcell)
        metadata_iter_null$Flowcell=factor(metadata_iter_null$Flowcell)
        
        metadata_iter_true$Sex=factor(metadata_iter_true$Sex)
        metadata_iter_null$Sex=factor(metadata_iter_null$Sex)
        
        dge_true <- DGEList(counts=reads_true)
        dge_true <- calcNormFactors(dge_true)
        design_true = model.matrix(~Sex+fraction_assigned+CD14+CD4+CD20+Admixture+SNP_true,data=metadata_iter_true)
        v_true=voom(dge_true,design_true,plot=FALSE)
        v_combat_true = ComBat(dat=as.matrix(v_true$E), batch=metadata_iter_true$Flowcell, mod=design_true, par.prior=TRUE)
        v_true$E=v_combat_true
        
        dge_null <- DGEList(counts=reads_null)
        dge_null <- calcNormFactors(dge_null)
        design_null = model.matrix(~Sex+fraction_assigned+CD14+CD4+CD20+Admixture+SNP_null,data=metadata_iter_null)
        v_null=voom(dge_null,design_null,plot=FALSE)
        v_combat_null = ComBat(dat=as.matrix(v_null$E), batch=metadata_iter_null$Flowcell, mod=design_null, par.prior=TRUE)
        v_null$E=v_combat_null
        
        gene_model=lm(v$E[iter,]~Sex+fraction_assigned+CD14+CD4+CD20+Admixture,data=metadata,weights=v$weights[iter,])
        gene_model_true=lm(v_true$E[iter,]~Sex+fraction_assigned+CD14+CD4+CD20+Admixture+SNP_true,data=metadata_iter_true,weights=v_true$weights[iter,])
        gene_model_null=lm(v_null$E[iter,]~Sex+fraction_assigned+CD14+CD4+CD20+Admixture+SNP_null,data=metadata_iter_null,weights=v_null$weights[iter,])
        
        #rel_impo=suppressWarnings(calc.relimp(gene_model,type=c("lmg","genizi","car")))
        rel_impo=suppressWarnings(calc.relimp(gene_model,type="last"))
        rel_impo_true=suppressWarnings(calc.relimp(gene_model_true,type="last"))
        rel_impo_null=suppressWarnings(calc.relimp(gene_model_null,type="last"))
        
        model$relimp[iter]=rel_impo$last["Admixture"]
        model$relimp_genotype_true[iter]=rel_impo_true$last["Admixture"]
        model$relimp_genotype_null[iter]=rel_impo_null$last["Admixture"]
        
    }
    model=model[hits_DE,]
    genotypes_ok=get(paste0(cond,"_true"))
    genotypes_ok=genotypes_ok[which(genotypes_ok$Gene %in% rownames(model)),]
    
    model$Delta_relimp_true=(model$relimp-model$relimp_genotype_true)/model$relimp
    model$Delta_relimp_null=(model$relimp-model$relimp_genotype_null)/model$relimp
    model=get_pvalues_and_qvalues(model)
    model$SNP=genotypes_ok$SNP
    model=model[,c(1:4,12,5:11)]
    colnames(model)=c("log2FC_popDE","t_popDE","P_value_popDE","Fdr_popDE","top_EQTL","PVE_base","PVE_after_EQTL","PVE_after_top_random_EQTL","Delta_PVE","Delta_PVE_random","P_value_Delta_PVE","Fdr_Delta_PVE")
    
    model=model[,c(1,3,4,9,11,12)]
    
    write.csv(model,paste0("Outputs/4_Delta_PVE/Tab_S6_Delta_PVE_",cond,".csv"))
    #write.table(model,paste0("Outputs/4_Delta_PVE/txt/Tab_S6_Delta_PVE_",cond,".txt"))

}

### Hits:
### CTL: 163, 269.

ctl=read.csv(paste0("Outputs/4_Delta_PVE/Tab_S6_Delta_PVE_CTL.csv"))
lps=read.csv(paste0("Outputs/4_Delta_PVE/Tab_S6_Delta_PVE_LPS.csv"))
gard=read.csv(paste0("Outputs/4_Delta_PVE/Tab_S6_Delta_PVE_GARD.csv"))

th=0.05

length(which(ctl$Fdr_Delta_PVE<th))
length(which(lps$Fdr_Delta_PVE<th))
length(which(gard$Fdr_Delta_PVE<th))

## 163, 46, 61

th=0.1

length(which(ctl$Fdr_Delta_PVE<th))
length(which(lps$Fdr_Delta_PVE<th))
length(which(gard$Fdr_Delta_PVE<th))

## 269, 117, 196

