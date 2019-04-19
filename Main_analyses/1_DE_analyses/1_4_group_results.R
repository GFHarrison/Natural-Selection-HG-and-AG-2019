{
    DE_CTL=read.table("Outputs/1_DE_analyses/1_1_popDE/results/popDE_CTL.txt")
    DE_LPS=read.table("Outputs/1_DE_analyses/1_1_popDE/results/popDE_LPS.txt")
    DE_GARD=read.table("Outputs/1_DE_analyses/1_1_popDE/results/popDE_GARD.txt")
    
    DR_LPS=read.table("Outputs/1_DE_analyses/1_2_popDR/results/popDR_LPS.txt")
    DR_GARD=read.table("Outputs/1_DE_analyses/1_2_popDR/results/popDR_GARD.txt")
    
    stim_LPS=read.table("Outputs/1_DE_analyses/1_3_stimDE/results/stimDE_LPS.txt")
    stim_GARD=read.table("Outputs/1_DE_analyses/1_3_stimDE/results/stimDE_GARD.txt")
    
    ## Subset and clean format
    
    ##
    stim_LPS=stim_LPS[,c("StimDE_log2FC","P_value","Fdr")]
    colnames(stim_LPS)=c("LPS_stimulation_logFC","LPS_stimulation_p_value","LPS_stimulation_FDR")
    
    stim_GARD=stim_GARD[,c("StimDE_log2FC","P_value","Fdr")]
    colnames(stim_GARD)=c("GARD_stimulation_logFC","GARD_stimulation_p_value","GARD_stimulation_FDR")
    
    DE_CTL=DE_CTL[,c("PopDE_log2FC","P_value","Fdr")]
    colnames(DE_CTL)=c("PopDE_CTL_logFC","PopDE_CTL_p_value","PopDE_CTL_FDR")
    
    DE_LPS=DE_LPS[,c("PopDE_log2FC","P_value","Fdr")]
    colnames(DE_LPS)=c("PopDE_LPS_logFC","PopDE_LPS_p_value","PopDE_LPS_FDR")
    
    DE_GARD=DE_GARD[,c("PopDE_log2FC","P_value","Fdr")]
    colnames(DE_GARD)=c("PopDE_GARD_logFC","PopDE_GARD_p_value","PopDE_GARD_FDR")
    
    DR_LPS=DR_LPS[,c("PopDR_effect_log2FC","P_value","Fdr")]
    colnames(DR_LPS)=c("PopDR_LPS_logFC","PopDR_LPS_p_value","PopDR_LPS_FDR")
    
    DR_GARD=DR_GARD[,c("PopDR_effect_log2FC","P_value","Fdr")]
    colnames(DR_GARD)=c("PopDR_GARD_logFC","PopDR_GARD_p_value","PopDR_GARD_FDR")
    
    
    ##let's be sure that the genes are properly ordered
    
    length(which(rownames(stim_LPS)!=rownames(stim_GARD)))
    length(which(rownames(stim_LPS)!=rownames(DE_CTL)))
    length(which(rownames(stim_LPS)!=rownames(DE_LPS)))
    length(which(rownames(stim_LPS)!=rownames(DE_GARD)))
    length(which(rownames(stim_LPS)!=rownames(DR_LPS)))
    length(which(rownames(stim_LPS)!=rownames(DR_GARD)))
    
    DR_LPS=DR_LPS[order(rownames(DR_LPS)),]
    DR_GARD=DR_GARD[order(rownames(DR_GARD)),]
    stim_LPS=stim_LPS[order(rownames(stim_LPS)),]
    stim_GARD=stim_GARD[order(rownames(stim_GARD)),]
    DE_CTL=DE_CTL[order(rownames(DE_CTL)),]
    DE_LPS=DE_LPS[order(rownames(DE_LPS)),]
    DE_GARD=DE_GARD[order(rownames(DE_GARD)),]
    
    length(which(rownames(stim_LPS)!=rownames(stim_GARD)))
    length(which(rownames(stim_LPS)!=rownames(DE_CTL)))
    length(which(rownames(stim_LPS)!=rownames(DE_LPS)))
    length(which(rownames(stim_LPS)!=rownames(DE_GARD)))
    length(which(rownames(stim_LPS)!=rownames(DR_LPS)))
    length(which(rownames(stim_LPS)!=rownames(DR_GARD)))
    
    
    results=cbind(DE_CTL,DE_GARD,DE_LPS,DR_GARD,DR_LPS,stim_GARD,stim_LPS)
    
}

write.csv(results,"Outputs/1_DE_analyses/table_S2.csv")

