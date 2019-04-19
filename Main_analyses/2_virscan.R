####################################################################################
### 1. Load and groom fdata, metadata & pvalues table for epitope-individual tests
####################################################################################

library(limma)
library(ggrepel)
library(ggplot2)


pvals=read.table("Inputs/2_virscan/Pvalues.txt")
fdata=read.table("Inputs/2_virscan/feature_data.txt")
metadata = read.table("Inputs/2_virscan/metadata.txt", header=TRUE)
types=read.table("Inputs/2_virscan/endemic_african_virus_types.txt",header=TRUE)

### This is necessary, and got lost when printing (printing does not accept colnames repeated)
old_colnames=colnames(pvals)
colnames(pvals)=metadata$Genotyping_ID

####################################################################################
####  2. Create table with counts of positive epitopes per individual per virus ####
####################################################################################

threshold=-log10(0.05)
counts=data.frame(matrix(ncol=length(unique(metadata$Genotyping_ID)),nrow=length(unique(fdata$Species))))
colnames(counts)=unique(metadata$Genotyping_ID)
rownames(counts)=unique(fdata$Species)

core_routine=function(i){
    virus=rownames(counts)[i]
    set_epitopes=fdata$Epitope[which(fdata$Species==virus)]
    values=pvals[which(rownames(pvals) %in% set_epitopes),which(colnames(pvals) %in% individual)]
    
    if(ncol(values)>1)
    {
        for(cols in 1:ncol(values))
        values[which(is.na(values[,cols])),cols]=0
        values=apply(values,1,min)
    }
    return(length(which(values>threshold)))
}

for(ind in c(1:ncol(counts)))
{
    print(ind)
    individual=colnames(counts)[ind]
    counts[,ind]=sapply(c(1:nrow(counts)),core_routine)
}

####################################################################################
####  3. Create species summary table: # of epitopes, mean of + epitopes
####################################################################################

species_summary=data.frame(Species=unique(fdata$Species),epitopes=NA)

for(ii in 1:nrow(species_summary))
{
    species_summary$epitopes[ii]=length(which(fdata$Species==species_summary$Species[ii]))
}

species_summary=species_summary[order(species_summary$Species),]
counts=counts[order(rownames(counts)),]

length(which(rownames(counts)!=species_summary$Species))

species_summary$means_counts=apply(counts,1,mean)

########################################################################################
####  4. Create matrix of relative deviation in seropositivity per individual and virus.
########################################################################################

relative_deviation_seropositivity=counts
for(ii in 1:nrow(counts))
{
    relative_deviation_seropositivity[ii,]=100*(counts[ii,]-species_summary$means_counts[ii])/species_summary$means_counts[ii]
}

########################################################################################
####  5. Model it against admixture.
########################################################################################

medias=apply(counts,1,mean)
relative_deviation_seropositivity=relative_deviation_seropositivity[which(medias>2),]
counts=counts[which(medias>2),]
medias=apply(counts,1,mean)

cols=metadata[which(!duplicated(metadata$Genotyping_ID)),]
design=model.matrix(~Admixture,data=cols)

fit <-lmFit(relative_deviation_seropositivity,design)
fit <- eBayes(fit)
fit_fdrs=data.frame(fit$p.value)
fit_fdrs[,2]=qvalue(fit$p.value[,2])$qvalue

betas=data.frame(fit$coefficients[order(rownames(fit$coefficients)),])
ps=data.frame(fit$p.value[order(rownames(fit$p.value)),])
fit_fdrs=fit_fdrs[order(rownames(fit_fdrs)),]

result=data.frame(logFC=betas$Admixture,p_value=ps$Admixture,Fdr=fit_fdrs$Admixture,Prevalence_media=medias)
rownames(result)=rownames(fit_fdrs)

########################################################################################
####  6. Do Volcano plot
########################################################################################


result$hit="A"
result$hit[which(result$Fdr<0.05 & result$logFC<0)]="Bakiga"
result$hit[which(result$Fdr<0.05 & result$logFC>0)]="Batwa"

result$Viral_species=rownames(result)
result$label=''
result$label[which(result$Fdr<0.05)]=result$Viral_species[which(result$Fdr<0.05)]
result=result[order(result$p_value),]
result$num_label=""
hits_indexes=which(result$hit!="A")
result$num_label[hits_indexes]=c(1:length(hits_indexes))


rownames(types)=types$Species
types=types[which(types$Species %in% rownames(fit_fdrs)),]
types$Species=factor(types$Species)
result_final=merge(result,types,by=0,all=FALSE)
rownames(result_final)=result_final$Row.names
result_final=result_final[order(result_final$p_value),c(11,2,3,4,6,9)]

rownames(result_final)[36]="SARS_related_coronavirus"
rownames(result_final)[19]="MERS_coronavirus"
result_final$hit=as.character(result_final$hit)
result_final$hit[which(result_final$hit=="A")]=""

#non, bakiga, batwa
colores=c("grey","#FF959F","#F21243")

vo=ggplot(result_final)  +
theme_classic()+
scale_fill_manual(values = colores) +
geom_text_repel(aes(x=logFC, y=-log10(p_value), label=num_label),force=1,size = 3, segment.color = 'grey50')+
geom_point(aes(x=logFC,y=-log10(p_value),fill=factor(hit), size=-log10(Fdr),alpha=1),stroke=1,color="white",shape=21)+
theme(
#legend.position="none",
axis.text.y   = element_text(size=14),
axis.title.y   = element_text(size=18),
axis.text.x   = element_text(size=14),
axis.title.x  = element_text(size=18),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.ticks.y = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1.5,linetype="solid"),
legend.title=element_text(size=22),
legend.text=element_text(size=24)
)+xlab("Relative variation in seropositivity (%)")+ylab("-log10(P value)")

system("mkdir -p Outputs/2_virscan")

pdf(paste0("Outputs/2_virscan/Volcano_Plot.pdf"),width=9,height=5)
print(vo)
dev.off()

write.table(result_final,paste0("Outputs/2_virscan/result.txt"))

############################################################################################
####  7. Get enrichment of double-stranded DNA viruses among virus more prevalent in Batwa
############################################################################################


DNAs=c("dsDNA")

DNA_assoc=length(which((result_final$Type %in% DNAs) & (result_final$Fdr<=0.05 & result_final$logFC>0)))
other_assoc=length(which(!(result_final$Type %in% DNAs) & (result_final$Fdr<=0.05 & result_final$logFC>0)))
DNA_no_assoc=length(which((result_final$Type %in% DNAs) & !(result_final$Fdr<=0.05 & result_final$logFC>0)))
other_no_assoc=length(which(!(result_final$Type %in% DNAs) & !(result_final$Fdr<=0.05 & result_final$logFC>0)))

mat=matrix(c(DNA_assoc,DNA_no_assoc,other_assoc,other_no_assoc),
nrow = 2,
dimnames = list(D_prev = c("Associations", "No_associations"),
Virus_type = c("DNA", "Other")))
fisher.test(mat, alternative = "two.sided")


## p-value = 0.006148
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
## 1.335128 8.663110
## sample estimates:
## odds ratio
## 3.331524




