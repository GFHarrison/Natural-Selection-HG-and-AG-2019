##############################################################
#### 1. Load dependencies ####################################
##############################################################

### Note: declare condition equal to "CTL, "LPS", or "GARD" to obtain cis-EQTLs for each condition.
condition="CTL"
## Number of random permutations (for empiric fdr corrections).
iterations=10

library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)

current=getwd()
setwd("Github_codes/common_functions")
source("permFDR.R")
setwd(current)

##########################################################
#### 2. Load input files & clean/declare output files ####
##########################################################

genepos = read.table("Inputs/3_EQTL_mapping/gene_positions.txt",header = TRUE, stringsAsFactors = FALSE)
snpspos = read.table(paste0("Inputs/3_EQTL_mapping/SNP_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
gtypes = read.table(paste0("Inputs/3_EQTL_mapping/genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)

system(paste0("mkdir -p Outputs/3_EQTL_mapping/temp_files/",condition))
system(paste0("mkdir -p Outputs/3_EQTL_mapping/",condition,"/raw_results"))

## Erase and create again the temp file.
system(paste0("rm -rf Outputs/3_EQTL_mapping/temp_files/"))
system(paste0("mkdir -p Outputs/3_EQTL_mapping/temp_files/",condition))

############################################################################################################
#### 3. Obtain clean input files:                                                                       ####
####   - Voomed Expression tables per condition, where first n principal components are regressed out.  ####
####   - Covariate tables for matrixEQTL, including:                                                    ####
####     2 first genotypes PCs, Flowcell,Sex,fraction of reads assigned and tissue composition          ####
############################################################################################################

#####
##### Clean global tables of expression and metadata
#####

metadata_whole = read.table(paste0("Inputs/1_DE_analyses/RNAseq_metadata.txt"))
reads_whole = read.table("Inputs/1_DE_analyses/RNAseq_reads_matrix.txt")

## Select only samples for which genotype data is available in the global metadata table, and order samples:
metadata_whole=metadata_whole[which(metadata_whole$EQTL_set==1),]
metadata_whole=metadata_whole[order(rownames(metadata_whole)),]

## Subset the corresponding columns in the reads matrix, and order samples and genes.
reads_whole=reads_whole[,which(colnames(reads_whole) %in% rownames(metadata_whole))]
reads_whole=reads_whole[order(rownames(reads_whole)),order(colnames(reads_whole))]

## Check that order of elements in the metadata_whole and reads_whole now are congruent, samples-wise
length(which(rownames(metadata_whole)!=colnames(reads_whole)))

## Obtain log(cpm), voomed expression data
design = model.matrix(~Condition +Condition:Admixture, data=metadata_whole)
dge <- DGEList(counts=reads_whole)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=FALSE)
voomed_reads_whole = as.data.frame(v$E)

#####
##### Filter expression and metadata per condition.
#####

## Filter metadata per condition
metadata=metadata_whole[which(metadata_whole$Condition==condition),]

## Clean factor variables, and mean center numeric ones.
metadata$Condition=factor(metadata$Condition)
metadata$Genotyping_ID=factor(metadata$Genotyping_ID)
metadata$Flowcell=factor(metadata$Flowcell)
metadata$CD14=metadata$CD14-mean(metadata$CD14)
metadata$CD4=metadata$CD4-mean(metadata$CD4)
metadata$CD20=metadata$CD20-mean(metadata$CD20)
metadata$fraction_assigned=metadata$fraction_assigned-mean(metadata$fraction_assigned)

## Filter voomed expression per condition
voomed_reads=voomed_reads_whole[,which(colnames(voomed_reads_whole) %in% rownames(metadata))]

## Check again the con=herence of samples order
length(which(colnames(voomed_reads)!=rownames(metadata)))

### Shift from sampleIDs to Genotyping_IDs (there only will be one sample per genotype in every analysis).
colnames(voomed_reads)=metadata$Genotyping_ID
rownames(metadata)=metadata$Genotyping_ID

### Recover alphabetical order with the new IDs.
voomed_reads=voomed_reads[,order(colnames(voomed_reads))]
metadata=metadata[order(metadata$Genotyping_ID),]
length(which(colnames(voomed_reads)!=rownames(metadata)))

###
### Build matrixEQTL input: expression tables: regressing out first n PCs from voomed_reads:
###

## We will remove 8,8 and 11 PCs for CTL, GARD and LPS, respectively
if(condition=="LPS")
{
    pc_set=c(1:11)
}else{
    pc_set=c(1:8)}

### Regress those out.
pca_rm <- function(input_data, pc_set) {
    pca = prcomp(t(input_data), na.action = na.omit)
    new = input_data
    new = apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
    new = t(new)
    colnames(new) = colnames(input_data)
    rownames(new) = rownames(input_data)
    return(new)
}
expression = pca_rm(voomed_reads, pc_set)

###
### Build matrixEQTL input: covariates tables: transpose metadata & add genotype 1st and 2nd PCs:
###

### Get genotype PC analysis and clean genotype data

metadata_individuals=metadata_whole[which(!duplicated(metadata_whole$Genotyping_ID)),]
metadata_individuals=metadata_individuals[order(metadata_individuals$Genotyping_ID),]
rownames(metadata_individuals)=metadata_individuals$Genotyping_ID

samples=colnames(gtypes)
gtypes_pca=data.frame(snp_id=rownames(gtypes),gtypes)

snpgdsCreateGeno("Outputs/3_EQTL_mapping/temp_files/GDS_genotypes.gds",
genmat = as.matrix(gtypes_pca[, samples]),
sample.id = unique(samples),
snp.id = gtypes_pca$snp_id,
snpfirstdim=TRUE)

snpgdsSummary("Outputs/3_EQTL_mapping/temp_files/GDS_genotypes.gds")

genofile <- snpgdsOpen("Outputs/3_EQTL_mapping/temp_files/GDS_genotypes.gds")

pca <- snpgdsPCA(genofile)
tab <- data.frame(sample.id = pca$sample.id,
PC1 = pca$eigenvect[,1],    # the first eigenvector
PC2 = pca$eigenvect[,2],
PC3 = pca$eigenvect[,3],
PC4 = pca$eigenvect[,4],
PC5 = pca$eigenvect[,5],
stringsAsFactors = FALSE)

### Declare covariates table.

pcs_genotypes=tab[which(tab$sample.id %in% rownames(metadata)),]
pcs_genotypes=pcs_genotypes[order(pcs_genotypes$sample.id),]
length(which(rownames(metadata) !=pcs_genotypes$sample.id))
metadata$PC1=pcs_genotypes$PC1
metadata$PC2=pcs_genotypes$PC2

covariates=t(model.matrix(~PC1+PC2+Flowcell+Sex+fraction_assigned+CD14+CD4+CD20,data=metadata))
covariates=covariates[2:nrow(covariates),]

expression=expression[which(rownames(expression) %in% genepos$Gene_ID),]

##################################################################################################
####  4. Before calling matrixEQTL subset individuals present in the condition to analyze and ####
####     remove genes and SNPs from genotype and expression tables                            ####
####     for which there is no available position (These wouldn't be tested for cis-EQTL,     ####
####     but had useful info to include in PC analyses performed above)                       ####
##################################################################################################

### Subset the genotypes file with the individuals present in each condition
### and the SNPs present in the SNP positions file.
genotypes=gtypes[,which(colnames(gtypes) %in% colnames(covariates))]
genotypes=genotypes[which(rownames(genotypes) %in% snpspos$snp),]
genotypes=genotypes[,order(colnames(genotypes))]
length(which(rownames(genotypes)!=snpspos$snp))

###############################################
#### 5. Check input data files congruence. ####
###############################################

## Samples-wise
length(which(rownames(metadata)!=colnames(covariates)))
length(which(rownames(metadata)!=colnames(expression)))
length(which(rownames(metadata)!=colnames(genotypes)))

## SNPs-wise
length(which(rownames(genotypes)!=snpspos$snp))

## Gene-wise
length(which(rownames(expression)!=genepos$Gene_ID))
## All 0, everything is coherent.

##############################################
#### 6. save matrixEQTL temp input files. ####
##############################################


snps_positions_file_name="Inputs/3_EQTL_mapping/SNP_positions.txt"
gene_positions_file_name="Inputs/3_EQTL_mapping/gene_positions.txt"
expression_file_name=paste0("Outputs/3_EQTL_mapping/temp_files/",condition,"/expression.txt")
covariates_file_name=paste0("Outputs/3_EQTL_mapping/temp_files/",condition,"/covariates.txt")
SNP_file_name=paste0("Outputs/3_EQTL_mapping/temp_files/",condition,"/genotypes.txt")

write.table(genotypes,SNP_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(expression,expression_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(covariates,covariates_file_name, quote=F, sep="\t", row.names=TRUE)

### In this loop iter=0 runs the actual analyses, iters 1 to iterations (10), run permutation instances for
### FDR correction

permuted_pvalues_folder=paste0("Outputs/3_EQTL_mapping/",condition,"/raw_results/")
for(iteration in 0:iterations){
    
    ##############################################################
    #### 7. Permute genotype data (only for iterations>0 #########
    ##############################################################
    
    if(iteration>0){
        
        cols<-colnames(genotypes)
        cols.perm<-sample(cols)
        if(iteration==1){
        random_individuals_df=data.frame(cols.perm)
        }else{
            random_individuals_df=cbind(random_individuals_df,cols.perm)
        }
        genotypes<-genotypes[,cols.perm]
        colnames(genotypes)<-cols
        write.table(genotypes,SNP_file_name, sep="\t", quote = FALSE)
    }
    
    ##############################################################
    #### 8. Prepare & Run Matrix EQTL ############################
    ##############################################################
    
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";     # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## Load covariates
    
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    
    ## Load genotype data
    
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA";
    snps$fileOmitCharacters = "-9" # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    useModel = modelLINEAR
    output_file_name_cis = tempfile()
    pvOutputThreshold_cis = 1
    pvOutputThreshold = 0;
    errorCovariance = numeric()
    cisDist = 1e5
    output_file_name = tempfile()
    output_file_name_cis = tempfile()
    
    me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
    
    ##############################################################
    #### 9. Write temporal output files ##########################
    ##############################################################
    
    unlink(output_file_name_cis);
    
    if(iteration==0){
        write.table(me$cis$eqtls, file = paste0("Outputs/3_EQTL_mapping/",condition,"/raw_results/result_original.txt"))}else{
            write.table(me$cis$eqtls, file = paste0("Outputs/3_EQTL_mapping/",condition,"/raw_results/result_permuted_",iteration,".txt"))}
        
}

##############################################################################################################
#### 10. Write Best SNP-gene associations files for Delta-PVE analyses & prepare files for FDR corrections ###
##############################################################################################################

## Selecting the top cis-SNP for each gene in true and permuted files.

for(iteration in 0:iterations)
{
    
    if(iteration==0){
        event=read.table(paste0("Outputs/3_EQTL_mapping/",condition,"/raw_results/result_original.txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,5]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,5]),]
    }else{
        event=read.table(paste0(permuted_pvalues_folder,"/result_permuted_",iteration,".txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,5]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,5]),]
    }
    if(iteration==0){
        original_best_EQTL=event.bestQTL
    }else{
        if(iteration==1)
        {
            Permutation_Input = event.bestQTL[4]
        }else{
            Permutation_Input=cbind(Permutation_Input,event.bestQTL[4])}
    }
}

##################################################
#### 11. Run FDR corrections and write outputs ###
##################################################

cis_eQTL_fdrs = permFDR(full_data= original_best_EQTL, full_column_id = "pvalue", perm_data= Permutation_Input, perm_column_ids="all", output_name= paste0("Outputs/3_EQTL_mapping/",condition))

write.table(cis_eQTL_fdrs$fdrs,file = paste0("Outputs/3_EQTL_mapping/",condition,"/results_best_SNPs.txt"), quote=FALSE)





