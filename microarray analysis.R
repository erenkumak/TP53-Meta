GSE160497

library(GEOquery)
gsm <- getGEO('GSE160497')



gse <- gsm[[1]]
gse
pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data
sampleInfo <- pData(gse)
sampleInfo
colnames(sampleInfo)
library(dplyr)
sampleInfo <- sampleInfo[,c("title", "geo_accession","source_name_ch1")]

colnames(sampleInfo)[3]="group"
library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write.csv(full_output, "gse_full_output.csv",row.names = F)
library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
## the column names are a bit ugly, so we will rename
colnames(design) <- sort(c("Dox","minusDox","Dox_Aldo","minusDox_Aldo","Dox_E2","Dox_E2_Aldo","Dox_RA","Dox_RA_Aldo"))
saveRDS(sampleInfo,"sampleInfo.RDS")
summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse1 <- gse[keep,]
saveRDS(gse1,"gse1.RDS")

data=exprs(gse1) # to make it normal csv file as we used in rnaseq
anno <- fData(gse1)
anno <- dplyr::select(anno,Symbol)
data2=merge(anno,data,by=0)

data3=data2[data2$Symbol!="",]
class(data3$Symbol)
data3=data3[,-1]
data3 <- data3 %>% group_by(Symbol) %>% dplyr::summarize_all((mean))
rownames(data3) <- data3$Symbol
data3=data.frame(data3)
write.csv(data3,"GSE160497_Array_filtered.csv",row.names = F)


fit <- lmFit(data3[,-1], design)
head(fit$coefficients)
contrasts <- makeContrasts(Dox_Aldo - Dox, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)


anno <- fData(gse1)
anno
anno <- dplyr::select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)


