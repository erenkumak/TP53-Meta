library(DESeq2)
Deseq2 <- function(experiment_title, control_group, treatment_group, file1 = NULL, file2 = NULL) {
  
  if(experiment_title=="Mutant data first batch"){ #changed
    data <- read.csv("data_files/mutant_data_first_batch.csv", header=TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("EV", "WT","WT","R175","R175","R273","R273"))
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    # Create a data set
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  if(experiment_title=="Mutant data combatseq"){ #changed combatseq
    data <- read.csv("data_files/mutant_data_combatseq.csv", header=TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    print(colnames(data))
    group <- factor(c("EV", "WT","WT","R175","R175","R273","R273", "WT", "EV", "R248W", "R248W"))
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    # Create a data set
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  if(experiment_title=="Mutant data second batch"){ #changed
    data <- read.csv("data_files/mutant_data_second_batch.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    colnames(data)=c("WT","EV","R248W_1","R248W_2")
    group <- factor(c("WT","EV" ,rep("R248W", 2)))
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  if(experiment_title=="Doxo First batch"){ #changed
    data <- read.csv("data_files/Doxo_data_firstBatch.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("EV_D_A5","EV_D","EV_A5","EV_A5" ,  
                      "EV","EV","R175_D_A5", "R175_D","R175_A5",
                      "R175", "R175" ,"R273_D_A5","R273_D_A5","R273_A5",
                      "R273_A5", "R273" ,"R273"))
    
    colnames(data)=c("EV_Doxo_A5","EV_Doxo_sc","EV_Serum_A5_1","EV_Serum_A5_2" ,  
                     "EV_Serum_sc_1","EV_Serum_sc_2","R175H_Doxo_A5", "R175H_Doxo_sc","R175H_Serum_A5",
                     "R175H_Serum_sc_1", "R175H_Serum_sc_2" ,"R273H_Doxo_A5_1","R273H_Doxo_A5_2","R273H_Serum_A5_1",
                     "R273H_Serum_A5_2", "R273H_Serum_sc_1" ,"R273H_Serum_sc_2")
  
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  if(experiment_title=="GSE301555"){ #changed
    data <- read.csv("data_files/siA5_siTP53_data.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("sc20",	"sc20",	"siA5",	"siA5",	"siP53",	"siP53",	"siA5_siP53",	"siA5_siP53"))
    
    colnames(data)=c("sc20_a",	"sc20_b",	"siA5_a",	"siA5_b",	"siP53_a",	"siP53_b",	"siA5_siP53_a",	"siA5_siP53_b")
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  if(experiment_title=="Doxo second batch"){ #changed
    data <- read.csv("data_files/Doxo_data_secondBatch.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("WT" ,"R175_A5" ,"WT_D"  ,  "WT_A5", "WT_D_A5" ,   "WT_A5", "WT"))
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  
  if(experiment_title=="Doxo data combatseq"){ #changed combatseq
    data <- read.csv("data_files/doxo_data_combatseq.csv", header=TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    print(colnames(data))
    group <- factor(c("EV_D_A5","EV_D","EV_A5","EV_A5" ,  
                      "EV","EV","R175_D_A5", "R175_D","R175_A5",
                      "R175", "R175" ,"R273_D_A5","R273_D_A5","R273_A5",
                      "R273_A5", "R273" ,"R273","WT" ,"R175_A5" ,"WT_D"  ,  "WT_A5", "WT_D_A5" ,   "WT_A5", "WT"))
    colnames(data)=c("EV_Doxo_A5","EV_Doxo_sc","EV_Serum_A5_1","EV_Serum_A5_2" ,  
                     "EV_Serum_sc_1","EV_Serum_sc_2","R175H_Doxo_A5", "R175H_Doxo_sc","R175H_Serum_A5",
                     "R175H_Serum_sc_1", "R175H_Serum_sc_2" ,"R273H_Doxo_A5_1","R273H_Doxo_A5_2","R273H_Serum_A5_1",
                     "R273H_Serum_A5_2", "R273H_Serum_sc_1" ,"R273H_Serum_sc_2", "WT_serum_Sc_2",                       
                     "R175_serum_A5","WT_doxo_Sc","WT_serum_A5_1" ,                      
                     "WT_doxo_A5", "WT_serum_A5_2","WT_serum_Sc_1")
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    # Create a data set
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
    
  }
  
  
  if(experiment_title=="GSE47042_Nutlin"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE47042_Nutlin.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("Control" ,"Control","Nutlin" , "Nutlin" ))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE173364_TP53_MUT"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE173364_TP53_MUT.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("TP53_KO","TP53_KO_R175H.OE",
                      "TP53_KO_R273H.OE","shp53","TP53_KO",                                 
                      "TP53_KO_R175H.OE","TP53_KO_R273H.OE","shp53",                                
                      "shp53","shp53_R175H.OE","shp53_R175H.OE",
                      "TP53_KO","TP53_KO_R175H.OE_","TP53_KO_R273H.OE",          
                      "shp53_R175H.OE","TP53_KO_six.hours.post.seeding","TP53_KO_six.hours.post.seeding"    ,           
                      "TP53_KO_R273H.OE_6h.post.seeding" ,"shp53_R273H.OE","shp53_R273H.OE",                  
                      "shp53_R273H.OE","TP53_KO_R273H.OE_6h.post.seeding"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE68248_siTP53"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE68248_siTP53.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("siTP53",  "siTP53", "siTP53","siCtrl","siCtrl" ,"siCtrl"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE139546"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE139546.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    colnames(data) <- c("X","CTR_siRNA","CTR_siRNA","p53_siRNA","p53_siRNA","iASPP_siRNA",
                        "iASPP_siRNA","p53_iASPP_siRNA","p53_iASPP_siRNA","CTRsi_Nutlin","CTRsi_Nutlin",
                        "p53si_Nutlin","p53si_Nutlin","iASPPsi_Nutlin","iASPPsi_Nutlin","p53si_iASPPsi_Nutlin","p53si_iASPPsi_Nutlin")
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("CTR_siRNA","CTR_siRNA","p53_siRNA","p53_siRNA","iASPP_siRNA",
                      "iASPP_siRNA","p53_iASPP_siRNA","p53_iASPP_siRNA","CTRsi_Nutlin","CTRsi_Nutlin",
                      "p53si_Nutlin","p53si_Nutlin","iASPPsi_Nutlin","iASPPsi_Nutlin","p53si_iASPPsi_Nutlin","p53si_iASPPsi_Nutlin"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE83635"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE83635.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("Control","Control","Control","siSET","siSET","siSET","CRISPR_p53","CRISPR_p53","CRISPR_p53_siSET","CRISPR_p53_siSET"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE84877"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE84877.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("Scrambled","Scrambled","Scrambled","si_p53","si_p53","si_p53","si_PKM2","si_PKM2","si_PKM2","siPKM2_sip53","siPKM2_sip53","siPKM2_sip53"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE74493"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE74493.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("siControl","siControl","siLATS1_2","siLATS1_2","sip53","sip53","siLATS1_2_sip53","siLATS1_2_sip53"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE47042"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE47042.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("Control","Control","Nutlin","Nutlin"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE124508"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE124508.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("DMSO","DMSO","DMSO","Idasanutlin","Idasanutlin","Idasanutlin"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE59841"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE59841.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("DMSO","Nutlin_low","Nutlin_high","DMSO","Nutlin_high","DMSO","Nutlin_low","Nutlin_high"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE86219"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE86219.csv", header=TRUE)
    data$X <- NULL
    data[,c(2,3,8)] <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("DMSO","DMSO","Nutlin","Nutlin"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE248670"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE248670.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("WT_DMSO","WT_DMSO","WT_DMSO","WT_Nutlin3A","WT_Nutlin3A","WT_Nutlin3A","GDS3_overexp_DMSO","GDS3_overexp_DMSO","GDS3_overexp_DMSO",
                      "GDS3_overexp_Nutlin3A","GDS3_overexp_Nutlin3A","GDS3_overexp_Nutlin3A"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="GSE89110"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <-  read.csv("data_files/GSE_new/GSE89110.csv", header=TRUE)
    data$X <- NULL
    data <- data[!duplicated(data$SYMBOL),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c("p53WT_shCtrl","p53WT_shCtrl","p53WT_shCtrl","p53WT_shPFKFB4","p53WT_shPFKFB4","p53WT_shPFKFB4",
                      "p53NULL_shCtrl","p53NULL_shCtrl","p53NULL_shCtrl","p53NULL_shPFKFB4","p53NULL_shPFKFB4","p53NULL_shPFKFB4"))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  if(experiment_title=="Your_Data"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    data <-  file1
    cond <- file2
    data <- data[!duplicated(data$Gene),]
    data <- na.omit(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    data <- data[rowSums(data != 0) > 0, ]
    group <- factor(c(cond$Condition))
    
    
    
    design=data.frame(row.names = c(colnames(data)), condition=group)
    
    dds=DESeqDataSetFromMatrix(data, colData = design, design = ~condition)
    dds <- DESeq(dds)
    res <- data.frame(DESeq2::results(dds, contrast= c("condition", treatment_group,control_group),cooksCutoff=FALSE ))
    names <- paste(treatment_group,  "_vs_", control_group,sep = "")
    colnames(res) <- c( paste( names , c("baseMean", "logFC", "lfcSE","stat", "p value", "p adj value")))
    result_df <- res
  }
  return(result_df)
}
