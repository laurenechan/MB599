# R script for differential gene expression analysis utilizing a linear model with quantile transformation and sample-quality weight normalization

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(biomaRt)
#library(gplots)
#library(RColorBrewer)

# DATA ORGANIZATION 

#read in counts and create DGElist object
counts <- read.table("~/RNAseq/all-results.tsv", sep='\t', header=TRUE, stringsAsFactors = FALSE, row.names=1) #row.names required for import into DGEList as matrix
counts
as.data.frame(counts) # just to confirm that it is a dataframe, edgeR will read in as matrix

#this next bit is me embarrassingly renaming every sample name
colnames(counts)[colnames(counts) == "sample_1_lane7"] <- "plus_12_7_1"
colnames(counts)[colnames(counts) == "sample_2_lane7"] <- "plus_12_7_2"
colnames(counts)[colnames(counts) == "sample_3_lane8"] <- "plus_12_8_3"
colnames(counts)[colnames(counts) == "sample_4_lane8"] <- "plus_12_8_4"
colnames(counts)[colnames(counts) == "sample_5_lane7"] <- "minus_12_7_1"
colnames(counts)[colnames(counts) == "sample_6_lane7"] <- "minus_12_7_2"
colnames(counts)[colnames(counts) == "sample_7_lane8"] <- "minus_12_8_3"
colnames(counts)[colnames(counts) == "sample_8_lane8"] <- "minus_12_8_4"
colnames(counts)[colnames(counts) == "sample_9_lane7"] <- "plus_18_7_1"
colnames(counts)[colnames(counts) == "sample_10_lane7"] <- "plus_18_7_2"
colnames(counts)[colnames(counts) == "sample_11_lane8"] <- "plus_18_8_3"
colnames(counts)[colnames(counts) == "sample_12_lane8"] <- "plus_18_8_4"
colnames(counts)[colnames(counts) == "sample_13_lane7"] <- "minus_18_7_1"
colnames(counts)[colnames(counts) == "sample_14_lane7"] <- "minus_18_7_2"
colnames(counts)[colnames(counts) == "sample_15_lane8"] <- "minus_18_8_3"
colnames(counts)[colnames(counts) == "sample_16_lane8"] <- "minus_18_8_4"
colnames(counts)[colnames(counts) == "sample_17_lane7"] <- "plus_24_7_1"
colnames(counts)[colnames(counts) == "sample_18_lane7"] <- "plus_24_7_2"
colnames(counts)[colnames(counts) == "sample_19_lane8"] <- "plus_24_8_3"
colnames(counts)[colnames(counts) == "sample_20_lane8"] <- "plus_24_8_4"
colnames(counts)[colnames(counts) == "sample_21_lane7"] <- "minus_24_7_1"
colnames(counts)[colnames(counts) == "sample_22_lane7"] <- "minus_24_7_2"
colnames(counts)[colnames(counts) == "sample_23_lane8"] <- "minus_24_8_3"
colnames(counts)[colnames(counts) == "sample_24_lane8"] <- "minus_24_8_4"
colnames(counts)

# create the DGEList list-based data object
dge <- DGEList(counts)
dge 
dim(dge) # numbers represent transcripts/genes and number of samples, 32254 and 24
colnames(dge) # checking to see if names changed like I intended
# create the 3 variables that we will test for or model our data with diet, time and sequencing lane
diet <- factor(c("plus", "plus", "plus", "plus", "minus", "minus", "minus", "minus", "plus", "plus", "plus", "plus", "minus", "minus", "minus", "minus", "plus", "plus", "plus", "plus", "minus", "minus", "minus", "minus"))
time <- factor(c("12", "12", "12", "12", "12", "12", "12", "12", "18", "18", "18", "18", "18", "18", "18", "18", "24", "24", "24", "24", "24", "24", "24", "24"))
lane <- factor(c("7", "7", "8", "8", "7", "7", "8", "8", "7", "7", "8", "8", "7", "7", "8", "8", "7", "7", "8", "8", "7", "7", "8", "8"))
dge$samples$diet <- diet
dge$samples$time <- time
dge$samples$lane <- lane
dge$samples # confirm that variables correctly apply to samples
dge$samples$group <- interaction(diet, time) # define the interaction we care about most and will test for
dge <- DGEList(counts, group=group) # apply that interaction to our DGEList

# use biomaRt to load zebrafish genome
mart=useMart(dataset="drerio_gene_ensembl",biomart = "ensembl")
# we have ensembl ids output from STAR, want entrezid for TopGO/KEGG and zfin for legibility 
symbols <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "zfin_id_symbol"), 
                 filter = "ensembl_gene_id", values = geneid, mart = ensembl)
head(symbols) # the list provided is for some reason in numerical order, not by transcript list you provide
geneid <- rownames(dge) # so you have to pull out the names
ord <- match(unlist(lapply(strsplit(geneid, split = ".", fixed = T), function(x)x[1])), symbols[,"ensembl_gene_id"]) # and match your newfound gene names to what is in your list
dge$genes <- symbols[ord,"entrezgene_id"] #then add the new entrezgene_id symbol to a new array $genes in the DGEList
dge # $genes should carry through all analysis ahead

# DATA EXPLORATION BEFORE ANALYSIS

# STEP 1, FILTER AND NORMALIZE
#filer out low expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, ,keep.lib.sizes=FALSE] # original gene set list was 32524, filtered down to 17368
# sometimes I get the error incorrect number of dimensions which only happens after I add the symbols to the DGEList
# erorr does not change ability to complete analysis

#TMM normalization
dge <- calcNormFactors(dge, method="TMM") # makes geometric mean of each effective library size the same, scaled up or down based on lib.size seen in dge$samples
dge$samples

#STEP 2, CLUSTER BY LOGFC SHOW SAMPLES AND LANE(BATCH) EFFECT
#lcpm <- cpm(dge, log=TRUE)
#par(mfrow=c(1,2))
#col.group <- group
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#col.group <- as.character(col.group)
#col.lane <- lane
#levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
#col.lane <- as.character(col.lane)
#plotMDS(lcpm, labels=group, col=col.group)
#title(main="A. Sample groups")
#plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
#title(main="B. Sequencing lanes")

# STEP 3, REMOVAL OF LOW EXPRESSED GENES BY DEFIED MEAN AND MEDIAN CUTOFF
#cpm <- cpm(dge)
#lcpm <- cpm(dge, log=TRUE)
#L <- mean(dge$samples$lib.size) * 1e-6
#M <- median(dge$samples$lib.size) * 1e-6
#c(L, M)
#lcpm.cutoff <- log2(10/M + 2/L)
#nsamples <- ncol(dge)
#col <- brewer.pal(nsamples, "Paired")
#par(mfrow=c(1,2))
#plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
#title(main="A. Raw data", xlab="Log-cpm")
#abline(v=lcpm.cutoff, lty=3)
#for (i in 2:nsamples){
#  den <- density(lcpm[,i])
#  lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", samplenames, text.col=col, bty="n")
#lcpm <- cpm(dge, log=TRUE)
#plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
#title(main="B. Filtered data", xlab="Log-cpm")
#abline(v=lcpm.cutoff, lty=3)
#for (i in 2:nsamples){
#  den <- density(lcpm[,i])
#  lines(den$x, den$y, col=col[i], lwd=2)
#}
#legend("topright", samplenames, text.col=col, bty="n")

# STEP 4, EXPRESSION DISTRIBUTIONS BEFORE AND AFTER NORMALIZATION BY TMM
#dge2 <- dge
#dge2$samples$norm.factors <- 1
#dge2$counts[,1] <- ceiling(dge2$counts[,1]*0.05)
#dge2$counts[,2] <- x2$counts[,2]*5

#par(mfrow=c(1,2))
#lcpm <- cpm(dge2, log=TRUE)
#boxplot(lcpm, las=2, col=col, main="")
#title(main="A. Example: Unnormalised data", ylab="Log-cpm")
#dge2 <- calcNormFactors(dge2)
#dge2$samples$norm.factors

#lcpm <- cpm(dge2, log=TRUE)
#boxplot(lcpm, las=2, col=col, main="")
#title(main="B. Example: Normalised data", ylab="Log-cpm")

# DIFFERENTIAL GENE EXPRESSION ANALYSIS

# START BY MAKING A DESIGN MATRIX
design <- model.matrix(~0 + group) # ~0 removes the intercept, compares all groups to each other
colnames(design) <- gsub("group", "", colnames(design)) #design matrix adds extra word group, this removes it
design # confirm comparisons are correct
colnames(design) # confirm group names are correct

# MAKE A CONTRAST MATRIX OF DESIRED CONTRASTS TO TEST
timecontr.matrix <- makeContrasts(
  PlusvMinus_12h = plus.12 - minus.12,
  PlusvMinus_18h = plus.18 - minus.18, 
  PlusvMinus_24h = plus.24 - minus.24,
  levels = colnames(design))
timecontr.matrix # originally thought this comparison was desired, can be used to test all group interactions at once
contrast12 <- makeContrasts(minus.12 - plus.12, levels=design) # decided to include individual contrasts between diet groups
contrast18 <- makeContrasts(minus.18 - plus.18, levels=design)
contrast24 <- makeContrasts(minus.24 - plus.24, levels=design)

# ANALYSIS BY VOOM TREND WITH SAMPLE QUALITY WEIGHT AND TMM NORMALIZATION
vwts <- voomWithQualityWeights(dge, design, normalization="TMM", plot=TRUE) 
# Each block is for each contrast to measure with individual outputs
# 12 hpf contrast, include in output the entrezgene id and zfin id
vfit <- lmFit(vwts)
vfit <- contrasts.fit(vfit, contrast12)
vfit <- eBayes(vfit)
tmp <- topTable(vfit,sort.by="P", n=Inf)
tmp$Gene <- rownames(tmp)
ord <- match(unlist(lapply(strsplit(tmp$Gene, split = ".", fixed = T), function(x)x[1])), symbols[,"ensembl_gene_id"])
tmp$NCBI <- symbols[ord,"entrezgene_id"]
tmp$ZFIN <- symbols[ord,"zfin_id_symbol"]
tmp <- tmp[,c("Gene","NCBI", "ZFIN","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp$adj.P.Val < 0.05))
write.table(tmp, file = "12hpf_DEG_VoomWQW.txt", row.names = F, sep = "\t", quote = F)

# 18 hpf contrast, include in output the entrezgene id and zfin id
vfit <- lmFit(vwts)
vfit <- contrasts.fit(vfit, contrast18)
vfit <- eBayes(vfit)
tmp <- topTable(vfit,sort.by="P", n=Inf)
tmp$Gene <- rownames(tmp)
ord <- match(unlist(lapply(strsplit(tmp$Gene, split = ".", fixed = T), function(x)x[1])), symbols[,"ensembl_gene_id"])
tmp$NCBI <- symbols[ord,"entrezgene_id"]
tmp$ZFIN <- symbols[ord,"zfin_id_symbol"]
tmp <- tmp[,c("Gene","NCBI", "ZFIN","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp$adj.P.Val < 0.05))
write.table(tmp, file = "18hpf_DEG_VoomWQW.txt", row.names = F, sep = "\t", quote = F)

# 24 hpf contrast, include in output the entrezgene id and zfin id
vfit <- lmFit(vwts)
vfit <- contrasts.fit(vfit, contrast24)
vfit <- eBayes(vfit)
tmp <- topTable(vfit,sort.by="P", n=Inf)
tmp$Gene <- rownames(tmp)
ord <- match(unlist(lapply(strsplit(tmp$Gene, split = ".", fixed = T), function(x)x[1])), symbols[,"ensembl_gene_id"])
tmp$NCBI <- symbols[ord,"entrezgene_id"]
tmp$ZFIN <- symbols[ord,"zfin_id_symbol"]
tmp <- tmp[,c("Gene","NCBI", "ZFIN","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp$adj.P.Val < 0.05))
write.table(tmp, file = "24hpf_DEG_VoomWQW.txt", row.names = F, sep = "\t", quote = F)

# all time contrast matrix, include in output the entrezgene id and zfin id
vfit <- lmFit(vwts)
vfit <- contrasts.fit(vfit, timecontr.matrix)
vfit <- eBayes(vfit)
tmp <- topTable(vfit,sort.by="P", n=Inf)
tmp$Gene <- rownames(tmp)
ord <- match(unlist(lapply(strsplit(tmp$Gene, split = ".", fixed = T), function(x)x[1])), symbols[,"ensembl_gene_id"])
tmp$NCBI <- symbols[ord,"entrezgene_id"]
tmp$ZFIN <- symbols[ord,"zfin_id_symbol"]
tmp <- tmp[,c("Gene","NCBI", "ZFIN","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp$adj.P.Val < 0.05))
write.table(tmp, file = "AllContrasts_DEG_VoomWQW.txt", row.names = F, sep = "\t", quote = F)

# ONTOLOGY & KEGG PATHWAY ANALYSIS
# I recommend performing the GO and KEGG analysis separately for each vfit for further independent testing and exploration
go <- goana.MArrayLM(vfit, geneid=vfit$genes, species="Dr")
kegg <- kegga.MArrayLM(vfit, geneid=vfit$genes, species="Dr")
topGO(go, sort="Up")
topGO(go, sort="Down")
topKEGG(kegg, sort = "Up")
topKEGG(kegg, sort = "Down")

> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.15.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] biomaRt_2.38.0 edgeR_3.24.3   limma_3.38.3  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.4.6         compiler_3.5.1       prettyunits_1.1.1    bitops_1.0-6         tools_3.5.1         
 [6] progress_1.2.2       statmod_1.4.34       digest_0.6.25        bit_1.1-15.2         RSQLite_2.2.0       
[11] memoise_1.1.0        lattice_0.20-41      pkgconfig_2.0.3      rlang_0.4.6          DBI_1.1.0           
[16] rstudioapi_0.11      curl_4.3             parallel_3.5.1       stringr_1.4.0        httr_1.4.1          
[21] S4Vectors_0.20.1     vctrs_0.3.0          IRanges_2.16.0       hms_0.5.3            locfit_1.5-9.4      
[26] stats4_3.5.1         bit64_0.9-7          grid_3.5.1           Biobase_2.42.0       R6_2.4.1            
[31] AnnotationDbi_1.44.0 XML_3.99-0.3         GO.db_3.7.0          blob_1.2.1           magrittr_1.5        
[36] org.Dr.eg.db_3.7.0   splines_3.5.1        BiocGenerics_0.28.0  stringi_1.4.6        RCurl_1.98-1.2      
[41] crayon_1.3.4        


