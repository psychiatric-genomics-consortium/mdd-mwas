# set seed so all randomly created distributions are the same
set.seed(1)

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0) {
	stop("Needs and argument\n");
} else {
	F_PDATA <- args[1]
}
cat(paste0("NOTE Using ", F_PDATA, " as input\n"))

check_file_die <- function(file) {
	if(!file.exists(file)){
		cat(paste0(file, " doesn't exist\n"))
		quit(save="no", status=1, runLast=FALSE)
	}
}

rename_file <- function(file) {
	if(file.exists(file)){
		file.rename(file, paste0(file, ".previous"))
	}
}

logit <- function(str) { cat(paste0(seq, "\t", date(), "\t", str, "\n"), file=NOTES, append=TRUE); seq <<- seq+1 }

############################################################

check_file_die(F_PDATA)
check_file_die(F_MVALS  <- "/disk/home/mberming/Konrad_corrected/output/ForAna_GRM_corrected_Mvalues.rds")
check_file_die(F_PCA    <- "/disk/home/mberming/Konrad_corrected/output/Konrad_Corrected_resid_mvals_relateds_PCA_301017.rds")
check_file_die(F_PROBES <- "/data/stradl/wave1/methylation-datasets/SNP_CH_probes")
NOTES                   <- "Summary-__TRAIT__.txt"
seq                     <- 1
rename_file(NOTES)
cat(paste0("n\tTime\tWhat\tRows\tCols\n"), file=NOTES, append=TRUE)

############################################################

pdata <- read.table(F_PDATA, sep=",", header=TRUE)
logit(paste0("pdata\t", nrow(pdata)))

# looking for outliers
jpeg("boxplot-__TRAIT__.jpeg")
length(boxplot(pdata$age, range=3)$out)
dev.off()

pdata<-pdata[is.na(pdata$__TRAIT__)==FALSE,]
logit(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))

table(pdata$__TRAIT__)

#load norm mvals object
system.time({mvals <- readRDS(F_MVALS)})[3]
logit(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

# Extracting the sample order from the mvals file. We will use this index to maintain this order in the phenotype and PC data frames.
ids <- data.frame(Sample_Sentrix_ID=colnames(mvals), N=c(1:length(colnames(mvals))))
logit(paste0("ids\t", nrow(ids), "\t", ncol(ids)))

#load PCA analysis of residualised M-values 
system.time({PCA <- readRDS(F_PCA)})[3]

#load list of probes affected by SNPs/that are predicted to cross-hybridise
probes_to_exclude <- read.table(F_PROBES, sep="\t", header=FALSE)
logit(paste0("probes\t", nrow(probes_to_exclude)))

pdata <- merge(pdata, ids, by="Sample_Sentrix_ID", all.x=TRUE)
logit(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))

##sorting the phenotype file
pdata <- pdata[order(pdata$N),]
logit(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))

mvals <- mvals[,which(colnames(mvals) %in% pdata$Sample_Sentrix_ID)]
logit(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

PCA <- as.data.frame(PCA$x)
PCA <- data.frame(PCA[ ,1:50])
PCA$Sample_Sentrix_ID <- row.names(PCA)
logit(paste0("PCA\t", nrow(PCA), "\t", ncol(PCA)))

##Merge PCs to phenotype data
pdata <- merge(pdata, PCA, by="Sample_Sentrix_ID", all.x=TRUE)
logit(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))
##sorting the phenotype file
pdata <- pdata[order(pdata$N),]
logit(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))

#Exclude SNP/cross-hybridising probes
mvals <- mvals[-which(rownames(mvals) %in% probes_to_exclude$V1),]
logit(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

#Remove rows containing NAs
row.has.na <- apply(mvals, 1, function(x){any(is.na(x))})
mvals <- mvals[-which(rownames(mvals) %in% names(which(row.has.na))),]
logit(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

#remove unwanted variables
rm(row.has.na, probes_to_exclude, PCA, ids)
gc()

# remove people without TRAIT info from meth matrix
# mvals_unrelateds_noSNPCH_noNA___TRAIT__ <- mvals
# mvals_unrelateds_noSNPCH_noNA___TRAIT__ <- mvals_unrelateds_noSNPCH_noNA[,-which(is.na(pdata_unrelateds$__TRAIT__))]

# logit(paste0("mvals_unrelateds_noSNPCH_noNA___TRAIT__\t", nrow(mvals_unrelateds_noSNPCH_noNA___TRAIT__), "\t", ncol(mvals_unrelateds_noSNPCH_noNA___TRAIT__)))

# subset phenodata to remove people without TRAIT info
# pdata_unrelateds___TRAIT__ <- pdata_unrelateds[-which(is.na(pdata_unrelateds$__TRAIT__)),]
# pdata_unrelateds___TRAIT__ <- pdata

# fit linear models to identify differentially methylated positions
# The next part of the script involves fitting linear models with 0, 20 or 50 PCs.
# These steps are the same for all traits. To adapt this part of the script for different traits,
# you just need to change the names of the pdata and mvals objects
# (i.e. pdata_unrelateds_BMI and mvals_unrelateds_noSNPCH_noNA_BMI) and change the first parameter
# passed to the function model.matrix to create the "design" variable.
# For example, to change BMI to a new trait, xxx, this line would become
# "design <- model.matrix(~pdata_unrelateds_xxx$xxx...)",
# where the xxx after the $ corresponds the column heading of interest in the pdata object.

pdata$ever_smoke[is.na(pdata$ever_smoke)] <- 5
design_20 <- model.matrix(~pdata$__TRAIT__ + pdata$age + pdata$sex + as.factor(pdata$ever_smoke) + pdata$pack_years + pdata$PC1 + pdata$PC2 + pdata$PC3 + pdata$PC4 + pdata$PC5 + pdata$PC6 + pdata$PC7 + pdata$PC8 + pdata$PC9 + pdata$PC10 + pdata$PC11 + pdata$PC12 + pdata$PC13 + pdata$PC14 + pdata$PC15 + pdata$PC16 + pdata$PC17 + pdata$PC18 + pdata$PC19 + pdata$PC20)

rm_NA_20 <- as.data.frame(cbind(pdata$__TRAIT__, pdata$age, pdata$sex, pdata$pack_years))

uniq_20 <- unique(which(is.na(rm_NA_20), arr.ind=TRUE)[,1])

if(length(uniq_20) == 0) {
	cat(paste0("DELETE\t0\n"))
} else {
	cat(paste0("DELETE\t", uniq_20, "\n"))
	mvals <- mvals[,-uniq_20]
}
logit(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

library(limma)
fit_20 <- lmFit(mvals, design_20)
efit_20 <- eBayes2(fit_20)	
TT_20 <- topTable(efit_20, coef=2, adjust='fdr', number=dim(mvals)[1])
TT_20$ID<-rownames(TT_20)

# Merge the annotation data
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
TT_20 <-merge(genes,TT_20 , by="ID", all.y=TRUE)

TT_20$df <- fit_20$df.residual

write.table(TT_20, file="STRADL_relateds___TRAIT___resid_mvals_20PC_EWAS_toptable.txt", sep="\t", col.names=NA)

log_p_20 <- -log10(TT_20$P.Value)
log_p_20 <- sort(log_p_20)
uni_pvalues <- runif(length(log_p_20))
log_uni_pvalues <- sort(-log10(uni_pvalues))

postscript(file="STRADL_relateds___TRAIT___resid_mvals_20PC_QQ.ps", paper='a4', width=11, height=8, pagecentre=TRUE)
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main="__TRAIT__ relateds 20PC QQ", cex.lab=1.5, cex.main=1.5)
abline(0,1)
dev.off()

jpeg(file="STRADL_relateds___TRAIT___resid_mvals_20PC_QQ.jpeg")
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main="__TRAIT__ relateds 20PC QQ", cex.lab=1.5, cex.main=1.5)
abline(0,1)
dev.off()

jpeg(file="STRADL_relateds___TRAIT___resid_mvals_20PC_QQ_fine.jpeg")
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main="__TRAIT__ relateds 20PC QQ", cex.lab=1.5, cex.main=1.5, col="red", pch='.')
abline(0,1)
dev.off()

jpeg(file="STRADL_relateds___TRAIT___resid_mvals_20PC_QQ_zoom.jpeg")
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main="__TRAIT__ relateds 20PC QQ", cex.lab=1.5, cex.main=1.5, col="red", pch='.', xlim=c(0,2), ylim=c(0,2))
abline(0,1)
dev.off()

IF_20 <- median(TT_20$t^2)/0.4549
cat(paste0("Inflation factor for 20 PC: ", IF_20, "\n"))
