# set seed so all randomly created distributions are the same
set.seed(1)

# find location of this script
all_args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", all_args[grep(file_flag, all_args)])
script_dir <- dirname(script_path)

.libPaths(file.path(script_dir, 'Rlibrary'))

library(optparse)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

source(file.path(script_dir, 'eBayes2.R'))

parse <- OptionParser()

option_list <- list(
make_option('--pdata', type='character', help="Phenotype CSV with columns 'id' and columns for phenotypes", action='store'),
make_option('--pheno', type='character', help='Name of phenotype column', action='store'),
make_option('--cov', type='character', help='Covariates file (rds)', action='store'),
make_option('--sentrix', type='character', help='Linker between GS IDs and Sentrix IDs (rds)', action='store'),
make_option('--mvals', type='character', help='Mvalue file', action='store'),
make_option('--pca', type='character', help='PCA file', action='store'),
make_option('--probes', type='character', help='Mvalue file', action='store'))



args = commandArgs(trailingOnly=TRUE)

opt <- parse_args(OptionParser(option_list=option_list), args=args)


if(is.null(opt[['pdata']])) {
 stop('Missing phenotype csv argument --pdata')
}
if(is.null(opt[['pheno']])) {
 stop('Missing phenotype argument --pheno')
}

if(is.null(opt[['cov']])) {
 stop('Missing covariates argument --cov')
}

if(is.null(opt[['sentrix']])) {
 stop('Missing Sentrix sample ID argument --sentrix')
}

if(is.null(opt[['mvals']])) {
 stop('Missing Mvalues argument --mvals')
}

if(is.null(opt[['pca']])) {
 stop('Missing PCA argument --pca')
}

if(is.null(opt[['probes']])) {
 stop('Missing probes argument --probes')
}

pheno=opt$pheno
F_PDATA=opt$pdata
F_COV=opt$cov
F_SENTRIX=opt$sentrix
F_MVALS=opt$mvals
F_PCA=opt$pca
F_PROBES=opt$probes


check_file_die <- function(file) {
	if(!file.exists(file)){
		stop(paste0(file, " doesn't exist\n"))
	}
}

rename_file <- function(file) {
	if(file.exists(file)){
		file.rename(file, paste0(file, ".previous"))
	}
}

logging <- function(str) { cat(paste0(seq, "\t", date(), "\t", str, "\n"), file=NOTES, append=TRUE); seq <<- seq+1 }

############################################################

check_file_die(F_PDATA)
check_file_die(F_MVALS)
check_file_die(F_COV)
check_file_die(F_SENTRIX)
check_file_die(F_PCA)
check_file_die(F_PROBES)
NOTES                   <- paste0("Summary-", pheno, ".txt")
seq                     <- 1
rename_file(NOTES)
cat(paste0("n\tTime\tWhat\tRows\tCols\n"), file=NOTES, append=TRUE)

############################################################

pdata <- read.table(F_PDATA, sep=",", header=TRUE)
logging(paste0("pdata\t", nrow(pdata)))

covariates <- readRDS(F_COV)
logging(paste0("cov\t", nrow(cov)))

# load linker file between GS IDs and Sentrix IDs
sentrix <- readRDS(F_SENTRIX)

# looking for outliers
#jpeg("boxplot-__TRAIT__.jpeg")
#length(boxplot(pdata$age, range=3)$out)
#dev.off()

pdata_covariates <- merge(merge(pdata,
                              sentrix, by='id'),
                        covariates, by='id', suffixes=c('.x', ''))

pdata_inputs <- pdata_covariates[,c('id', pheno, 'Sample_Sentrix_ID', names(covariates)[-1])]

pdata_na_rows <- which(rowSums(is.na(pdata_inputs)) > 0)

pdata_complete <- pdata_inputs[-pdata_na_rows,]

logging(paste0("pdata complete\t", nrow(pdata_complete), "\t", nrow(pdata_complete)))




#load norm mvals object
system.time({mvals <- readRDS(F_MVALS)})[3]
logging(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

# Extracting the sample order from the mvals file. We will use this index to maintain this order in the phenotype and PC data frames.
ids <- data.frame(Sample_Sentrix_ID=colnames(mvals), N=c(1:length(colnames(mvals))))
logging(paste0("ids\t", nrow(ids), "\t", ncol(ids)))

#load PCA analysis of residualised M-values 
system.time({PCA <- readRDS(F_PCA)})[3]

#load list of probes affected by SNPs/that are predicted to cross-hybridise
probes_to_exclude <- read.table(F_PROBES, sep="\t", header=FALSE)
logging(paste0("probes\t", nrow(probes_to_exclude)))

pdata_sentrix <- merge(pdata_complete, ids, by="Sample_Sentrix_ID", all.x=TRUE)
logging(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))

##sorting the phenotype file
pdata_sorted <- pdata_sentrix[order(pdata_sentrix$N),]
logging(paste0("pdata\t", nrow(pdata_sorted), "\t", ncol(pdata_sorted)))

mvals_sentrix <- mvals[,which(colnames(mvals) %in% pdata_sorted$Sample_Sentrix_ID)]
logging(paste0("mvals\t", nrow(mvals_sentrix), "\t", ncol(mvals_sentrix)))

logging(paste0("PCA\t", nrow(PCA), "\t", ncol(PCA)))

##Merge PCs to phenotype data
pdata_pcs <- merge(pdata_sorted, PCA, by="Sample_Sentrix_ID", all.x=TRUE)
logging(paste0("pdata\t", nrow(pdata), "\t", ncol(pdata)))
##sorting the phenotype file
pdata_pcs_sorted <- pdata_pcs[order(pdata_pcs$N),]
logging(paste0("pdata\t", nrow(pdata_pcs_sorted), "\t", ncol(pdata_pcs_sorted)))

#Exclude SNP/cross-hybridising probes
mvals_excl <- mvals_sentrix[-which(rownames(mvals_sentrix) %in% probes_to_exclude$V1),]
logging(paste0("mvals\t", nrow(mvals_excl), "\t", ncol(mvals_excl)))

#Remove rows containing NAs
row.has.na <- apply(mvals_excl, 1, function(x){any(is.na(x))})
mvals_qc <- mvals_excl[-which(rownames(mvals_excl) %in% names(which(row.has.na))),]
logging(paste0("mvals\t", nrow(mvals_qc), "\t", ncol(mvals_qc)))


# remove people without TRAIT info from meth matrix
# mvals_unrelateds_noSNPCH_noNA___TRAIT__ <- mvals
# mvals_unrelateds_noSNPCH_noNA___TRAIT__ <- mvals_unrelateds_noSNPCH_noNA[,-which(is.na(pdata_unrelateds$__TRAIT__))]

# logging(paste0("mvals_unrelateds_noSNPCH_noNA___TRAIT__\t", nrow(mvals_unrelateds_noSNPCH_noNA___TRAIT__), "\t", ncol(mvals_unrelateds_noSNPCH_noNA___TRAIT__)))

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

# put formula together for phenotype, covariates, and 20 PCs
model_formula <- as.formula(paste('~', pheno, '+', paste(names(covariates)[-1], collapse=' + '), '+', paste0('PC', 1:20, collapse=' + ')))

design_20 <- model.matrix(model_formula, data=pdata_pcs_sorted)


logging(paste0("mvals\t", nrow(mvals), "\t", ncol(mvals)))

# verify data ordering

if(! all(pdata_pcs_sorted$Sample_Sentrix_ID == names(mvals_qc))) {
        stop('Design matrix ordering does not match Mvalues')
}

fit_20 <- lmFit(mvals_qc, design_20)
efit_20 <- eBayes2(fit_20)	
TT_20 <- topTable(efit_20, coef=2, adjust='fdr', number=dim(mvals)[1])
TT_20$ID<-rownames(TT_20)

# Merge the annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
TT_20 <-merge(genes,TT_20 , by="ID", all.y=TRUE)

TT_20$df <- fit_20$df.residual

write.table(TT_20, file=paste0("STRADL_relateds_", pheno, "_resid_mvals_20PC_EWAS_toptable.txt"), sep="\t", col.names=NA)

log_p_20 <- -log10(TT_20$P.Value)
log_p_20 <- sort(log_p_20)
uni_pvalues <- runif(length(log_p_20))
log_uni_pvalues <- sort(-log10(uni_pvalues))

postscript(file=paste0("STRADL_relateds_", pheno, "_resid_mvals_20PC_QQ.ps"), paper='a4', width=11, height=8, pagecentre=TRUE)
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main=paste(pheno, "relateds 20PC QQ"), cex.lab=1.5, cex.main=1.5)
abline(0,1)
dev.off()

jpeg(file=paste0("STRADL_relateds_", pheno, "_resid_mvals_20PC_QQ.jpeg"))
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main=paste(pheno, "relateds 20PC QQ"), cex.lab=1.5, cex.main=1.5)
abline(0,1)
dev.off()

jpeg(file=paste0("STRADL_relateds_", pheno, "_resid_mvals_20PC_QQ_fine.jpeg"))
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main=paste(pheno, "relateds 20PC QQ"), cex.lab=1.5, cex.main=1.5, col="red", pch='.')
abline(0,1)
dev.off()

jpeg(file=paste0("STRADL_relateds_", pheno, "_resid_mvals_20PC_QQ_zoom.jpeg"))
plot (log_uni_pvalues,log_p_20, xlab="Expected", ylab="Observed", main=paste(pheno,  "relateds 20PC QQ"), cex.lab=1.5, cex.main=1.5, col="red", pch='.', xlim=c(0,2), ylim=c(0,2))
abline(0,1)
dev.off()

IF_20 <- median(TT_20$t^2)/0.4549
logging(paste0("Inflation factor for 20 PC: ", IF_20, "\n"))
