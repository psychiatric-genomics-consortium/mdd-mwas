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
make_option('--out', type='character', help='Output prefix', action='store'),
make_option('--cov', type='character', help='Covariates file (rds)', action='store'),
make_option('--sentrix', type='character', help='Linker between GS IDs and Sentrix IDs (rds)', action='store'),
make_option('--mvals', type='character', help='Mvalue file', action='store'),
make_option('--pca', type='character', help='PCA file', action='store'),
make_option('--prune', action='store_true', default=FALSE,  help='Prune probes with missing samples'),
make_option('--exclude', type='character', help='Probes to exclude', action='store'))



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

chr <- as.integer(Sys.getenv('SGE_TASK_ID'))

pheno=opt$pheno
F_PDATA=opt$pdata
F_COV=opt$cov
F_SENTRIX=opt$sentrix
F_PCA=opt$pca
F_PROBES=opt$exclude

out <- opt$out


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

log_file <- paste0(out, ".log")
logging <- function(str) { cat(paste0(paste0(str, collapse=''), '\n'), file=log_file, append=TRUE) }
 
logging('MDD MWAS')

logging(c("Started: ", date()))

logging(c('Phenotype data: ', F_PDATA))
logging(c('Phenotype: ', pheno))
logging(c('Covariates data: ', F_COV))
logging(c('Sentrix ID linker: ', F_SENTRIX))
logging(c('Mvalues pattern: ', opt$mvals))
logging(c('Principal components: ', F_PCA))
logging(c('Probes to exclude: ', F_PROBES))
logging('')
logging(c('Job ID: ', Sys.getenv('JOB_ID')))
logging(c('Job standard error file: ', Sys.getenv('SGE_STDERR_PATH')))
logging(c('Job standard out file: ', Sys.getenv('SGE_STDOUT_PATH')))
logging(c('MWAS script file: ', script_path))
logging('')

############################################################

check_file_die(F_PDATA)
check_file_die(F_COV)
check_file_die(F_SENTRIX)
check_file_die(F_PCA)
if(!is.null(F_PROBES)) {
  check_file_die(F_PROBES)
}
############################################################



pdata <- read.table(F_PDATA, sep=",", header=TRUE)

logging(c('Input sample: ', nrow(pdata)))

covariates <- readRDS(F_COV)


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

pdata_complete_rows <- which(rowSums(is.na(pdata_inputs)) == 0)

pdata_complete <- pdata_inputs[pdata_complete_rows,]

logging(c('Complete phenotype and covariates: ', nrow(pdata_complete)))

#load PCA analysis of residualised M-values 
system.time({PCA <- readRDS(F_PCA)})[3]




#load list of probes affected by SNPs/that are predicted to cross-hybridise
if(!is.null(F_PROBES)) {
  probes_to_exclude <- read.table(F_PROBES, sep="\t", header=FALSE)
} else {
  probes_to_exclude <- data.frame(V1=NULL)
}

logging(c('Probes excluded: ', nrow(probes_to_exclude)))


mwas <- function(chr) {


  cat(paste('MWAS for chromosome', chr, '\n'))
  # get path for this chromosome's mvalues 
  F_MVALS=gsub('%', chr, opt$mvals)

  #load norm mvals object
  system.time({mvals <- readRDS(F_MVALS)})[3]
  
  # Extracting the sample order from the mvals file. We will use this index to maintain this order in the phenotype and PC data frames.
  ids <- data.frame(Sample_Sentrix_ID=colnames(mvals), N=c(1:length(colnames(mvals))))
  
  pdata_sentrix <- merge(pdata_complete, ids, by="Sample_Sentrix_ID", all.x=TRUE)
  
  
  ##sorting the phenotype file
  pdata_sorted <- pdata_sentrix[order(pdata_sentrix$N),]
  
  
  ##Merge PCs to phenotype data
  pdata_pcs <- merge(pdata_sorted, PCA, by="Sample_Sentrix_ID", all.x=TRUE)
  
  ##sorting the phenotype file
  pdata_pcs_sorted <- pdata_pcs[order(pdata_pcs$N),]
  
  mvals_sentrix <- mvals[,which(colnames(mvals) %in% pdata_sorted$Sample_Sentrix_ID)]

  rm(mvals)
  gc()

  #Exclude SNP/cross-hybridising probes
  mvals_excl <- mvals_sentrix[-which(rownames(mvals_sentrix) %in% probes_to_exclude$V1),]

  rm(mvals_sentrix)
  gc()

  #Remove rows containing NAs
  if(opt$prune) {
    row.has.na <- apply(mvals_excl, 1, function(x){any(is.na(x))})
    mvals_qc <- mvals_excl[-which(rownames(mvals_excl) %in% names(which(row.has.na))),]
  } else {
    mvals_qc <- mvals_excl
  }

  rm(mvals_excl)
  gc()
  
  # put formula together for phenotype, covariates, and 20 PCs
  model_formula <- as.formula(paste('~', pheno, '+', paste(names(covariates)[-1], collapse=' + '), '+', paste0('PC', 1:20, collapse=' + ')))
  
  design_20 <- model.matrix(model_formula, data=pdata_pcs_sorted)

  if(chr == 1) {
          logging(c('Model: ', model_formula))
          logging(c('Design matrix: ', paste('~', paste(dimnames(design_20)[[2]], collapse=' + '))))
          logging(c('MWAS sample size: ', nrow(design_20)))
  }

  
  # verify data ordering
  
  if(! all(pdata_pcs_sorted$Sample_Sentrix_ID == names(mvals_qc))) {
          stop('Design matrix ordering does not match Mvalues')
  }
  
  fit_20 <- lmFit(mvals_qc, design_20)

  fit_20$N <- rowSums(!is.na(mvals_qc))

  rm(design_20, mvals_qc)
  gc()

  return(fit_20)

}

 

cat(paste(date(), 'Starting MWAS', '\n'))
# load model fits from each chromosome
fits_chr <- lapply(1:22, mwas)

# merge into a single MArrayLM
cat(paste(date(), 'Merging MWAS', '\n'))
names(fits_chr[[1]])

fits <- 
new("MArrayLM",
list(
coefficients=do.call(rbind, lapply(fits_chr, function(m) m[['coefficients']])),
rank=fits_chr[[1]]$rank,
qr=fits_chr[[1]]$qr,
df.residual=do.call(c, lapply(fits_chr, function(m) m[['df.residual']])),
sigma=do.call(c, lapply(fits_chr, function(m) m[['sigma']])),
cov.coefficients=fits_chr[[1]]$cov.coefficients,
stdev.unscaled=do.call(rbind, lapply(fits_chr, function(m) m[['stdev.unscaled']])),
pivot=fits_chr[[1]]$pivot,
Amean=do.call(c, lapply(fits_chr, function(m) m[['Amean']])),
method=fits_chr[[1]]$method,
design=fits_chr[[1]]$design,
N=do.call(c, lapply(fits_chr, function(m) m[['N']]))))


cat(paste(date(), 'eBayes analysis', '\n'))
efit <- eBayes2(fits)	

TT <- topTable(efit, coef=2, adjust='fdr', number=length(fits$Amean))

TT$ID<-rownames(TT)

logging(c('MWAS probes: ', nrow(TT)))

# Merge the annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
TT_anno <-merge(genes, TT, by="ID", all.y=TRUE)

TT_anno$df <- fits$df.residual

# linear models

fits_n <- data.frame(ID=names(fits$N),  N=fits$N)

pheno_betas <- data.frame(ID=row.names(fits$coefficients), beta=fits$coefficients[,2])

TT_linear <- TT_anno;
TT_linear$beta <- fits$coefficients[as.character(TT_linear$ID),pheno];
TT_linear$se <- TT_linear$beta / TT_linear$t
TT_linear$N <- fits$N[as.character(TT_linear$ID)];

efit_betas <- plyr::adply(efit$coefficients, 1, function(x) data.frame(param=names(x), beta=x))

names(efit_betas)[1] <- 'ID'

efit_t <- plyr::adply(efit$t, 1, function(x) data.frame(param=names(x), t=x))

names(efit_t)[1] <- 'ID'

efit_p <- plyr::adply(efit$p.value, 1, function(x) data.frame(param=names(x), p=x))

names(efit_p)[1] <- 'ID'

efit_linear <- efit_betas
efit_linear$t <- efit_t$t
efit_linear$p <- efit_p$p
efit_linear$se <- efit_linear$beta / efit_linear$t


tt_output_file <- paste0(out, ".toptable.txt")
lm_output_file <- paste0(out, ".linear.txt")

logging(c('Top-ranked genes: ', tt_output_file))
logging(c('Linear model output: ', lm_output_file))
cat(paste(date(), 'Writing results', '\n'))
write.table(TT_linear, file=tt_output_file, sep="\t", row.names=FALSE, quote=F)
write.table(efit_linear[,c('ID', 'param', 'beta', 'se', 't', 'p')], file=lm_output_file, sep="\t", row.names=FALSE, quote=F)

logging(c('Finished: ', date()))
