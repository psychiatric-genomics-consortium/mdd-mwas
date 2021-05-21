# EWAS using generalized estimating equations


# find location of this script
all_args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", all_args[grep(file_flag, all_args)])
script_dir <- dirname(script_path)

rlib <- file.path(script_dir, 'Rlibrary')

library(optparse, lib.loc=rlib)
library(gee, lib.loc=rlib)
library(bigmemory)
library(doMC)
library(dplyr, lib.loc=rlib)


option_list <- list(
make_option('--pheno', type='character', help="Phenotype CSV with columns 'id' and columns for phenotypes", action='store'),
make_option('--pheno-name', type='character', help='Name of phenotype column', action='store', dest='pheno_name'),
make_option('--out', type='character', help='Output prefix', action='store'),
make_option('--cov', type='character', help='Covariates file (rds)', action='store'),
make_option('--smoking-name', type='character', help='Name of smoking covariate', action='store', dest='smoking_name'),
make_option('--sentrix', type='character', help='Linker between GS IDs and Sentrix IDs (rds)', action='store'),
make_option('--fam', type='character', help='PLINK fam file', action='store'),
make_option('--betas', type='character', help='Normed betas file pattern', action='store'),
make_option('--exclude', type='character', help='Probes to exclude', action='store'))


args = commandArgs(trailingOnly=TRUE)

opt <- parse_args(OptionParser(option_list=option_list), args=args)

if(is.null(opt[['pheno']])) {
 stop('Missing phenotype csv argument --pheno')
}
if(is.null(opt[['pheno_name']])) {
 stop('Missing phenotype name argument --pheno-name')
}

if(is.null(opt[['cov']])) {
 stop('Missing covariates argument --cov')
}

if(is.null(opt[['smoking_name']])) {
 stop('Missing smoking covariate name argument --smoking-name')
}
if(is.null(opt[['sentrix']])) {
 stop('Missing Sentrix sample ID argument --sentrix')
}

if(is.null(opt[['betas']])) {
 stop('Missing Betas argument --betas')
}

if(is.null(opt[['fam']])) {
 stop('Missing FAM file argument --fam')
}

logging('STRADL EWAS using GEE')

logging(c("Started: ", date()))
logging(c('Phenotype data: ', opt$pheno))
logging(c('Phenotype: ', opt$pheno_name))
logging(c('Covariates data: ', opt$cov))
logging(c('Smoking covarite: ', opt$smoking_name))
logging(c('Sentrix ID linker: ', opt$sentrix))
logging(c('Fam file: ', opt$fam))
logging(c('Betas pattern: ', opt$betas))
logging(c('Probes to exclude: ', opt$exclude))
logging('')
logging(c('Job ID: ', Sys.getenv('JOB_ID')))
logging(c('Job standard error file: ', Sys.getenv('SGE_STDERR_PATH')))
logging(c('Job standard out file: ', Sys.getenv('SGE_STDOUT_PATH')))
logging(c('EWAS script file: ', script_path))
logging('')

pheno <- read.csv(opt$pheno)

pheno_name <- opt$pheno_name

covar <- readRDS(opt$cov)

sentrix <- readRDS(opt$sentrix)

fam <- read.table(opt$fam, col.names=c('FID', 'IID', 'father', 'mother', 'sex_code', 'phenotype'))

if(!is.null(opt$exclude)) {
  probes_exclude <- read.table(opt$exclude, stringsAsFactors=F)$V1
} else {
  probes_exclude <- c()
}

sample_data <- pheno %>%
inner_join(covar, by='id') %>%
inner_join(fam %>% select(FID, IID), by=c('id'='IID')) %>%
inner_join(sentrix, by='id')

geewas <- function(chr) {

        betas <- as.data.frame(t(readRDS(gsub('%', chr, opt$betas))))

        probes <- colnames(betas)

        analyze_probes <- probes[which(! probes %in% probes_exclude)]

        betas$Sample_Sentrix_ID <- rownames(betas)

        samples_betas <- sample_data %>%
        inner_join(betas, by='Sample_Sentrix_ID') %>%
        arrange(FID)

        rm(betas)
        gc()

        beta_big <- as.big.matrix(samples_betas[,analyze_probes])

        samples <- samples_betas %>%
        select(-one_of(probes))

        k <- 4
        registerDoMC(k)
        
        m1_predictor_formula <- as.formula(paste('. ~', paste(c(pheno_name, names(covar)[which(!names(covar) %in% 'id')]), collapse='+')))

        
        models1 <- 
        foreach(prob=analyze_probes) %dopar% {
                model1_formula <- update(m1_predictor_formula, paste(prob, '~.'))
                samples_cpgi <- samples
                samples_cpgi[,prob] <- beta_big[,prob]

                model1 <- gee(model1_formula, id=FID, family='gaussian', corstr='exchangeable', data=samples_betas, maxiter=100, na.action=na.omit)
                return(model1)
        }

        return(models1)

}
