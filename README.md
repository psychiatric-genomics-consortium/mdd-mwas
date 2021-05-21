# MWAS
MWAS scripts

Scripts to perform MWAS on cohort data. 

**Main scripts that are needed for conducting MDD MWAS**

mdd_mwas

mdd_mwas.R

mdd_mwas.sh

ebayes2.R

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**INSTRUCTIONS**

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**DEFAULT FILES AND ARGUMENTS**

The following arguments in the pipeline are default. You will need to change the script to reflect the paths to these files. 

FILE: **mdd_mwas**
Lines: 8-13; 46-51; 60-65

**Paths to change:** all paths to change are noted as "/path/to/[respective file]" in **mdd_mwas**. 

1. Covariates: RDS FILE; there is a default covariate file, although you can override this with an argument in the pipeline while you run MWAS (e.g. --cov /path/to/covariate_file). It is good to have a base covariate file with e.g. cohort IDs, age, and sex, in the pipeline script.
     
     **FORMAT:** the covariate file should have "id" for cohort ID and covariates to include in MWAS (any name)
     
2. Sentrix ID linker: RDS FILE;  this is a file with 2 columns, one with cohort IDs and one with DNA methylation IDs. 
    
    **FORMAT:** the ID linker file should have two columns with following column names: "id" (cohort IDs) and "Sample_Sentrix_ID" (DNAm IDs). We recommend renaming your respective cohort and DNAm IDs with these to automate the pipeline
     
3. M-value chromosome files: RDS FILES;  DNAm data should be chromosome-separated, if possible. The path here should be: /path/to/mvalue/files/mvalues.chr%.rds ; the m-value files can be renamed, although they should have the chromosome number somewhere in the file name.
     
     **FORMAT:** all m-value DNAm files should have CpGs as rows and participants as columns (participant IDs should be DNAm IDs)
     
4. PCA file: RDS FILE; this should be a path to a file containing 20 principal components calculated for the dataset. *If PCs have already been accounted for during DNAm pre-processing, please let us know and we will provide a modified pipeline that does not take PCs into account*.
     
     **FORMAT:** PCA file should contain a column with DNAm IDs and first 20 principal components, named "PC 1:20"
     
5. Prune pipeline argument: the default for pruning probes with missing samples; this can be changed to "--no-prune", which will retain probes with missing samples.

6. Probes to exclude: TEXT FILE;  this should point to a file containing a list of CpG sites (no header) that are cross-hybridising or polymorphic. This will remove these CpGs prior to running MWAS. If you wish to retain these CpGs, supply an empty file here.
      
      **FORMAT:** textfile containing a list of CpGs to remove (no header)

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**INPUT FILES**

The following arguments in the pipeline are free text input. You will need input the paths and file names for these based on the phenotype and covariate files you are analysing.

FILE: **mdd_mwas**
Lines: 6:8, 14; 44:46, 52; 58:60, 66 (DO NOT change anything in the main script; this is just for information purposes)

1. Phenotype file: .csv format
       
      **Format:** .csv file containing cohort ID ("id"), as well as a list of phenotypes to analyse (the script will take one phenotype at a time)
2. "pheno" argument: the main phenotype to analyse in "Phenotype file"
3. Covariate file: pleasre refer to point 1 above in "Default files and arguments"
       
      **FORMAT:** the covariate file should have "id" for cohort ID and covariates to include in MWAS (any name)
4. "Out": specifies the location where MWAS results and analysis log file will be output to

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**PIPELINE FORMAT**

1. First, add your local path to your .bash_profile

PATH=$PATH: /path/to/local/bin

2. Update your path with the shell command:

source ~/.bash_profile

3. Load R; PLEASE NOTE: the pipeline runs with R version 3.3.2, therefore all .rds files should be saved in this version

The pipeline can be started directly from the command line, as follows:

mdd_mwas --pdata /path/to/phenotype-file.csv --pheno [phenotype name in phenotype-file.csv] --cov /path/to/covariate-file.rds [optional, if you do not want to prune probes with missing samples] --no-prune --out /path/to/output/directory/phenotype-mwas [the output file does not require any extension, as it will output 3 files]

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**Other information**

TEMPLATE.R has the token \____TRAIT\___ which should be replaced by the trait to be analysed

eBayes2.R is an updated ebayes function by Dr Rosie Walker to get round integer overflow problem (essentially adds as.numeric() )

split_chr_betas, split_chr_mvalues, and DNAm_separate_files are different ways to separate DNA methylation files based on chromosomes (the pipeline runs the MWAS chromosome-by-chromosome).


