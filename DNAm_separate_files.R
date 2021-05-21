## Separate DNA methylation data into chromosome-separate files
# Miruna Carmen Barbu - 23/04/2021
# For more information on the script and files included, please read the 'README DNAm separate files.txt' file

# Information on file formats ----------------------------------------------------------
# Included with script: EPIC and 450K array annotation files
# DNAm data: DNAm participant ID (columns) and CpG (rows) format, in .rds file

# R packages ----------------------------------------------------------
library(dplyr) # Data management
library(optparse) # Data management

# Parse arguments ----------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
parse <- OptionParser()
option_list <- list(
  make_option('--methdata', type='character', help="Location for DNAm data file", action='store'),
  make_option('--annotfile', type='character', help="Annotation file for 450K array or EPIC array", action='store')
)
args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)

methyl.data = opt$methdata
annotat.file = opt$annotfile
# Create an output directory in the current file location to store DNAm data
dir.create("./out")

# Load data into R ---------------------------------------------------------------
# DNA methylation data
meth.data <- readRDS(methyl.data)
# Annotation file
annot.file <- readRDS(annotat.file)

# Pre-process the annotation file ---------------------------------------------------------------
# Separate the file "annot.file" into separate chromosome files
list.df <- as.data.frame(annot.file) # turn to dataframe
list.chrom <- split(list.df, list.df$chr) # split the dataset into a list of datasets based on ann$chr

# Remove unneccessary data ---------------------------------------------------------------
rm(annot.file)
rm(list.df)

# Subset ---------------------------------------------------------------
# Take DNAm file, and subset CpGs based on those corresponding to different chromosomes in list_df
n <- 1:24 # Change according to number of elements in list
meth.sep = lapply(n, function(x) subset(meth.data, rownames(meth.data) %in% list.chrom[[x]]$Name)) # Subset
names(meth.sep) = names(list.chrom)

# Save ---------------------------------------------------------------
for (df in names(meth.sep)) saveRDS(meth.sep[[df]], file=paste0("./out/mvalues.", df, '.rds'))
