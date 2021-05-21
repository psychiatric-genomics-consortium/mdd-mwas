
# set up install location for required libraries
rlibrary <- file.path('/exports/eddie/scratch', Sys.getenv('USER'), 'Rlibrary', paste(R.version$major, R.version$minor, sep='.'))

dir.create(rlibrary, showWarnings=FALSE)
.libPaths(rlibrary)

# only install library if it isn't already installed
source("https://bioconductor.org/biocLite.R")


if(!library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, lib=rlibrary, logical.return=TRUE)) {biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", ask=FALSE, lib=rlibrary, dependencies=TRUE)}

# probe annotations
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

args <- commandArgs(TRUE)

mvalues_file <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings=FALSE)

mvalues <- readRDS(mvalues_file)

chromosomes <- c(1:22, 'X', 'Y')

for(i in chromosomes) {

  print(paste('Subsetting chromsome', i))
  chr_locations <- subset(Locations, chr %in% paste0('chr', i))
  
  chr_probes <- rownames(chr_locations)
  
  chr_probes_idx <- which(rownames(mvalues) %in% chr_probes)

  mvalues.chrN <- mvalues[chr_probes_idx,]

  print(paste0('Writing object for chromsome ', i, ' (', length(chr_probes_idx), ' probes)'))

  saveRDS(mvalues.chrN, file.path(out_dir, paste0('mvalues.chr', i, '.rds')))

  rm(mvalues.chrN)

  gc()

}
