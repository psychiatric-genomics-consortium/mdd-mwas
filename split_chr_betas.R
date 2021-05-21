
# find location of this script
all_args <- commandArgs(trailingOnly = FALSE)
file_flag <- "--file="
script_path <- sub(file_flag, "", all_args[grep(file_flag, all_args)])
script_dir <- dirname(script_path)

rlib = file.path(script_dir, 'Rlibrary')

# probe annotations
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, lib.loc=rlib)

args <- commandArgs(TRUE)

betas_file <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings=FALSE)

load(betas_file)

chromosomes <- c(1:22, 'X', 'Y')

for(i in chromosomes) {

  print(paste('Subsetting chromsome', i))
  chr_locations <- subset(Locations, chr %in% paste0('chr', i))
  
  chr_probes <- rownames(chr_locations)
  
  chr_probes_idx <- which(rownames(x) %in% chr_probes)

  betas.chrN <- x[chr_probes_idx,]

  print(paste0('Writing object for chromsome ', i, ' (', length(chr_probes_idx), ' probes)'))

  saveRDS(betas.chrN, file.path(out_dir, paste0('betas.chr', i, '.rds')))

  rm(betas.chrN)

  gc()

}

saveRDS(x, file.path(out_dir, 'norm_betas_NAs_replaced.rds'))
