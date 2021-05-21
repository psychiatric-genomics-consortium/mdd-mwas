README file for creating chromosome-separated DNAm data
----------------------------------------------------------

The current package contains the Rscript and files necessary to separate DNAm files by chromosome in order to fit the University of Edinburgh's MWAS pipeline. The package contains:

----------------------------------------------------------

Rscript file: DNAm_separate_files.R
Please load R version 3.6.3

The file requires two arguments, --methdata and --annotfile; the following is an example:
Rscript DNAm_separate_files.R --methdata /path/to/methylation-data.rds --annotfile /path/to/annotation-file.rds

Please note the file creates an output directory named "out" in the current directory, where the script runs. Chromosome-separated files will be stored here, in the following format: "mvalues.chr%.rds" (where % = chr1:22).

----------------------------------------------------------

The package contains two DNAm annotation files in .rds format, depending on which array your DNAm data is in (450K array or EPIC array); please supply either the 450K or EPIC array to argument "--annotfile" when running the Rscript

The two files are derived from the BiocManager IlluminaHumanMethylation packages:

# 450K array

# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# IlluminaHumanMethylation450kmanifest # check manifest file
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# annot.file450 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# EPIC array
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
# IlluminaHumanMethylationEPICmanifest # check manifest file
# library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# annot.file <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

----------------------------------------------------------

For more information about the script, please contact: Dr Miruna C Barbu (mbarbu@ed.ac.uk)
