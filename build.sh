#! /bin/sh

# install required packages to script directory on Eddie

module load igmm/apps/R/3.4.1

module load openmpi/1.10.1

module load igmm/libs/tiff/3.9.3

mkdir Rlibrary

Rscript -e 'install.packages("optparse", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("dplyr", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("truncnorm", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("Rmpi", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE, configure.args=c(Rmpi="--with-mpi=/exports/applications/apps/SL7/openmpi/1.10.1"))'

Rscript -e 'install.packages("doMPI", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("OpenMx", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("RCurl", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("bitops", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("matrixStats", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("foreach", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("iterators", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("HardyWeinberg", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("gee", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("doMC", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("bigmemory", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

#curl -O https://cran.r-project.org/src/contrib/rtiff_1.4.5.tar.gz

#PKG_CPPFLAGS=-I/exports/igmm/software/pkg/el7/libs/tiff/3.9.3/include PKG_LDFLAGS=-L/exports/igmm/software/pkg/el7/libs/tiff/3.9.3/libs R CMD INSTALL --library=Rlibrary rtiff_1.4.5.tar.gz

Rscript -e 'install.packages("qgraph", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("MultiPhen", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'install.packages("semPlot", lib="Rlibrary", repos="https://cran.rstudio.com/", dependencies=TRUE)'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("limma", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocGenerics", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("minfi", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("GenomeInfoDb", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("S4Vectors", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("IRanges", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("GenomicRanges", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("Biobase", ask=F, lib="Rlibrary")'

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", ask=F, lib="Rlibrary")'


