Residualise m/beta-values
================
X Shen
23 August, 2021

-----

This protocol is created for (1) covarying and (2) creating methylation
PCs using m/beta-values. You could choose to use either of the utilities
or both.

The function uses multiple files or a single file.

Contact: xueyi.shen@ed.ac.uk

-----

### Required R packages

  - dplyr

  - plyr

  - FactoMineR

  - optparse

### Function preparation

  - Download the R function from
    [here](https://github.com/psychiatric-genomics-consortium/mdd-mwas/blob/main/methylation%20PCs/resid_pc_mvalue.R)

-----

### File preparation

#### Methylation data

  - **rds** format

  - Row names: CpG probe names (typically start with ‘cg’). Column
    names: Subject ID

  - Data should look like (italic cells are either column/row names):

|              | *Subject1* | *Subject2* | *Subject3* | *…* |
| :----------: | :--------: | :--------: | :--------: | :-: |
| *cg08730728* |  \-0.0012  | \-0.07376  |  0.47836   |  …  |
| *cg03070533* |  0.13523   | \-0.81432  |  0.24235   |  …  |
| *cg21588135* |  \-0.7958  |  0.000257  | \-0.00532  |  …  |
|     *…*      |     …      |     …      |     …      |  …  |

#### Covariates

  - **rds** format

  - First column should be participant IDs, consistent with methylation
    IDs.

  - Suggested covariates: age, sex, batch, assessment centre (please
    include as many as you can)

Data should look like (italic cells are either column/row names):

|   *ID*   | *age* | *sex*  | *batch* |
| :------: | :---: | :----: | :-----: |
| Subject1 |  45   |  Male  |    3    |
| Subject2 |  63   | Female |    1    |
| Subject3 |  60   | Female |    1    |
|    …     |   …   |   …    |    …    |

-----

### Run the script

#### Example scripts (running from bash command line)

For a single file:

``` bash
Rscript resid_pc_mvalue.R \
--meth data/mvalue.chr22.rds \
--covar covs.rds \
--residualise y \
--makePC y \
--out data/
```

For a list of files (large memory required; example file names:
data/mvalues.chr1.rds, data/mvalues.chr2.rds, …):

``` bash
Rscript resid_pc_mvalue.R \
--meth data/mvalue.chr \
--covar covs.rds \
--residualise y \
--makePC y \
--out data/
```

#### Process(es)

  - –residualise: y/n (y = residualise methylation data against
    covariates and save residualised data)

  - –makePC: y/n (y = make methylation PCs)

#### Output(s)

In the given output directory:

  - A list of files for residualised data (end with *\_resid.rds*), if
    ‘resid’ is provided for the ‘process’ option

  - A file called *methPCs.rds*, if ‘pc’ is provided for the ‘process’
    option.
