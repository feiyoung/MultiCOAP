# MultiCOAP
High-dimensional covariate-augmented multi-study nonlinear factor model

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/MultiCOAP)](https://cran.r-project.org/package=MultiCOAP)
[![](https://cranlogs.r-pkg.org/badges/MultiCOAP?color=orange)](https://cran.r-project.org/package=MultiCOAP)
[![](https://cranlogs.r-pkg.org/badges/grand-total/MultiCOAP?color=orange)](https://cran.r-project.org/package=MultiCOAP)
<!-- badges: end -->


We introduce a novel modeling strategy for analyzing single-cell multi-omics data sourced from various studies. Our approach is driven by the surge in single-cell multi-omics research in biological and medical fields, exemplified by the RNA-protein simultaneous sequencing case-control study on human peripheral blood mononuclear cells (PBMC). To capture the nonlinear relationships between diverse molecular types, including genes and proteins, and to identify both study-shared and study-specific features while accounting for the count nature of the data, we propose a multi-study nonlinear factor model incorporating study-shared and study-specific factors, along with augmented omics variables as covariates. 



Check out the  [Package Website](https://feiyoung.github.io/MultiCOAP/index.html) for a more complete description of the methods and analyses. 

# Installation
"MultiCOAP" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
## Method 1：
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/MultiCOAP")

## Method 2: install from CRAN
install.packages("MultiCOAP")

```



## Tutorials
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Simulated data](https://feiyoung.github.io/MultiCOAP/articles/MultiCOAPsimu.html)

## Simulated codes
For the codes in simulation study, check the `simuCodes` directory of the repo.


## News

MultiCOAP version 1.1 released! (2024-02-25) 


