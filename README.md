pmd: Paired Mass Distance Analysis for GC/LC-MS Based Non-Targeted Analysis and Reactomics Analysis
================

[![CRAN status](http://www.r-pkg.org/badges/version/pmd)](https://cran.r-project.org/package=pmd) [![Download counter](http://cranlogs.r-pkg.org/badges/pmd)](https://cran.r-project.org/package=pmd) [![](https://cranlogs.r-pkg.org/badges/grand-total/pmd)](https://cran.r-project.org/package=pmd) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Build status](https://api.travis-ci.org/yufree/pmd.svg?branch=master)](https://travis-ci.org/yufree/pmd)

Paired mass distance (PMD) analysis proposed in Yu, Olkowicz and Pawliszyn (2018) for gas/liquid chromatography–mass spectrometry (GC/LC-MS) based non-targeted analysis. PMD analysis including GlobalStd algorithm and structure/reaction directed analysis. GlobalStd algorithm could found independent peaks in m/z-retention time profiles based on retention time hierarchical cluster analysis and frequency analysis of paired mass distances within retention time groups. Structure directed analysis could be used to find potential relationship among those independent peaks in different retention time groups based on frequency of paired mass distances. Reactomics analysis could also be performed to build PMD network, assign sources and make biomarker reaction discovery. GUIs for PMD analysis is also included as 'shiny' applications.


Installation
------------

You can install package from this GitHub repository:

``` {r}
devtools::install_github("yufree/pmd")
```

Or find a stable version from CRAN:

``` {r}
install.packages('pmd')
```

Usage
-----

[PMD Analysis Tutorial](https://yufree.github.io/pmd/articles/globalstd.html)

[Reactomics Analysis Tutorial](https://yufree.github.io/pmd/articles/reactomics.html)

To perform GlobalStd algorithem, use the following code:

``` {r}
library(pmd)
data("spmeinvivo")
pmd <- getpaired(spmeinvivo, rtcutoff = 10, ng = 10)
std <- getstd(pmd)
```
To perform structure/reaction directed analysis, use the following code:

``` {r}
sda <- getsda(std)
```

To perform GlobalStd algorithem along with structure/reaction directed analysis, use the following code:

``` {r}
result <- globalstd(spmeinvivo)
```

To use the shiny application within the package, use the following code:

```{r}
runPMD()
```

To check the pmd reaction database:

```{r}
# all reaction
data("omics")
View(omics)
# kegg reaction
data("keggrall")
View(keggrall)
# literature reaction for mass spectrometry
data("sda")
View(sda)
```

To check the HMDB pmd database:

```{r}
data("hmdb")
View(hmdb)
```
