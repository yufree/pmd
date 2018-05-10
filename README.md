pmd: Paired Mass Distance(PMD) analysis for GC/LC-MS based non-targeted analysis
================

[![CRAN status](http://www.r-pkg.org/badges/version/pmd)](https://cran.r-project.org/package=pmd) [![Download counter](http://cranlogs.r-pkg.org/badges/pmd)](https://cran.r-project.org/package=pmd) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Build status](https://api.travis-ci.org/yufree/pmd.svg?branch=master)](https://travis-ci.org/yufree/pmd)

`pmd` provides functions for Paired Mass Distance(PMD) analysis including GlobalStd algorithm and structure/reaction directed analysis. GlobalStd algorithm could found independant peaks in m/z-retention time profiles based on retention time hierarchical cluster analysis and frequency analysis of paired mass distances within retention time groups. Structure directed analysis could be used to find potential relationship among those independant peaks in different retention time groups based on freqency of paired mass Distances. A GUI for PMD analysis is also included as a shiny application.

Installation
------------

You can install package from this GitHub repository:

``` {r}
devtools::install_github("yufree/pmd")
```

Usage
-----

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
