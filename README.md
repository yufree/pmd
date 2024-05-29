pmd: Paired Mass Distance Analysis for GC/LC-MS Based Non-Targeted Analysis and Reactomics Analysis
================

[![CRAN status](http://www.r-pkg.org/badges/version/pmd)](https://cran.r-project.org/package=pmd) [![Download counter](http://cranlogs.r-pkg.org/badges/pmd)](https://cran.r-project.org/package=pmd) [![](https://cranlogs.r-pkg.org/badges/grand-total/pmd)](https://cran.r-project.org/package=pmd) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/) 

Paired mass distance (PMD) analysis proposed in Yu, Olkowicz and Pawliszyn (2018) and PMD based reactomics proposed in Yu and Petrick (2020) for gas/liquid chromatography–mass spectrometry (GC/LC-MS) based non-targeted analysis. PMD analysis including GlobalStd algorithm and structure/reaction directed analysis. GlobalStd algorithm could found independent peaks in m/z-retention time profiles based on retention time hierarchical cluster analysis and frequency analysis of paired mass distances within retention time groups. Structure directed analysis could be used to find potential relationship among those independent peaks in different retention time groups based on frequency of paired mass distances. Reactomics analysis could also be performed to build PMD network, assign sources and make biomarker reaction discovery. GUIs for PMD analysis is also included as 'shiny' applications.


Installation
------------

You can install package from this GitHub repository:

``` {r}
remotes::install_github("yufree/pmd")
```

Or find a stable version from CRAN:

``` {r}
install.packages('pmd')
```

Usage
-----

- [PMD Analysis Tutorial](https://yufree.github.io/pmd/articles/globalstd.html)

- [Reactomics Analysis Tutorial](https://yufree.github.io/pmd/articles/reactomics.html)

- For MS only data such as FT-ICR or MS imaging data, check the section of "Reactomics analysis for MS only data" in [[Reactomics Analysis Tutorial]] or this blog [post](https://yufree.cn/en/2024/05/29/reactomics-analysis-for-ms-only-data/).

- [PMD analysis paper](https://www.sciencedirect.com/science/article/abs/pii/S0003267018313047). This paper proposed structure/reaction directed analysis with PMD.

- [Reactomics paper](https://www.nature.com/articles/s42004-020-00403-z). This paper contained the concepts of PMD based reactomics, applications and data mining of reaction database and compounds database.

- [Slides](http://yufree.github.io/presentation/reactomics/pres-asms.html). This is the slides for ASMS 2020 Reboot and here is the [video](https://youtu.be/-mT3HcVygHE) of presentation. Press "P" and you will see the notes for each slide with details. Another full version of reactomics presentation for one hour presentation could be found [here](http://yufree.github.io/presentation/reactomics/pres). I will not update the conference presentation while I will add new contents for the full version of reactomics presentation whenever I have new results.


To perform GlobalStd algorithm, use the following code:

``` {r}
library(pmd)
data("spmeinvivo")
pmd <- getpaired(spmeinvivo)
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

To cite related papers:

```{r}
citation('pmd')
```
