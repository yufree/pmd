# pmd 0.2.7

- CRAN

# pmd 0.2.6

- Add `getpmddf` to show pmd details with given m/z only data and m/z group information(optional, segmentation m/z group for spatial reactomics analysis)
- Update reactomics vignette with section "Reactomics analysis for MS only data" and showcase the quantitative analysis for certain PMD
- add support for mass only data for `getchain`

# pmd 0.2.5

- Add KEGG reaction class and enzyme number to the keggrall database
- Fix url of demo data

# pmd 0.2.4

- add support for mass only data for quantitative reactomics analysis
- add quantitative methods for dynamic pmds
- add support for multiple pmds in getpmd function
- update reactomics vignette to add more details for quantitative analysis of PMD
- update getchain to handle large data
- add mass defect filter for `getrda` and `getpaired` to retain reaction related PMDs
- add parameter for `getrda` for pmd sets and mass defect table
- change `getchain` corcutoff to 0.6

# pmd 0.2.3

- add MaConDa database
- fix NULL default value issue in shiny apps

# pmd 0.2.2

- spell check
- goodpractice package check

# pmd 0.2.1

- CRAN

# pmd 0.2.0

- add option to skip sda in GlobalStd algorithm and set default to F
- organize the R files
- add vignette for reactomics analysis
- add correlation directed analysis function
- modified getcorcluster function to find independent peaks
- add vignette section for reduced independent peaks selection in GlobalStd algorithm
- fix the issue for getchain with multiple masses
- fix the correlation issue in pos/neg linkage function
- Output within RT clusters high frequencies PMD(s) as message for user to check
- Change default ng to NULL in getpared function for automated generate parameter based on data
- update with citation of cc paper
- fix the bug in pmd ms/ms annotation due to the change of enviGCMS package
- fix the order issue from CRAN

# pmd 0.1.9

- CRAN

# pmd 0.1.8

- update kegg/hmdb database
- update getsda to use largest average distance to find pmd frequency cutoff, more robust to large dataset

# pmd 0.1.7

- add function for pmd ms/ms annotation
- add function to read in msp file as database
- detach rcdk package
- add function to link pos/neg by pmd

# pmd 0.1.6

- rewrite getchain to speed up
- add shiny application pmdnet to perform PMD network analysis
- add support for formula in getchain to find compounds ions
- remove frequency cutoff in getsda and use PMD network clusters analysis to determine the cutoff
- add merge feature for getcluster and such methods could be used to further reduce the GlobalStd peaks
- remove hmdbp data since pmd network analysis could cover this topic
- remove the dependence of group for quantitative paired peaks
- improve shiny application for sda analysis

# pmd 0.1.5

- CRAN

# pmd 0.1.4

- isotope selection improved to get rid of 1&2 issue 
- fix top issue in getsda
- fix peak index issue in getcluster
- add corcutoff for getpmd

# pmd 0.1.3

- speed up GlobalStd by mapply
- add message for getrda
- add hmdb pmd analysis results as dataset
- add digits for mass accuracy
- fix the ms1 larger issue in getpmd
- add correlations in getpmd
- add support for quantitative paired peaks list selection for specific reaction
- add support for target pmd and compound analysis for reaction chain

# pmd 0.1.2

- update vignettes
- change default ng value into auto-detection
- add top option to limit sda output
- add support for GlobalStd based targeted analysis
- add support to extract specific pmd across different retention time groups
- add PCA similarity factor function from EvolQG package
- add support to export std peaks based on correlation within retention time group
- add support to export index for peaks with highest intensity in peaks cluster
- add support to use intensity data to refine GlobalStd results
- add support to generate sda analysis for mass list only #5
- remove multi chargers with a strict rule #4
- add parameter selection part in vignette #3
- add support for peaks cluster output #2
- add support for formula generation in enviGCMS package #1

# pmd 0.1.1

- CRAN with ACA [publication](https://doi.org/10.1016/j.aca.2018.10.062)

# pmd 0.1.0

- new package for paired mass distance analysis
