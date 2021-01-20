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
- Change default ng to NULL in getpared function for automately generate parameter based on data
- update with citation of cc paper
- fix the bug in pmd ms/ms annotation due to the change of enviGCMS package
- fix the order issue from CRAN

# pmd 0.1.9

- CRAN

# pmd 0.1.8

- update kegg/hmdb database
- update getsda to use largest average distance to find pmd freqency cutoff, more robust to large dataset

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
- add merge feature for getcluster and such methods could be used to furthor reduce the GlobalStd peaks
- remove hmdbp data since pmd network analysis could cover this topic
- remove the dependance of group for quantitative paired peaks
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
- add correlationship in getpmd
- add support for quantitative paired peaks list selection for specific reaction
- add support for target pmd and compound analysis for reaction chain

# pmd 0.1.2

- update vignettes
- change default ng value into auto-detection
- add top option to limit sda output
- add suport for GlobalStd based targeted analysis
- add suport to extract specific pmd across different retention time groups
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
