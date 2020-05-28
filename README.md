## About this repository
This is set of files to support "Probability of Fit Failure with Reuse of N95 Respirators"

The file n95_report.R does most of the heavy lifting and caches the longer calculation results. 
It assumes it is run in the directory contain the data file "Final Pilot Data.csv"
Required R packages are "magrittr", "dplyr", "interval", "cgam", "icenReg", "km.ci", "flexsurv" , "caret", "tableone", "rmarkdown", "knitr", "bookdown", "PropCIs"
I have run it in R 3.5.1 and 3.6.3.

The file n95_report.Rmd is an rmarkdown file that generates the web appendix supplementary results.
It looks for the cached results from n95_report.R and runs it if not present.
It is designed to be invoked with R -e 'rmarkdown::render("n95_report.Rmd") '
It depends on a markdown and LaTeX install with lmodern.sty, pandoc-citeproc

