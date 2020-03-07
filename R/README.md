## Building R packages
- [How to build an R package](https://andrewbtran.github.io/NICAR/2018/workflow/docs/01-workflow_intro.html?utm_content=buffer858fd&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)
- run this after updading code: library(roxygen2); roxygenise()

## Installing older versions of an R package
Two examples of how to install older package in R
* Colin Fay https://colinfay.me/docker-r-reproducibility/.
    * install.packages('remotes')
    * remotes::install_version('Seurat', '2.3.4')  # ←- however, I’m finding this to be flaky. Sometimes it finds 2.3.4 and other times it does not find the package.
* From Seurat website: Replace '2.3.0' with your desired version
    * devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))
    * library(Seurat)

## Switching R on Mac RStudio
- See this [documentation](https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-RStudio-desktop)
- Do not use `ln -s`; it is too trick to do properly
- Best thing is to `Run the installer from CRAN for the R version you want to be current`
- [Version 3.5.3](https://cran.r-project.org/bin/macosx/el-capitan/base/)

### tidyverse
* evalutaion
    * https://github.com/cwickham/quotation - !! bang bang
    * https://edwinth.github.io/blog/dplyr-recipes/
    * [Lynda Tidyverse](https://www.lynda.com/R-tutorials/Non-standard-evaluation-programming-tidyverse/586672/649013-4.html?autoplay=true)

## Learning Bioconductor
- [BiocWorkshops](https://bioconductor.github.io/BiocWorkshops/)

## Learning R
- [swirl](https://swirlstats.com/)
