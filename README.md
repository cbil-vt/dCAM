
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dCAM

<!-- badges: start -->

<!-- badges: end -->

The goal of dCAM is to …

## Installation

You can install the development version of dCAM like so:

``` r
install.packages("devtools")
pak::pkg_install("cbil-vt/dCAM")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dCAM)
#> my Java class path: [1] "C:\\Users\\freem\\AppData\\Local\\R\\cache\\R\\renv\\library\\dCAMPkg-6f7bcdf8\\windows\\R-4.6\\x86_64-w64-mingw32\\rJava\\java"                 
#> [2] "C:\\Users\\freem\\AppData\\Local\\R\\cache\\R\\renv\\library\\dCAMPkg-6f7bcdf8\\windows\\R-4.6\\x86_64-w64-mingw32\\dCAM\\java"                  
#> [3] "C:\\Users\\freem\\AppData\\Local\\R\\cache\\R\\renv\\library\\dCAMPkg-6f7bcdf8\\windows\\R-4.6\\x86_64-w64-mingw32\\dCAM\\java\\CornerDetect.jar"
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
