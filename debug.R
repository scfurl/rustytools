#debug.R
#remotes::install_github("extendr/rextendr")

rm(list=ls())

library(rextendr)
library(rustytools)

roxygen2::roxygenise()


rextendr::clean()
rextendr::document()

getconsensus("NNNNNNAATGNNNNGGGNNN", 0)
getconsensus("NNNnnnAATGaaaaGGGNNN", 0)
getconsensus("AATGNNNGGG", 0)


# Run once to configure package to use pkgdown
#usethis::use_pkgdown()
usethis::use_pkgdown_github_pages()
# Run to build the website

