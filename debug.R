#debug.R
#remotes::install_github("extendr/rextendr")

rm(list=ls())

library(rextendr)
library(rustytools)

roxygen2::roxygenise()


rextendr::clean()
rextendr::document()

get_consensus("NNNNNNAATGNNNNGGGNNN")
get_consensus("NNNnnnAATGaaaaGGGNNN")
get_consensus("AATGNNNGGG")


# Run once to configure package to use pkgdown
#usethis::use_pkgdown()
usethis::use_pkgdown_github_pages()
# Run to build the website

