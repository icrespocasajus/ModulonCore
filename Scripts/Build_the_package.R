# creates description and namespace files
usethis::use_description()
usethis::use_namespace()

# this will setup the folders needed for the data and raw-data
usethis::use_data_raw()
#usethis::use_data(data, overwrite = TRUE)# This line should be included in an RScript to generate the data; a second script (data.R) has to be created, 
# executed and kept in data-raw; a new directory named data with the required files is created

# Create R directory
base::dir.create("R")

# creates Package-level documentation so you can run ?nameofpackage
usethis::use_package_doc()

# created README.Rmd for Github landing page
# an .Rbuildignore file gets created
usethis::use_readme_rmd()
# Creating an R package
#https://sahirbhatnagar.com/rpkg/

# creates license file
usethis::use_mit_license("Sahir Bhatnagar")

# creates news file
usethis::use_news_md()

# setup continuous integration via travis-ci
usethis::use_travis()

# sets up testing infrastructure
usethis::use_testthat()

# this will setup the folders needed for the data and raw-data
usethis::use_data_raw()


# this will generate documentation for your functions
pacman::p_load(sinew)
sinew::makeOxyFile("R/ModulonCore_functions.R")

# See available functions
pacman::p_functions("ModulonCore")

devtools::check()# WARNING! Be aware that using check from RStudio build panel is going to wipe out the doc folder; vignettes may stop working properly

usethis::use_vignette(name = "Introduction")

tools::buildVignettes(dir = ".", tangle=TRUE)
dir.create("inst")
dir.create("inst/doc")
file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)

# Commint and re-build

browseVignettes('ModulonCore')

vignette('Introduction','ModulonCore')

# Dependencies NOT WORKING!
usethis::use_package("limma", type = "Imports")
usethis::use_package("igraph", type = "Imports")
usethis::use_package("SANTA", type = "Imports")


# Check the package
devtools::check()