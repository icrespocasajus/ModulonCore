
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ModulonCore

<!-- badges: start -->
<!-- badges: end -->

The goal of ModulonCore is to identify the modulon regulatory cores or
connected components with statistically supported over-connectivity

## Installation

You can install the development version of ModulonCore from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("icrespocasajus/ModulonCore")
```

## Example

This is a basic example which shows you how to find and explore the
regulatory cores of a given set of modulons

``` r
library(ModulonCore)

# For this example we will use the TILs dataset included in ModulonCore

# Build the graph
graph = build.graph(network.TILs)
# Transform the graph into an igraph object
igraph = build.igraph(graph)
# Look for connected components
cc=find.connected.components(net=igraph,modulons = modulons.TILs)
# Check the statistical over-connectivity of the connected components
set.seed(123)
results.connectivity=connectivity( 
 net=igraph,
 cc=cc,
 min.size='5',
 dist.method="shortest.paths",
 nperm=100,
 vertex.attr="query",
 edge.attr="Weight",
 significance=0.05)
```

We can collect the identified regulatory cores as follows:

``` r
# Collect the regulatory cores IDs
regulatory.cores=core(results.connectivity=results.connectivity)
regulatory.cores
#>  [1] "2__cc.2"  "2__cc.3"  "2__cc.11" "2__cc.22" "2__cc.27" "2__cc.37"
#>  [7] "2__cc.42" "2__cc.46" "3__cc.3"  "4__cc.1"  "5__cc.1"
```

The regulatory cores can be printed out using the core.out() function:

``` r
core.out(net=graph,cc=cc,regulatory.cores=regulatory.cores,dir="./cores")
```

Plot the results of the connectivity analysis for the regulatory core of
modulon 5

``` r
# Connectivity analysis plot
plot(results.connectivity[['5']][['cc.1']])
```

<img src="man/figures/README-plot connectivity modulon 5 cc.1-1.png" width="100%" />
Show the pvalue (empirical p-value calculated using a Z-test) of the connectivity analysis
for the regulatory core of modulon 5

``` r
# pvalue
results.connectivity[['5']][['cc.1']]$pval
#> [1] 5.906808e-08
```

Plot the results of the connectivity analysis for the regulatory core of
modulon 3

``` r
plot(results.connectivity[['3']][['cc.3']])
```

<img src="man/figures/README-plot connectivity modulon 3 cc.3-1.png" width="100%" />

Show the pvalue (empirical p-value calculated using a Z-test) of the connectivity analysis
 for the regulatory core of modulon 3

``` r
# pvalue
results.connectivity[['3']][['cc.3']]$pval
#> [1] 1.471661e-06
```

## Author

<img src="man/figures/README-isaaccrespo_WEB_big.jpg" width="30%" />

Isaac Crespo, phD
Senior Computational Scientist, CHUV \| Department of
Oncology \| George Coukos group
Ludwig Institute for Cancer Research \|
Lausanne Branch
AGORA, Bugnon 25A, 1005 Lausanne, 4th floor, Room 026
<isaaccrespo@hotmail.com>
