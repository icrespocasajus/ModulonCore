
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
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(SANTA)
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
#> computing graph distance matrix...
#> done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
#> computing graph distance matrix... done
#> computing graph distance bins... done
#> computing the clustering of the 'query' weights using 100 permutations...  done
```

We can collect the identified regulatory cores as follows:

``` r
# Collect the regulatory cores IDs
regulatory.cores=core(results.connectivity=results.connectivity)
#> [1] "Modulon:  2    Connected component:  cc.2    pvalue:  8.23547602166656e-09"
#> [1] "Modulon:  2    Connected component:  cc.3    pvalue:  2.95335396689315e-08"
#> [1] "Modulon:  2    Connected component:  cc.11    pvalue:  1.91201113434279e-07"
#> [1] "Modulon:  2    Connected component:  cc.22    pvalue:  2.63423192496079e-13"
#> [1] "Modulon:  2    Connected component:  cc.27    pvalue:  5.29418975606127e-07"
#> [1] "Modulon:  2    Connected component:  cc.37    pvalue:  9.96985726777097e-11"
#> [1] "Modulon:  2    Connected component:  cc.42    pvalue:  1.85924805953426e-23"
#> [1] "Modulon:  2    Connected component:  cc.46    pvalue:  2.18297793149058e-14"
#> [1] "Modulon:  3    Connected component:  cc.3    pvalue:  1.47166073436368e-06"
#> [1] "Modulon:  4    Connected component:  cc.1    pvalue:  1.28168023531566e-08"
#> [1] "Modulon:  5    Connected component:  cc.1    pvalue:  5.90680782426239e-08"
regulatory.cores
#>  [1] "2__cc.2"  "2__cc.3"  "2__cc.11" "2__cc.22" "2__cc.27" "2__cc.37"
#>  [7] "2__cc.42" "2__cc.46" "3__cc.3"  "4__cc.1"  "5__cc.1"
```

The regulatory cores can be printed out using the core.out() function:

``` r
core.out(net=graph,cc=cc,regulatory.cores=regulatory.cores,dir="./cores")
#> Warning in dir.create(dir): './cores' already exists
```

Plot the results of the connectivity analysis for the regulatory core of
modulon 5

``` r
# Explore results for modulon 5

# Connectivity analysis plot
plot(results.connectivity[['5']][['cc.1']])
```

<img src="man/figures/README-plot connectivity modulon 5 cc.1-1.png" width="100%" />
Show the pvalue (empiric z-score) of the connectivity analysis for the
regulatory core of modulon 5

``` r
# pvalue
results.connectivity[['5']][['cc.1']]$pval
#> [1] 5.906808e-08
```

Plot the results of the connectivity analysis for the regulatory core of
modulon 3

``` r
# Explore results for modulon 3
plot(results.connectivity[['3']][['cc.3']])
```

<img src="man/figures/README-plot connectivity modulon 3 cc.3-1.png" width="100%" />

Show the pvalue (empiric z-score) of the connectivity analysis for the
regulatory core of modulon 3

``` r
# pvalue
results.connectivity[['3']][['cc.3']]$pval
#> [1] 1.471661e-06
```

## Author

Isaac Crespo

<img src="man/figures/README-pressure-1.png" width="100%" />

Isaac Crespo, phD Senior Computational Scientist, CHUV \| Department of
Oncology \| George Coukos group Ludwig Institute for Cancer Research \|
Lausanne Branch AGORA, Bugnon 25A, 1005 Lausanne, 4th floor, Room 38
<isaaccrespo@hotmail.com>
