
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ModulonCore

## Title

ModulonCore

## Introduction

The goal of ModulonCore is to identify and characterize the modulon
regulatory cores or connected components with statistically supported
over-connectivity.  
The pipeline to achieve this goal can be described in 3 steps: 1)
modulon network assembly. The partial networks of the regulon
constituent elements are put together to assemble the modulon network,
where TFs are connected by transcriptional regulatory events inferred
from the data using the SCENIC regulon inference pipeline; 2) Connected
components identification. The modulon subnetwork (including only TF-TF
interactions) is further explored to find all the connected components;
3) Connectivity analysis. Every connected component is evaluated in
order to assess if their overconnectivity is statically supported using
the SANTA algorithm. The SANTA algorithm tests if the connectivity
within the connected component has a higher density than expected by
chance in random subnetworks of the same size subset from the global
TF-TF network derived from the assembly of all regulon partial networks.

<img src="./man/figures/ModulonCore_flowchart.jpg" width="100%" height="100%" style="display: block; margin: auto;" />

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

Loading the required packages:

``` r
library(ModulonCore)
library(igraph)
```

STEP 1: Modulon network assembly

``` r
# For this example we will use the TILs dataset included in ModulonCore
# Build the graph
graph = build.graph(network.TILs)
# Transform the graph into an igraph object
igraph = build.igraph(graph)
```

STEP 2: Connected components identification

``` r
# Look for connected components
cc=find.connected.components(net=igraph,modulons = modulons.TILs)
cc.table = tab.cc(cc)
head(cc.table)
#>          Modulon Connected.Component Size        Elements
#> 1__cc.7        1                cc.7    2     Atf6b, Bmyc
#> 1__cc.12       1               cc.12    2    Tcf7, Zfp973
#> 1__cc.13       1               cc.13    2     Rfx2, Trp53
#> 1__cc.16       1               cc.16    2 Gtf2ird1, Kdm5a
#> 1__cc.18       1               cc.18    2   Nr1d1, Zfp524
#> 1__cc.24       1               cc.24    2    Pbx3, Zfp143
```

STEP 3: Connectivity analysis

``` r
# Check the statistical over-connectivity of the connected components
set.seed(123)
results.connectivity = connectivity( 
 net=igraph,
 cc=cc,
 min.size='5',
 dist.method="shortest.paths",
 nperm=100,
 vertex.attr="query",
 edge.attr="Weight",
 significance=0.05)
```

The identified regulatory cores can be retrieved as follows:

``` r
# Collect the regulatory cores IDs
regulatory.cores = core(results.connectivity=results.connectivity)
regulatory.cores
#>  [1] "2__cc.2"  "2__cc.3"  "2__cc.11" "2__cc.22" "2__cc.27" "2__cc.37"
#>  [7] "2__cc.42" "2__cc.46" "3__cc.3"  "4__cc.1"  "5__cc.1"
```

This information can be added to the connected components table as
follows:

``` r
# Collect the regulatory cores IDs
cc.table$Core = ifelse(rownames(cc.table) %in% regulatory.cores,T,F)
head(cc.table)
#>          Modulon Connected.Component Size        Elements  Core
#> 1__cc.7        1                cc.7    2     Atf6b, Bmyc FALSE
#> 1__cc.12       1               cc.12    2    Tcf7, Zfp973 FALSE
#> 1__cc.13       1               cc.13    2     Rfx2, Trp53 FALSE
#> 1__cc.16       1               cc.16    2 Gtf2ird1, Kdm5a FALSE
#> 1__cc.18       1               cc.18    2   Nr1d1, Zfp524 FALSE
#> 1__cc.24       1               cc.24    2    Pbx3, Zfp143 FALSE
```

The regulatory core subnetworks can be printed out using the core.out()
function:

``` r
regulatory.core.networks=core.out(net=graph,cc=cc,regulatory.cores=regulatory.cores,dir="./cores")
head(regulatory.core.networks[['5__cc.1']])
#>    Source Interaction  Target
#> 1   Nfat5          -> Bhlhe40
#> 2     Rel          ->    Irf8
#> 3     Rel          ->  Pou2f2
#> 4     Rel          ->    Relb
#> 5     Rel          ->   Nfkb2
#> 6 Bhlhe40          ->   Hif1a
```

Plot the results of the connectivity analysis for the regulatory core of
modulon 5

``` r
# Explore results for modulon 5

# Connectivity analysis plot
plot(results.connectivity[['5']][['cc.1']])
```

<img src="man/figures/README-plot connectivity modulon 5 cc.1-1.png" width="100%" style="display: block; margin: auto;" />
The connectivity of the connected component cc.1 from modulon 5 is
evaluated using the area under the Knet-function (AUK) by comparison of
the connectivity of the regulatory core against random subnetworks of
the same size derived from the global transcriptional network inferred
with SCENIC. An empirical p-value is calculated using a Z-test

Show the pvalue (empirical p-value calculated using a Z-test) of the
connectivity analysis for the regulatory core of modulon 5

``` r
# pvalue
results.connectivity[['5']][['cc.1']]$pval
#> [1] 5.906808e-08
```

Display the network of a regulatory core:

``` r
# Extract the subnetwork
net = igraph::subgraph(igraph,vids = cc[['5']][['cc.1']])
net = igraph::simplify(net,remove.multiple = T)
l <- igraph::layout_in_circle(net)
# Calculate node outdegree to use define the node size
deg <- igraph::degree(net, mode="out")
plot(net, layout=l,vertex.size=5 + deg*2,vertex.label.color="black", vertex.label.dist=c(-3,3), edge.arrow.size=.4,vertex.color = adjustcolor("SkyBlue2", alpha.f = .5))
```

<img src="man/figures/README-plot network modulon 5 cc.1-1.png" width="100%" style="display: block; margin: auto;" />

The node size is related with the outdegree; bigger nodes have a greater
outdegree.

## Author

Isaac Crespo

<img src="./man/figures/README-isaaccrespo_WEB_big.jpg" width="30%" />

Isaac Crespo, phD  
Senior Computational Scientist  
CHUV \| Department of Oncology \| George Coukos group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne, 4th floor, Room 026  
<isaaccrespo@hotmail.com>
