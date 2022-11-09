suppressMessages(require(igraph))
suppressMessages(require(SANTA))

MODE = c("weak", "strong")
DIR = c("./")
MINSIZE=c('5')
WEIGHT=c('1')




#' @title Build the graph
#' @description Find all regulatory regulatory events between TFs and add weights for the regulatory interactions
#' @param net A dataframe with a network encoded in 3 columns: 'Source','Interaction','Target'.
#' @param weight A vector with the weight of each interaction. If none is provided, it will be assumed the same weight (1) for all the interactions
#' @return A dataframe with the network encoded in 3 columns:  'Source','Target','Weight'.
#' @details The resulting network only contains regulatory events between TFs
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  build.graph(network.TILs)
#'  }
#' }
#' @rdname build.graph
#' @export 
build.graph = function(net,weight=WEIGHT){
  network.SCENIC = net
  #Remove duplicated edges
  network.SCENIC = unique(network.SCENIC)
  #Remove self-loops
  network.SCENIC = network.SCENIC[!(network.SCENIC$Source == network.SCENIC$Target),]
  all.TF = unique(network.SCENIC$Source)
  #all.targets = unique(network.SCENIC$Target)
  all.targets = intersect(unique(network.SCENIC$Target),all.TF)
  all.targets=all.TF
  network.SCENIC = network.SCENIC[network.SCENIC$Target %in% all.TF,]
  network=data.frame(Source = network.SCENIC$Source,Target=network.SCENIC$Target,Weight = as.numeric(weight))
  return(network)
}



#' @title Build the igraph
#' @description Covert the input network into an igraph object
#' @param net A dataframe with the network encoded in 3 columns:  'Source','Target','Weight'.
#' @return A network as an igraph object
#' @details The resulting network only contains regulatory events between TFs
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  graph = build.graph(network.TILs)
#'  build.igraph(graph)
#'  }
#' }
#' @rdname build.igraph
#' @export 
build.igraph = function(net){
  network = net
  # Network in igraph format
  network.igraph = igraph::graph_from_edgelist(as.matrix(network[c(1:2)]))
  igraph::E(network.igraph)$Weight = network$Weight
  network.igraph = igraph::simplify(network.igraph,remove.multiple = T,remove.loops = T,edge.attr.comb='first')
  return(network.igraph)
}

#' @title Connected Components Identification
#' @description Find all connected components within modulons or any cluster of genes
#' @param net A network as an igraph object
#' @param modulons A list with as many elements as modulons/clusters containing the constituent elements.
#' @return A list with as many elements as modulons/clusters containing  all the connected components
#' @details The names of connected components are encoded as cc.X, being 'X' the ID number of the connected component
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  graph = build.graph(network.TILs)
#'  igraph = build.igraph(graph)
#'  find.connected.components(net=igraph,modulons = modulons.TILs)
#'  }
#' }
#' @rdname find.connected.components
#' @export 
find.connected.components <- function(net, modulons,mode="weak",dir=DIR) {
  
  # Query
  TF.AUC.clusters.TILs = modulons
  
  # Extract all modulon subnetworks as igraph object
  modulon.subnetworks.igraph = lapply(TF.AUC.clusters.TILs,function(x){
    subnetwork = igraph::subgraph(net,vids = intersect(names(igraph::V(net)),x))
  })
 
  # Extract all modulon subnetworks connected components
  modulon.subnetworks.connected.components = lapply(modulon.subnetworks.igraph,function(x){
    res = igraph::components(x, mode = c(mode))
    membership.df = as.data.frame(res$membership)
    size.df = as.data.frame(res$csize)
    connected.components.list = split(rownames(membership.df),membership.df[1])
    connected.components.list = connected.components.list[size.df[names(connected.components.list),1] > 1]
    return(connected.components.list)
  })
  # Rename connected components
  modulon.subnetworks.connected.components = lapply(modulon.subnetworks.connected.components,function(x){
    names(x)=paste('cc.',names(x),sep = '')
    return(x)
  })
  saveRDS(file=paste(dir,"Modulon_Connected_Components.Rds",sep = ''),modulon.subnetworks.connected.components)
  return (modulon.subnetworks.connected.components)
}



#' @title SANTA analysis
#' @description Evaluate the statistical support for connected components over-connectivity
#' @param cc A list with all the connected components of each modulon
#' @return A list with as many elements as modulons containing the results of the SANTA analysis
#' @details The results can be plotted as follows: plot(x[[][[]]])
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  graph = build.graph(network.TILs)
#'  igraph = build.igraph(graph)
#'  cc=find.connected.components(net=igraph,modulons = modulons.TILs)
#'  connectivity( net=igraph, cc=cc, min.size='5',
#'  dist.method="shortest.paths",
#'  nperm=100,
#'  vertex.attr="query",
#'  edge.attr="Weight",
#'  significance=0.05)
#'  }
#' }
#' @rdname connectivity
#' @export 
connectivity <- function(net,cc, min.size,dist.method,nperm,vertex.attr,edge.attr,significance) {
  modulon.subnetworks.connected.components=cc
  modulon.subnetworks.connected.components.filtered = lapply(modulon.subnetworks.connected.components,function(x){
    x[lapply(x,length)>=as.integer(min.size)]
  })
  vertex_names = igraph::get.vertex.attribute(net,'name')
  SANTA.INPUT = modulon.subnetworks.connected.components.filtered
  modulon.subnetworks.connected.components.SANTA = lapply(SANTA.INPUT,function(x){
    tmp = lapply(x,function(y){
      query = y
      igraph::V(net)$query = 0
      
      for (j in 1: length(query)){
        for (i in 1: length(igraph::V(net))){
          if (identical(as.character(query[j]),vertex_names[i])){
            igraph::V(net)$query[i] = 1;
            print (paste(as.character(query[j]),",equals to ",vertex_names[i]))
          }
          else{
            #print("nope");
          }
        }
      }
      
      g.clustered = SANTA::Knet(net, dist.method=dist.method, nperm=nperm, vertex.attr=vertex.attr,edge.attr=edge.attr)
      return(g.clustered)
    })
    return(tmp)
  })
  
  return (modulon.subnetworks.connected.components.SANTA)

}


#' @title Regulatory core identification
#' @description Evaluate the results of the connectivity analysis and collect connected components with statistically supported over-connectivity
#' @param results.connectivity A list with the output generated by the function connectivity().
#' @return A list with all the regulatory cores of each modulon
#' @details Modulons may have none, one or more than one regulatory cores
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  graph = build.graph(network.TILs)
#'  igraph = build.igraph(graph)
#'  cc=find.connected.components(net=igraph,modulons = modulons.TILs)
#'  results.connectivity=connectivity( net=igraph, cc=cc, min.size='5',
#'  dist.method="shortest.paths",
#'  nperm=100,
#'  vertex.attr="query",
#'  edge.attr="Weight",
#'  significance=0.05)
#'  core(results.connectivity=results.connectivity)
#' }
#' @rdname core
#' @export 
core <- function(results.connectivity,significance=0.05) {
  modulon.subnetworks.connected.components.SANTA = results.connectivity
  # Regulatory Cores
  # Check SANTA statistical significance
  regulatory.cores = c()
  for(i in 1:length(names(modulon.subnetworks.connected.components.SANTA))){
    modulon.tmp = names(modulon.subnetworks.connected.components.SANTA)[i]
    if(length(modulon.subnetworks.connected.components.SANTA[[modulon.tmp]]) == 0){next}
    for(j in 1:length(names(modulon.subnetworks.connected.components.SANTA[[modulon.tmp]]))){
      cc.tmp = names(modulon.subnetworks.connected.components.SANTA[[modulon.tmp]])[j]
      # print(as.numeric(modulon.subnetworks.connected.components.SANTA[[modulon.tmp]][[cc.tmp]]['pval'] ))
      pval.tmp = as.numeric(modulon.subnetworks.connected.components.SANTA[[modulon.tmp]][[cc.tmp]]['pval'] )
      #print(paste('Modulon: ',modulon.tmp,'   Connected component: ',cc.tmp,'   pvalue: ',pval.tmp))
      if(pval.tmp <= significance){
        regulatory.cores=c(regulatory.cores,paste(modulon.tmp,cc.tmp,sep = ('__')))
      }
    }
  }
  return (regulatory.cores)
  
}


#' @title Print out regulatory core analysis results
#' @description Generate files with regulatory core networks and constituent elements.
#' @param net Network built with build.graph()
#' @param cc Connected components generated with find.connected.components()
#' @param regulatory.cores Regulatory cores as the output of core()
#' @param dir Output directory
#' @return Files with regulatory core networks and constituent elements printed out to the output directory, and a list with all the regulatory core networks
#' @details The name of each file contains the corresponding modulon and connected component
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #'  graph = build.graph(network.TILs)
#'  igraph = build.igraph(graph)
#'  cc=find.connected.components(net=igraph,modulons = modulons.TILs)
#'  results.connectivity=connectivity( net=igraph, cc=cc, min.size='5',
#'  dist.method="shortest.paths",
#'  nperm=100,
#'  vertex.attr="query",
#'  edge.attr="Weight",
#'  significance=0.05)
#'  regulatory.cores=core(results.connectivity=results.connectivity)
#'  regulatory.core.networks=core.out(net=graph,cc=cc,regulatory.cores=regulatory.cores,dir="./cores")
#'  }
#' }
#' @rdname core.out
#' @export 
core.out = function(net,cc,regulatory.cores,dir = "./") {
  SANTA.INPUT=cc
  dir.create(dir)
  regulatory.core.networks = list()
  for(i in 1:length(regulatory.cores)){
    regulatory.core.tmp = regulatory.cores[i]
    modulon.tmp = strsplit2(regulatory.core.tmp, '__')[,1]
    cc.tmp = strsplit2(regulatory.core.tmp, '__')[,2]
    regulatory.core.elements.tmp=SANTA.INPUT[[modulon.tmp]][[cc.tmp]]
    write(file=paste(dir,"/Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.nodes.txt',sep = ''),regulatory.core.elements.tmp)
    network.tmp = net[((net$Source %in% regulatory.core.elements.tmp) & (net$Target %in% regulatory.core.elements.tmp)),]
    network.tmp = data.frame(Source=network.tmp[,1],Interaction='->',Target=network.tmp[,2])
    write.table(file=paste(dir,"/Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.network.txt',sep = ''),network.tmp,quote = F,row.names = F,sep = "\t")
    regulatory.core.networks[[regulatory.core.tmp]] = network.tmp
  }
  return(regulatory.core.networks)
}



#' @title Range between 0 and 1
#' @description Normalize absolute values of a numeric vector
#' @param x A numeric vector
#' @return Normalized values (x-min(x))/(max(x)-min(x))
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' range01(c('-50','-5','5','50','100')
#'  }
#' }
#' @rdname range01
#' @export 
range01 <- function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))}


#' @title Split Composite Names
#' @description Split a vector of composite names into a matrix of simple names.
#' @param x character vector
#' @param split character to split each element of vector on, see strsplit
#' @param ... other arguments are passed to strsplit
#' @return Normalized values (x-min(x))/(max(x)-min(x))
#' @details This function is the same as strsplit except that the output value is a matrix instead of a list. 
#' The first column of the matrix contains the first component from each element of x, the second column contains the second components etc.
#'  The number of columns is equal to the maximum number of components for any element of x.
#'   The motivation for this function in the limma package is handle input columns which are composites of two or more annotation fields.
#' @examples 
#' \dontrun{
#' if(interactive()){
#' x = c("2__cc.1","2__cc.2"); strsplit2(x,split="__")
#'  }
#' }
#' @rdname strsplit2
#' @export 
strsplit2 = function (x, split, ...) 
{
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i]) 
      out[i, 1:nc[i]] <- s[[i]]
  }
  out
}


#' @title Connected components table generation
#' @description Generate a table with all the connected components and their constituent elements
#' @param cc A list with all the connected components of each modulon.
#' @return A table of 4 columns
#' @details 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  cc = cc.TILs
#'  tab.cc(cc)
#'  }
#' }
#' @rdname tab.cc
#' @export 
tab.cc <- function(cc,cores=NULL) {
  all.cc = c()  
  for(i in 1:length(names(cc))) {
    modulon.tmp = names(cc)[i]
    cc.tmp = names(cc[[modulon.tmp]])
    all.cc = c(all.cc,paste(modulon.tmp,cc.tmp,sep = '__'))
    
  }   
  query.df = data.frame(Modulon=ModulonCore::strsplit2(all.cc,'__')[,1],Connected.Component=ModulonCore::strsplit2(all.cc,'__')[,2])
  rownames(query.df)=all.cc
  
  all.sizes = c()    
  cc.string = c()
  for(i in 1:nrow(query.df )){
    mod.tmp = query.df[i,'Modulon']
    cc.tmp = query.df[i,'Connected.Component']
    cc.string = c(cc.string,paste(sort(cc[[mod.tmp]][[cc.tmp]]),sep = ', ',collapse = ', '))
    all.sizes = c(all.sizes,length(cc[[mod.tmp]][[cc.tmp]]))
  }
  table = query.df
  table$Size = all.sizes
  table$Elements = cc.string
  return(table)
}