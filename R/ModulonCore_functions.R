suppressMessages(require(igraph))
suppressMessages(require(SANTA))
suppressMessages(require(limma))

MODE = c("weak", "strong")
DIR = c("./")
MINSIZE=c('5')
WEIGHT=c('1')


#' @export 
build.graph = function(net,weight=WEIGHT){
  network.SCENIC = net
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


#' @export 
build.igraph = function(net){
  network= net
  # Network in igraph format
  network.igraph = igraph::graph_from_edgelist(as.matrix(network[c(1:2)]))
  E(network.igraph)$Weight = network$Weight
  return(network.igraph)
}

#' @title Connected Components Identification
#' @description Find all connected components within in modulon
#' @param net A dataframe with the network encoded in three columns with the following column names: 'Source', 'Interaction' and 'Target.
#' @param modulons A list.
#' @return A list with as many elements as modulons containing the constituent transcripton factors of all the connected components of each modulon
#' @details The names of the outpuT list are composed by....
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname jaccard
#' @export 
find.connected.components <- function(net, modulons,mode=MODE,dir=DIR) {
  
  mode <- match.arg(mode) # default = "weak"
  
  # Query
  TF.AUC.clusters.TILs=modulons
  
  # Extract all modulon subnetworks as edge list
  modulon.subnetworks = lapply(TF.AUC.clusters.TILs,function(x){
    subnetwork = network[(network$Source %in% x) & (network$Target %in% x),]
  })
  
  # Extract all modulon subnetworks as igraph object
  modulon.subnetworks.igraph = lapply(modulon.subnetworks,function(x){
    g = igraph::graph_from_edgelist(as.matrix(x[c(1:2)]))
    E(g)$Weights = x$Weight
    V(g)$query = 0
    g
  })
  # Extract all modulon subnetworks connected components
  modulon.subnetworks.connected.components = lapply(modulon.subnetworks.igraph,function(x){
    res=igraph::components(x, mode = c(mode))
    membership.df = as.data.frame(res$membership)
    connected.components.list = split(rownames(membership.df),membership.df[1])
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
#' @description Find all regulatory cores within in modulon
#' @param cc A list with all the connected components of each modulon
#' @return A list with as many elements as modulons containing the constituent transcripton factors of all the connected components of each modulon
#' @details The names of the outpuT list are composed by....
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname jaccard
#' @export 
connectivity <- function(cc, min.size=MINSIZE,dist.method,nperm,vertex.attr,edge.attr,significance) {
  modulon.subnetworks.connected.components=cc
  min.size=match.arg(min.size)
  modulon.subnetworks.connected.components.filtered = lapply(modulon.subnetworks.connected.components,function(x){
    x[lapply(x,length)>=as.integer(min.size)]
  })
  # SANTA analysis
  #dist.method="shortest.paths"
  #nperm=100
  #vertex.attr="query"
  #edge.attr="Weight"
  #significance=0.05
  vertex_names = igraph::get.vertex.attribute(network.igraph,'name')
  SANTA.INPUT = modulon.subnetworks.connected.components.filtered
  modulon.subnetworks.connected.components.SANTA = lapply(SANTA.INPUT,function(x){
    tmp = lapply(x,function(y){
      query = y
      V(network.igraph)$query = 0;
      
      for (j in 1: length(query)){
        for (i in 1: length(V(network.igraph))){
          if (identical(as.character(query[j]),vertex_names[i])){
            V(network.igraph)$query[i] = 1;
            print (paste(as.character(query[j])," equals to ",vertex_names[i]))
          }
          else{
            #print("nope");
          }
        }
      }
      
      g.clustered = SANTA::Knet(network.igraph, dist.method=dist.method, nperm=nperm, vertex.attr=vertex.attr,edge.attr=edge.attr)
      return(g.clustered)
    })
    return(tmp)
  })
  
  return (modulon.subnetworks.connected.components.SANTA)

}


#' @title Regulatory core
#' @description Find all regulatory cores within in modulon
#' @param cc A list with all the connected components of each modulon
#' @return A list with as many elements as modulons containing the constituent transcripton factors of all the connected components of each modulon
#' @details The names of the outpuT list are composed by....
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname jaccard
#' @export 
core <- function(modulon.subnetworks.connected.components.SANTA,significance=0.05) {
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
      print(pval.tmp)
      if(pval.tmp <= significance){
        regulatory.cores=c(regulatory.cores,paste(modulon.tmp,cc.tmp,sep = ('__')))
      }
    }
  }
  return (regulatory.cores)
  
}


#' @title Print regulatory core
#' @description Generate files with connected component networks and connected component elements
#' @param net A list with all the connected components of each modulon
#' @return A list with as many elements as modulons containing the constituent transcripton factors of all the connected components of each modulon
#' @details The names of the outpuT list are composed by....
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname jaccard
#' @export 
print.core = function(net,SANTA.INPUT,regulatory.cores,dir = DIR) {
  for(i in 1:length(regulatory.cores)){
    regulatory.core.tmp = regulatory.cores[i]
    modulon.tmp = limma::strsplit2(regulatory.core.tmp, '__')[,1]
    cc.tmp = limma::strsplit2(regulatory.core.tmp, '__')[,2]
    regulatory.core.elements.tmp=SANTA.INPUT[[modulon.tmp]][[cc.tmp]]
    write(file=paste(DIR,"Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.nodes.txt',sep = ''),regulatory.core.elements.tmp)
    network.tmp = net[((net$Source %in% regulatory.core.elements.tmp) & (net$Target %in% regulatory.core.elements.tmp)),]
    write.table(file=paste(dir,"Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.network.txt',sep = ''),network.tmp,quote = F,row.names = F,sep = "\t")
  }
  
}
