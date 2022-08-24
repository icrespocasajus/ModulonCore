suppressMessages(require(igraph))
suppressMessages(require(SANTA))
suppressMessages(require(limma))

MODE = c("weak", "strong")
DIR = c("./")
MINSIZE=c('5')
WEIGHT=c('1')




#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @param weight PARAM_DESCRIPTION, Default: WEIGHT
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname build.graph
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



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[igraph]{graph_from_edgelist}}, \code{\link[igraph]{E}}
#' @rdname build.igraph
#' @export 
#' @importFrom igraph graph_from_edgelist E
build.igraph = function(net){
  network = net
  # Network in igraph format
  network.igraph = igraph::graph_from_edgelist(as.matrix(network[c(1:2)]))
  igraph::E(network.igraph)$Weight = network$Weight
  return(network.igraph)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @param modulons PARAM_DESCRIPTION
#' @param mode PARAM_DESCRIPTION, Default: MODE
#' @param dir PARAM_DESCRIPTION, Default: DIR
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[igraph]{subgraph}}, \code{\link[igraph]{V}}, \code{\link[igraph]{component_distribution}}
#' @rdname find.connected.components
#' @export 
#' @importFrom igraph subgraph V components
find.connected.components <- function(net, modulons,mode=MODE,dir=DIR) {
  
  mode <- match.arg(mode) # default = "weak"

  # Query
  TF.AUC.clusters.TILs=modulons
  
  # Extract all modulon subnetworks as igraph object
  modulon.subnetworks.igraph = lapply(TF.AUC.clusters.TILs,function(x){
    subnetwork = igraph::subgraph(igraph,vids = intersect(names(igraph::V(net)),x))
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



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @param cc PARAM_DESCRIPTION
#' @param min.size PARAM_DESCRIPTION
#' @param dist.method PARAM_DESCRIPTION
#' @param nperm PARAM_DESCRIPTION
#' @param vertex.attr PARAM_DESCRIPTION
#' @param edge.attr PARAM_DESCRIPTION
#' @param significance PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[igraph]{vertex_attr}}, \code{\link[igraph]{V}}
#'  \code{\link[SANTA]{Knet}}
#' @rdname connectivity
#' @export 
#' @importFrom igraph get.vertex.attribute V
#' @importFrom SANTA Knet
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
      igraph::V(net)$query = 0;
      
      for (j in 1: length(query)){
        for (i in 1: length(V(net))){
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param results.connectivity PARAM_DESCRIPTION
#' @param significance PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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
      print(paste('Modulon: ',modulon.tmp,'   Connected component: ',cc.tmp,'   pvalue: ',pval.tmp))
      if(pval.tmp <= significance){
        regulatory.cores=c(regulatory.cores,paste(modulon.tmp,cc.tmp,sep = ('__')))
      }
    }
  }
  return (regulatory.cores)
  
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net PARAM_DESCRIPTION
#' @param cc PARAM_DESCRIPTION
#' @param regulatory.cores PARAM_DESCRIPTION
#' @param dir PARAM_DESCRIPTION, Default: './'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[limma]{strsplit2}}
#' @rdname core.out
#' @export 
#' @importFrom limma strsplit2
core.out = function(net,cc,regulatory.cores,dir = "./") {
  SANTA.INPUT=cc
  dir.create(dir)
  for(i in 1:length(regulatory.cores)){
    regulatory.core.tmp = regulatory.cores[i]
    modulon.tmp = limma::strsplit2(regulatory.core.tmp, '__')[,1]
    cc.tmp = limma::strsplit2(regulatory.core.tmp, '__')[,2]
    regulatory.core.elements.tmp=SANTA.INPUT[[modulon.tmp]][[cc.tmp]]
    write(file=paste(dir,"Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.nodes.txt',sep = ''),regulatory.core.elements.tmp)
    network.tmp = net[((net$Source %in% regulatory.core.elements.tmp) & (net$Target %in% regulatory.core.elements.tmp)),]
    write.table(file=paste(dir,"Regulatory.Core.Modulon.",modulon.tmp,'_',cc.tmp,'.network.txt',sep = ''),network.tmp,quote = F,row.names = F,sep = "\t")
  }
  
}




#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname jaccard
#' @export 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname redundancy
#' @export 
redundancy <- function(a, b) {
  intersection = length(intersect(a, b))
  min = min(length(a),length(b))
  return (intersection/min)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname range01
#' @export 
range01 <- function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))}


