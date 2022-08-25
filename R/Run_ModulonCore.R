net = readRDS(file="../Results/network.TILs.Rds")
modulons=readRDS(file="../Results/TF.AUC.clusters.TILs.Rds")
dir= "../"


network = build.graph(net=net)
network.igraph = build.igraph(net = network)

cc=find.connected.components(net=network.igraph,modulons=modulons) 



modulon.subnetworks.connected.components.SANTA=connectivity( cc, min.size='5',
              dist.method="shortest.paths",
              nperm=100,
              vertex.attr="query",
              edge.attr="Weight",
              significance=0.05)


regulatory.cores=core(modulon.subnetworks.connected.components.SANTA,significance=0.05)

print.core(net = network,SANTA.INPUT=cc,regulatory.cores=regulatory.cores) 