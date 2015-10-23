######################
####  6 May 2015 ####
######################
##
## Nous implementons l'idée d'un demi cercle fractal pour résoudre les
## problèmes liés aux multifurcations
##
## Nous utilisons un code proche de celui de fractal.phylotrees.R
##
##
## le début du demi cercle a pour coordonnées 0,0
## 


require(ape)
require(phangorn)
require(scales)

rad<-function(ang) {
    return((ang*pi)/180)
}
init.angle<-90
opens<-180 ##controls the opening of the flowers (in degrees)


edges<-read.table("edges.txt")
species<-scan(what="character", file="species.txt", sep="\n")
nodes<-scan(what="character", file="nodes.txt", sep="\n")
NN.size<-scan(file="sizes.txt")


nbt<-length(species)

if (length(edges[edges[,1]==min(edges[,1]),2])==1) {
    edges<-edges[2:nrow(edges),]
    edges[,1]<-edges[,1]-1
    edges[edges[,2]>nbt,2]<-edges[edges[,2]>nbt,2]-1
    nodes<-nodes[2:length(nodes)]
    NN.size<-NN.size[2:length(NN.size)]
}
nbn<-length(nodes)

racine<-nbt+1
##we go through the edges rows and decide the x and y position of each node each time

## we first compute the number of descendants of each node, in the order they appear in edges[,1]
node.refs<-unique(edges[,1])
##
##
node.size<-c(rep(1,nbt),NN.size)
##node.size<-sqrt(node.size) ##could be log to be even smoother

##DOES NOT WORK FOR BIG TREES: node.size<-c(rep(1,nbt),unlist(lapply(node.refs, function(x,T) length(Descendants(T,x)[[1]]), T=tr)))


##node.size<-log(node.size)+1


##for equal stuff,let's fake everything is symmetrical
##node.size<-rep(1, nbt+nbn)


nodes.all<-c(node.refs, 1:nbt)
COORD<-matrix("na", ncol=5, nrow=length(nodes.all))
COORD<-as.data.frame(COORD)
colnames(COORD)<-c("istip","x","y","alpha","d")
COORD$istip<-c(rep("yes",nbt),rep("no",nbn)) ##'istip' will contain at the end yes or no depending whether the node is a tip or not. 
COORD$x<-0 ##x coordiante of nodes
COORD$y<--30 ##y coordinates of nodes
COORD$alpha<-0 ##angle of the half circle
COORD$alpha[racine]<-init.angle
COORD$d<-0 # ray of the half circle.
COORD$d[racine]<-60 ##the is the ray of the initial half circle
COORD$nbdesc<-0 ##number of descendants at each node
COORD$dist2root<-0
cpt<-0



dist2dad<-function(alpha, D) {
    ##gives the distance to the parental node (dist) according to the
    ##angle (representing the importance in terms of nb of desc.)
    ##newD is the new ray of the half-circle
    newD<-(D*sin(rad(alpha))/cos(rad(alpha)))/(1+(sin(rad(alpha))/cos(rad(alpha))))
    dist<-D-newD
    return(c(dist,newD))
}

nodes.all2<-nodes.all[2:length(nodes.all)]

pb<-txtProgressBar(min=0, max=length(nodes.all2), style=3)
    
while(length(nodes.all2)>0) {
    setTxtProgressBar(pb, length(nodes.all2))
    nod.current<-nodes.all2[1]    
    dad<-edges[edges[,2]==nod.current,1]
    brothers<-edges[edges[,1]==dad,2]
    nbdesc.bro<-node.size[brothers]
    ##we can order the brothers by their nb of descendants.
    ##brothers<-brothers[order(nbdesc.bro)]
    ##nbdesc.bro<-sort(nbdesc.bro)

    D.dad<-COORD$d[dad]
    Angles.brothers<-opens*(nbdesc.bro/sum(nbdesc.bro))/2
    distances.bro<-sapply(Angles.brothers, dist2dad, D=D.dad)
    ##first row of distances.bro is the absolute distance to the paternal node
    ##second row is the new ray of the half circle
    dists.bro<-distances.bro[1,]
    D.bro<-distances.bro[2,]
    
    Angles.of.desc<-cumsum(rep(opens*(nbdesc.bro/sum(nbdesc.bro))/2, each=2))
    Angles.of.desc2<-Angles.of.desc[seq(1, 2*length(brothers)-1, by=2)]
    alpha.bro<-Angles.of.desc2-((opens/2)-COORD$alpha[dad])
    
    x.dad<-COORD$x[dad]
    y.dad<-COORD$y[dad]
    x.bro<-x.dad+dists.bro*cos(rad(alpha.bro))
    y.bro<-y.dad+dists.bro*sin(rad(alpha.bro))
    
    COORD[brothers,]$x<-x.bro
    COORD[brothers,]$y<-y.bro
    COORD[brothers,]$alpha<-alpha.bro
    COORD[brothers,]$d<-D.bro
    ##and we update the distance to root for the current brothers
    COORD$dist2root[brothers]<-COORD$dist2root[dad]+1
    nodes.all2<-setdiff(nodes.all2, brothers)
}
close(pb)

COORD$names<-c(species, nodes)
COORD$nbdesc<-node.size


###We compute the Normalized dist2root:
COORD$dist2root.norm<-COORD$dist2root/max(COORD$dist2root)


halfCircle <- function(x,y,r,start=0,end=pi+start,nsteps=30,plot=TRUE,...){
   rs <- seq(start,end,len=nsteps)
   xc <- x+r*cos(rs)
   yc <- y+r*sin(rs)
   if (plot) polygon(xc,yc,...)
   return(cbind(xc,yc))
}


##ellipse that, with the halfcircle, will form the mediator-shaped clade polygon.
ellipse<-function(x,y,r, alpha, nsteps=30, plot=TRUE,...) {
    start<-0
    end=pi+start
    rs <- seq(start,end,len=nsteps)
    a<-r
    b<-r/4
    xs<-a*cos(rs)
    ys<-b*sin(rs)
    ##rotation 
    xs2<-x+(xs*cos(alpha)-ys*sin(alpha))
    ys2<-y+(xs*sin(alpha)+ys*cos(alpha))
    if (plot) polygon(xs2,ys2,...)
    return(cbind(xs2,ys2))
}

HalfCircPlusEllips<-function(x,y,r,alpha, start=0, end=pi+start,nsteps=30,plot=FALSE,...) {
    circ<-halfCircle(x,y,r,start,end, nsteps, plot,...)
    elli<-ellipse(x,y,r,alpha,nsteps,plot,...)
    ALL<-rbind(circ,elli)
    ALL
}


###we add $id to the COORD matrix
COORD$id<-1:(nbt+nbn)


##We rermove underscores in COORD$names
COORD$names<-gsub("_"," ",COORD$names)


X<-COORD$x[(nbt+1):nrow(COORD)]
Y<-COORD$y[(nbt+1):nrow(COORD)]
r<-COORD$d[(nbt+1):nrow(COORD)]
start<-rad(COORD$alpha[(nbt+1):nrow(COORD)])-pi/2
alphas<-rad(COORD$alpha[(nbt+1):nrow(COORD)])+pi/2
##UNCOMENT NEXT LINE TO PLOT
plot(COORD$x, COORD$y,type="n", xlim=c(-100,100), ylim=c(-100,100), frame=FALSE, axes=FALSE, xlab="", ylab="")
Polygons<-NULL
ID1<-nbt+nbn
ID<-ID1
len<-length(X)
pb<-txtProgressBar(min=0, max=len, style=3)
for (i in 1:len) {
    setTxtProgressBar(pb, len-i)
    ID<-ID+1
    Polygons<-rbind(Polygons, cbind(ID,HalfCircPlusEllips(X[i], Y[i], r[i],alphas[i],start[i], plot=TRUE, border=NA, col=alpha("red",0.06))))
    ##CHANGE TO  plot=TRUE TO PLOT
}
close(pb)

Polygons<-as.data.frame(Polygons)
IDend<-(ID+nrow(Polygons))
Polygons$nod<-(ID+1):IDend

##UNCOMENT NEXT LINE TO PLOT
segments(COORD$x[edges[,1]],COORD$y[edges[,1]],COORD$x[edges[,2]], COORD$y[edges[,2]], lwd=0.8)
##UNCOMENT NEXT LINE TO PLOT TEXT (UGLY)
##text(COORD$x[COORD$istip=="yes"], COORD$y[COORD$istip=="yes"], labels=COORD$names[COORD$istip=="yes"], cex=1)
##this is the segment leading to the root: goes from the root to a point that is 10% under the root (relative to initial ray of the biggest halfcircle)

##the distance to the root is called preRootDist
preRootDist<-6
alpha<-rad(init.angle+180)
x.pre<-COORD$x[nbt+1]+preRootDist*cos(alpha)
y.pre<-COORD$y[nbt+1]+preRootDist*sin(alpha)
##UNCOMENT NEXT LINE TO PLOT
##segments(COORD$x[nbt+1], COORD$y[nbt+1], x.pre, y.pre)


###GET PolyGon Center
##for this we take the x,y of the node, the d and the alpha. We divide x and y by 2 and rotate 90 degrees arounbd original point.
  ##rotation 
polygon.x<-cos(rad(COORD$alpha))*(COORD$d/3)+COORD$x
polygon.y<-sin(rad(COORD$alpha))*(COORD$d/3)+COORD$y
COORD<-cbind(COORD, polygon.x, polygon.y)


###NUMBER OF NECESSARY ZOOM LEVEL TO SEE EVERYTING
##minid is the min d values in COORD$d (the ray of the smallest half circle)
##100 is the ray of the biggest one; We want a zoom level so that the smallest half circle ends up as big as the biggest one.
minid<-min(COORD$d)
nbzoom<-ceiling(log2(100/minid))
##this number will be written as a comment in the osm file
print(paste("INFO: minimum ", nbzoom, " zoom levels required", sep=""))
##we calculate this min zoom (called "visible") for each COORD$d. It corresponds to the zoom necessary to see each the diamleter equal to 30% of the max dimameter at zoom 0
zoomview<-ceiling(log2(30/COORD$d))
zoomview[zoomview<0]<-0
COORD<-cbind(COORD, zoomview)

##############################
##############################
#####   WRITE OSM FILE   #####
##############################
##############################

###header;
print("writing file tree.osm with osm-formatted data suitable for osm2pgsql and mapnik...")

fileout<-"tree.osm"

cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", file=fileout)
cat("<osm version=\"0.6\" generator=\"HomeMade R Generator for Phylogeny\" copyright=\"OpenStreetMap and D.M. de Vienne\" attribution=\"http://www.openstreetmap.org/copyright\" license=\"http://opendatacommons.org/licenses/odbl/1-0/\">\n", file=fileout, append=TRUE)
cat(paste("<!-- <The number of required Zoom Levels to view this TreeMap is: ",nbzoom," > -->\n",sep=""), file=fileout, append=TRUE)
cat (" <bounds minlat=\"-90\" minlon=\"-180\" maxlat=\"90\" maxlon=\"180\"/>\n", file=fileout, append=TRUE)

##nodes (tips and nodes);

nodeinfo<-paste(" <node id=\"",COORD$id,"\" visible=\"true\" lat=\"",COORD$y,"\" lon=\"",COORD$x,"\">", sep="")
##tag1<-paste("  <tag k=\"nodal\" v=\"",COORD$dist2root.norm,"\"/>", sep="")
tag1<-paste("  <tag k=\"nodal\" v=\"",COORD$dist2root.norm,"\"/>", sep="")
tag2<-paste("  <tag k=\"name\" v=\"",COORD$names,"\"/>", sep="")
tag3<-paste("  <tag k=\"tip\" v=\"",COORD$istip,"\"/>", sep="")
tag4<-paste("  <tag k=\"zoomview\" v=\"",COORD$zoomview,"\"/>", sep="")
last<-rep(" </node>",nbt+nbn)
nodeblock<-paste(nodeinfo,tag1,tag2,tag3,tag4,last,sep="\n")
cat(nodeblock, sep="\n",file=fileout, append=TRUE)

##nodes (for polygons)
nodeinfo2<-paste(" <node id=\"",Polygons$nod,"\" visible=\"true\" lat=\"",Polygons$yc,"\" lon=\"",Polygons$xc,"\"/>", sep="")
cat(nodeinfo2, sep="\n",file=fileout, append=TRUE)


##branches;
IDendnew<-IDend+nrow(edges)
IDS.ways<-(IDend+1):IDendnew
wayinfo<-paste(" <way id=\"",IDS.ways,"\" visible=\"true\">", sep="")
nd1<-paste("  <nd ref=\"", edges[,1],"\"/>",sep="")
nd2<-paste("  <nd ref=\"", edges[,2],"\"/>",sep="")
tag1<-paste("  <tag k=\"branch\" v=\"",COORD$dist2root.norm[edges[,1]],"\"/>",sep="")
tag2<-paste("  <tag k=\"name\" v=\"",paste("&#8592;  ",COORD$names[edges[,1]]," - ",COORD$names[edges[,2]], "  &#8594;", sep="     "),"\"/>",sep="")
tag3<-paste("  <tag k=\"zoomview\" v=\"",COORD$zoomview[edges[,1]],"\"/>",sep="")
last2<-rep(" </way>",nrow(edges))
branchblock<-paste(wayinfo, nd1,nd2,tag1,tag2,tag3,last2,sep="\n")
cat(branchblock, sep="\n",file=fileout, append=TRUE)

##polygons
wayinfoclade<-paste(" <way id=\"",unique(Polygons$ID),"\" visible=\"true\">", sep="")
nodesbypolyg<-matrix(Polygons$nod, ncol=60, byrow=TRUE)
nodesbypolyg<-cbind(nodesbypolyg, nodesbypolyg[,1])
makenodtags<-function(arr) {
    ndbl<-paste(paste("  <nd ref=\"", arr,"\"/>",sep=""), collapse="\n")
    ndbl
}
ndBlocks<-apply(nodesbypolyg,1,makenodtags)
tag1.1<-paste("  <tag k=\"clade\" v=\"",COORD$dist2root.norm[COORD$istip=="no"],"\"/>",sep="")
tag2.1<-paste("  <tag k=\"name\" v=\"",COORD$names[COORD$istip=="no"],"\"/>",sep="")
tag3.1<-paste("  <tag k=\"zoomview\" v=\"",COORD$zoomview[COORD$istip=="no"],"\"/>",sep="")    
last3<-rep(" </way>",nbn)
allcladeinfo<-paste(wayinfoclade, ndBlocks, tag1.1, tag2.1, tag3.1, last3, sep="\n")
cat(allcladeinfo, sep="\n",file=fileout, append=TRUE)

##finaly we add the pre-root segment
lastIDnod<-IDendnew+1
##the node
cat(paste(" <node id=\"",lastIDnod,"\" visible=\"true\" lat=\"",y.pre,"\" lon=\"",x.pre,"\"/>", sep=""), file=fileout, append=TRUE)
##and the branch:
lastID<-lastIDnod+1
cat(paste(" <way id=\"",lastID,"\" visible=\"true\">\n", sep=""), file=fileout, append=TRUE)
cat(paste("  <nd ref=\"", COORD$id[nbt+1],"\"/>\n",sep=""), file=fileout, append=TRUE)
cat(paste("  <nd ref=\"", lastIDnod,"\"/>\n",sep=""), file=fileout, append=TRUE)
cat(paste("  <tag k=\"branch\" v=\"0\"/>\n",sep=""), file=fileout, append=TRUE)
cat(paste("  <tag k=\"name\" v=\"\"/>\n",sep=""), file=fileout, append=TRUE)
cat(" </way>\n", file=fileout, append=TRUE)



##WE now give the centers of polyg with their names and zoomview.
##nodes2 (for polygon NAMES (center of polyg))
Ultim.ID<-(lastIDnod+1):(lastIDnod+1+nbn)
polyg.names<-paste(" <node id=\"",Ultim.ID,"\" visible=\"true\" lat=\"",COORD$polygon.y[COORD$istip=="no"],"\" lon=\"",COORD$polygon.x[COORD$istip=="no"],"\">", sep="")
tag1.2<-paste("  <tag k=\"cladecenter\" v=\"",COORD$dist2root.norm[COORD$istip=="no"],"\"/>",sep="")
tag2.2<-paste("  <tag k=\"name\" v=\"",COORD$names[COORD$istip=="no"],"\"/>",sep="")
tag3.2<-paste("  <tag k=\"zoomview\" v=\"",COORD$zoomview[COORD$istip=="no"],"\"/>",sep="")    
last<-rep(" </node>",nbn)
polygNameBlock<-paste(polyg.names,tag1.2,tag2.2,tag3.2,last,sep="\n")
cat(polygNameBlock, sep="\n",file=fileout, append=TRUE)



cat("</osm>\n", file=fileout, append=TRUE)


print("DONE")




##############################
##############################
###   WRITE GEOJSON FILE   ###
##############################
##############################
## we create 1 json for pointys (names)
## and another one for polygons? (later maybe)
COORD<-COORD[1:100,]
OK<-paste("{\n\"geometry\":{\n\"type\":\"Point\",\n\"coordinates\":[\n", COORD$x,",\n",COORD$y,"\n]\n},\n\"type\":\"Feature\",\n\"id\":",COORD$id,",\n\"properties\":{\n\"istip\":\"",COORD$istip,"\",\n\"zoomview\":",COORD$zoomview,",\n\"nbdesc\":",COORD$nbdesc,",\n\"names\":\"",COORD$names,"\"\n}\n}\n",sep="")
start<-paste("{\n\"type\": \"FeatureCollection\",\n\"features\": [\n", sep="")
end<-"\n]\n}\n"
out<-"tree.geojson"
cat(start, file=out)
cat(OK, sep=",\n",file=out, append=TRUE)
cat(end, file=out, append=TRUE)


