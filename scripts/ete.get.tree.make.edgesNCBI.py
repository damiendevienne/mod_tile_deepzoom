#!/usr/bin/python

# This code is more complete than the previous one.
# We get the tree from the taxid given as argument and get all information from NCBI from here.
# Maybe also we could compute the coordinates their?

import sys
import numpy as np
from ete2 import Tree
from ete2 import NCBITaxa
import matplotlib.pyplot as plt

ncbi = NCBITaxa()



##t = ncbi.get_topology(ncbi.get_descendant_taxa(9443)) ##this is where we will change the number
##t.x = 5.0;
##t.y = 8.660254;
##t.alpha = 60.0;
##t.ray = 20.0;
##start = 10777;#this represents the smallest id value. Useful when combinaing multiple "flowers"

groupnb = sys.argv[1]; ##will be written

print sys.argv[1];
start = int(float(sys.argv[2]));
print "Downloading tree..."
if (sys.argv[1]=="1"):
##    t = ncbi.get_topology(ncbi.get_descendant_taxa(9443))
    t = ncbi.get_topology(ncbi.get_descendant_taxa(2157)) ##Archaea
    print "Archaeal tree downloaded from local database"
    t.x = 6.0;
    t.y = 9.660254-10.0;
    t.alpha = 30.0;
    t.ray = 10.0;
    start = start;
    osm = open("tree.osm", "w")
    
if (sys.argv[1]=="2"):
##    t = ncbi.get_topology(ncbi.get_descendant_taxa(9443))
    t = ncbi.get_topology(ncbi.get_descendant_taxa(2759)) ##Euka
    print "Eukaryotic tree downloaded from local database"
    t.x = -6.0;
    t.y = 9.660254-10.0;
    t.alpha = 150.0;
    t.ray = 10.0;
    start = start;
    osm = open("tree.osm", "a") ##a is for appending

if (sys.argv[1]=="3"): 
    t = ncbi.get_topology(ncbi.get_descendant_taxa(2)) ##Bacteria
##    t = ncbi.get_topology(ncbi.get_descendant_taxa(9443))
    print "Bacterial tree downloaded from local database"
    t.x = 0.0;
    t.y = -11.0;
    t.alpha = 270.0;
    t.ray = 10.0;
    start = start;
    osm = open("tree.osm", "a") ##a is for appending

    

t.zoomview = np.ceil(np.log2(30/t.ray));

#specis and node ids
nbsp = len(t)
spid = start
ndid = start + nbsp
rootnb = ndid+1
maxZoomView=0


#functions
def rad(deg):
    return((deg*np.pi)/180);
#ellipses
def halfCircle(x,y,r,start,end,nsteps):
    rs = np.linspace(start,end,num=nsteps)
    xc = x+r*np.cos(rs)
    yc = y+r*np.sin(rs)
    return(xc,yc)

def ellipse(x,y,r, alpha, nsteps):
    start=0
    end=np.pi+start
    rs = np.linspace(start,end,num=nsteps)
    a = r
    b = float(r)/4
    xs = a*np.cos(rs)
    ys = b*np.sin(rs)
    ##rotation 
    xs2 = x+(xs*np.cos(alpha)-ys*np.sin(alpha))
    ys2 = y+(xs*np.sin(alpha)+ys*np.cos(alpha))
    return(xs2,ys2)

def HalfCircPlusEllips(x,y,r,alpha, start, end,nsteps):
        circ = halfCircle(x,y,r,start,end, nsteps)
        elli = ellipse(x,y,r,alpha,nsteps)
        return (np.concatenate((circ[0], elli[0])),np.concatenate((circ[1], elli[1])))

#test compute xy and ray for first flower


#write osm elements
#header
##this is only written for first case:
if (sys.argv[1]=="1"):
    osm.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<osm version=\"0.6\" generator=\"HomeMade python Generator using ete toolkit\" copyright=\"OpenStreetMap and D.M. de Vienne\" attribution=\"http://www.openstreetmap.org/copyright\" license=\"http://opendatacommons.org/licenses/odbl/1-0/\">\n<bounds minlat=\"-90\" minlon=\"-180\" maxlat=\"90\" maxlon=\"180\"/>\n")
    ##we create luca's node
    osm.write("<node id=\"1000000000\" visible=\"true\" lat=\"-4.226497\" lon=\"0.000000\" >\n")
    osm.write(" <tag k=\"sci_name\" v=\"LUCA\"/>\n")
    osm.write(" <tag k=\"common_name\" v=\"Last Universal Common Ancestor\"/>\n")
    osm.write(" <tag k=\"rank\" v=\"root\"/>\n")
    osm.write(" <tag k=\"nbdesc\" v=\"3000000\"/>\n")
    osm.write(" <tag k=\"tip\" v=\"no\"/>\n")
    osm.write(" <tag k=\"zoomview\" v=\"1\"/>\n")
    osm.write("</node>\n")
    
def writeosmNode(node):
    ##we write INFO FOR EACH NODE. Clades will be delt with later on. We put less info than for the json file
    osm.write("<node id=\"%d\" visible=\"true\" lat=\"%f\" lon=\"%f\">\n" % (node.id,node.y,node.x))
    osm.write(" <tag k=\"taxid\" v=\"%d\"/>\n" % (node.taxid))
    sci_name = node.sci_name
    sci_name = sci_name.replace("<","&lt;")
    sci_name = sci_name.replace(">","&gt;")
    common_name = node.common_name
    common_name = common_name.replace("<","&lt;")
    common_name = common_name.replace(">","&gt;")
    osm.write(" <tag k=\"sci_name\" v=\"%s\"/>\n" % (sci_name))
    osm.write(" <tag k=\"common_name\" v=\"%s\"/>\n" % (common_name))
    osm.write(" <tag k=\"rank\" v=\"%s\"/>\n" % (node.rank))
    osm.write(" <tag k=\"nbdesc\" v=\"%d\"/>\n" % (node.nbdesc))
    if node.is_leaf():
        osm.write(" <tag k=\"tip\" v=\"yes\"/>\n")
    else:
        osm.write(" <tag k=\"tip\" v=\"no\"/>\n")
    osm.write(" <tag k=\"zoomview\" v=\"%s\"/>\n" % (node.zoomview))
    osm.write("</node>\n")

    
def writeosmWays(node, id):
    osm.write("<way id=\"%d\" visible=\"true\">\n" % (id))
    osm.write(" <nd ref=\"%d\"/>\n" % (node.id))
    osm.write(" <nd ref=\"%d\"/>\n" % (node.up.id))
    osm.write(" <tag k=\"branch\" v=\"1\" />\n")
    osm.write(" <tag k=\"zoomview\" v=\"%d\" />\n" % (node.up.zoomview))
    osm.write(" <tag k=\"ref\" v=\"%s\" />\n" % (groupnb))
    #we add here a name to the branches
    Upsci_name = node.up.sci_name
    Upsci_name = Upsci_name.replace("<","&lt;")
    Upsci_name = Upsci_name.replace(">","&gt;")
    Upcommon_name = node.up.common_name
    Upcommon_name = Upcommon_name.replace("<","&lt;")
    Upcommon_name = Upcommon_name.replace(">","&gt;")
    Downsci_name = node.sci_name
    Downsci_name = Downsci_name.replace("<","&lt;")
    Downsci_name = Downsci_name.replace(">","&gt;")
    Downcommon_name = node.common_name
    Downcommon_name = Downcommon_name.replace("<","&lt;")
    Downcommon_name = Downcommon_name.replace(">","&gt;")
    left = Upsci_name +  " " + Upcommon_name;
    right = Downsci_name + " " + Downcommon_name;
    if (node.x >= node.up.x): #we are on the right 
        wayName = "&#8636;  " + left + "     -     " + right + "  &#8640;"
    else: #we are on the left
        wayName = "&#8637;  " + right + "     -     " + left + "  &#8641;"
    osm.write(" <tag k=\"name\" v=\"%s\" />\n" % (wayName))
    osm.write("</way>\n")

    
def writeosmpolyg(node, ids):    
    polyg = HalfCircPlusEllips(node.x,node.y,node.ray,rad(node.alpha) + np.pi/2, rad(node.alpha) - np.pi/2, rad(node.alpha) + np.pi/2, 30)
    polygcenter = (np.mean(polyg[0]),np.mean(polyg[1]));
    #we first write the nodes
    for i in range(0,60):       
        osm.write("<node id=\"%d\" visible=\"true\" lat=\"%f\" lon=\"%f\"/>\n" % (ids[i],polyg[1][i],polyg[0][i]))
    osm.write("<way id=\"%d\" visible=\"true\">\n" % (ids[60]))
    for i in range(0,60):
        osm.write(" <nd ref=\"%d\"/>\n" % (ids[i]))
    osm.write(" <tag k=\"ref\" v=\"%s\"/>\n" % (groupnb))
    osm.write(" <tag k=\"clade\" v=\"1\" />\n")
    osm.write(" <tag k=\"taxid\" v=\"%d\"/>\n" % (node.taxid))
    osm.write(" <tag k=\"sci_name\" v=\"%s\"/>\n" % (node.sci_name))
    osm.write(" <tag k=\"common_name\" v=\"%s\"/>\n" % (node.common_name))
    osm.write(" <tag k=\"rank\" v=\"%s\"/>\n" % (node.rank))
    osm.write(" <tag k=\"nbdesc\" v=\"%d\"/>\n" % (node.nbdesc))
    osm.write(" <tag k=\"zoomview\" v=\"%d\" />\n" % (node.zoomview))
    osm.write("</way>\n")
    #and we now write the clade center.
    osm.write("<node id=\"%d\" visible=\"true\" lat=\"%f\" lon=\"%f\">\n" % (ids[61],polygcenter[1],polygcenter[0]))
    osm.write(" <tag k=\"cladecenter\" v=\"1\" />\n")
    osm.write(" <tag k=\"taxid\" v=\"%d\"/>\n" % (node.taxid))
    osm.write(" <tag k=\"sci_name\" v=\"%s\"/>\n" % (node.sci_name))
    osm.write(" <tag k=\"common_name\" v=\"%s\"/>\n" % (node.common_name))
    osm.write(" <tag k=\"rank\" v=\"%s\"/>\n" % (node.rank))
    osm.write(" <tag k=\"nbdesc\" v=\"%d\"/>\n" % (node.nbdesc))
    osm.write(" <tag k=\"zoomview\" v=\"%d\" />\n" % (node.zoomview))
    osm.write("</node>\n")
    #we add a way on which we will write the rank
    osm.write("<way id=\"%d\" visible=\"true\">\n" % (ids[62]))
    for i in range(35,45):
        osm.write(" <nd ref=\"%d\"/>\n" % (ids[i]))
    osm.write(" <tag k=\"rankname\" v=\"1\" />\n")
    osm.write(" <tag k=\"sci_name\" v=\"%s\"/>\n" % (node.sci_name))
    osm.write(" <tag k=\"zoomview\" v=\"%d\" />\n" % (node.zoomview))
    osm.write(" <tag k=\"nbdesc\" v=\"%d\"/>\n" % (node.nbdesc))
    osm.write(" <tag k=\"rank\" v=\"%s\"/>\n" % (node.rank))
    osm.write("</way>\n")

#write json elements
print "Tree traversal..."
for n in t.traverse():
    tot = 0.0;
    if n.is_leaf():
        spid = spid +1
        n.id = spid
    else:
        ndid = ndid+1
        n.id = ndid
    child = n.children;    
    for i in child: #new
        tot = tot + np.sqrt(len(i)); #new
    nbdesc = len(n);
    #add parenthesis to the common name
    if n.common_name!='':
        n.common_name = "(" + n.common_name + ")"
    n.nbdesc = nbdesc;
    nbsons = len(child);
    angles = [];
    ray = n.ray;
    for i in child:
        #i.ang = 180*(len(i)/float(nbdesc))/2;
        i.ang = 180*(np.sqrt(len(i))/tot)/2; #using sqrt we decrease difference between large and small groups
        angles.append(i.ang);
        i.ray = (ray*np.sin(rad(i.ang))/np.cos(rad(i.ang)))/(1+(np.sin(rad(i.ang))/np.cos(rad(i.ang))));
        i.dist = ray - i.ray;
    ang = np.repeat(angles, 2);
    ang = np.cumsum(ang);
    ang = ang[0::2];
    ang = [i-(90-n.alpha) for i in ang];
    cpt = 0
    for i in child:
        i.alpha = ang[cpt];
        i.x = n.x + i.dist*np.cos(rad(i.alpha));
        i.y = n.y + i.dist*np.sin(rad(i.alpha));
        i.zoomview = np.ceil(np.log2(30/i.ray))
        if i.zoomview <= 0:
            i.zoomview = 0
        if maxZoomView<i.zoomview:
            maxZoomView = i.zoomview
        cpt = cpt+1;
    #we write node info
    writeosmNode(n)

#now we redo a loop and write everything to external files;
print "Writing to external file..."
for n in t.traverse():
    if n.is_root()==False:
        ndid = ndid+1
        writeosmWays(n, ndid)    
    if n.is_leaf()==False:    
        indexes = np.linspace(ndid + 1,ndid+63,num=63)
        writeosmpolyg(n, indexes)
        ndid = ndid+63

    

        
##we add the line to LUCA (node id=1000000000)
ndid=ndid+1
osm.write("<way id=\"%d\" visible=\"true\">\n" % ndid)
osm.write(" <nd ref=\"%d\"/>\n" % rootnb)
osm.write(" <nd ref=\"1000000000\"/>\n")
osm.write(" <tag k=\"branch\" v=\"1\" />\n")
osm.write(" <tag k=\"zoomview\" v=\"0\" />\n")
osm.write(" <tag k=\"ref\" v=\"%s\" />\n" % groupnb)
osm.write("</way>\n")

##osm tail (only for last)
if (sys.argv[1]=="3"):
    osm.write("</osm>\n")

osm.close()        
print "DONE!"
print ndid;
print spid;
print ("Max zoom view : %d" % (maxZoomView)); 
