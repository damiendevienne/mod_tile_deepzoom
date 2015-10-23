#Scripts for the conversion of tree data into geographic-like ones

##ete.make.edgesNCBI.py
This script takes as input a Newick tree (obtained from ete ncbirequest in my case) and returns a collection of files that describe the tree and are used in the demi-cercles+json.R function

##demi-cercles+json.R
This script takes as input the files created by ete.make.edgesNCBI.py and returns (i) an xml file called tree.osm that will be used with osm2pgsql to feed the PostGreSQL database and (ii) a json file called tree.geojson that will be used by leaflet to display information on trees while browsing. Both files contain roughly the same information.


