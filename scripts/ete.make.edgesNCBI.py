#!/usr/bin/python
#we modify a bit the code to deal with the NHX format of the tree
import sys

from ete2 import Tree
t = Tree(sys.argv[1])

# we name all nodes
nbt = len(t)
i = nbt ## will name the nodes
j = 0 ## will name the leaves

edges = open("edges.txt", "w")
species = open("species.txt", "w")
nodes = open("nodes.txt", "w")
sizes = open("sizes.txt", "w")
taxid = open("taxid.txt", "w")

for node in t.traverse(strategy="preorder"):
    if (node.is_leaf()==False):
        i = i+1
        if hasattr(node,"sci_name"):
            nodes.write("%s \n" % (node.sci_name))
            taxid.write("%s \n" % (node.taxid))
        else:
            nodes.write("NA \n")
            taxid.write("NA \n")
        node.name = i
        tips = len(node.get_leaves())
        sizes.write("%d \t" % (tips))
    else :
        j = j+1
        if hasattr(node,"sci_name"):
            species.write("%s \n" % (node.sci_name))
        else:
            species.write("%s \n" % (node.name))
        node.name = j
        
for node in t.iter_descendants():
#    print node.up.name , " \t " , node.name
    edges.write("%d \t %d \n" % (node.up.name, node.name))



edges.close()
species.close()
nodes.close()
sizes.close()
