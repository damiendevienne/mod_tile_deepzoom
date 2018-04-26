mod_tile_deepzoom
========
A program to efficiently render and serve map tiles. 
This is an old fork of a 2016 version that allows reaching zoom level 45 (or more?) on the maps. 

It was modified for and is use in Lifemap (lifemap.univ-lyon1.fr/explore.html) to explore the Tree of Life (42 zoom levels needed)

Compiling
=========

./autogen.sh
./configure
make
sudo make install
sudo make install-mod_tile


