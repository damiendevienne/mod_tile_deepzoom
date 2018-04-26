mod_tile_deepzoom
========
A program to efficiently render and serve map tiles even for very high zoom levels. This is a 2016 fork of a original mod_tile.  

It was modified for, and is use in, Lifemap (lifemap.univ-lyon1.fr/explore.html) to explore the Tree of Life (42 zoom levels were needed)

modifications performed
=========
The only modification made here was to allow tiles with very large x and y coordinates to be rendered (>2*10^9). This required changing the "type" of x and y variables (giving coordinates of tiles) from int to long. This seems simple but took me some time. If it can be of any use to others, that's great!

Compiling
=========
./autogen.sh
./configure
make
sudo make install
sudo make install-mod_tile




See mod_tile github for more instructions https://github.com/openstreetmap/mod_tile/





