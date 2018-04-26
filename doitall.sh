make uninstall
make clean
./configure
make
sudo make install
sudo cp src/.libs/mod_tile.so /usr/lib/apache2/modules/
