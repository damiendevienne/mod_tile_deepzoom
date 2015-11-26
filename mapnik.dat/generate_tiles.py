#!/usr/bin/env python
from math import pi,cos,sin,log,exp,atan
from subprocess import call
import sys, os
import numpy as np
from Queue import Queue

import threading

import mapnik
custom_fonts_dir = '/usr/share/fonts'
mapnik.register_fonts(custom_fonts_dir)



DEG_TO_RAD = pi/180
RAD_TO_DEG = 180/pi

# Default number of rendering threads to spawn, should be roughly equal to number of CPU cores available
NUM_THREADS = 3

##we create a dictionnary with tuples that will contain the empty tiles, so that their
##descendants are not computed at all.
EmptyTiles={};
EmptyTiles[(-1)]={}
EmptyTiles[(-1)][(0)]={}
EmptyTiles[(-1)][(0)][(0)]=0
EmptyTiles[(-2)]={}
EmptyTiles[(-2)][(0)]={}
EmptyTiles[(-2)][(0)][(0)]=0

def minmax (a,b,c):
    a = max(a,b)
    a = min(a,c)
    return a

class GoogleProjection:
    def __init__(self,levels=18):
        self.Bc = []
        self.Cc = []
        self.zc = []
        self.Ac = []
        c = 256
        for d in range(0,levels):
            e = c/2;
            self.Bc.append(c/360.0)
            self.Cc.append(c/(2 * pi))
            self.zc.append((e,e))
            self.Ac.append(c)
            c *= 2
                
    def fromLLtoPixel(self,ll,zoom):
         d = self.zc[zoom]
         e = round(d[0] + ll[0] * self.Bc[zoom])
         f = minmax(sin(DEG_TO_RAD * ll[1]),-0.9999,0.9999)
         g = round(d[1] + 0.5*log((1+f)/(1-f))*-self.Cc[zoom])
         return (e,g)
     
    def fromPixelToLL(self,px,zoom):
         e = self.zc[zoom]
         f = (px[0] - e[0])/self.Bc[zoom]
         g = (px[1] - e[1])/-self.Cc[zoom]
         h = RAD_TO_DEG * ( 2 * atan(exp(g)) - 0.5 * pi)
         return (f,h)



class RenderThread:
    def __init__(self, tile_dir, mapfile, q, printLock, maxZoom):
        self.tile_dir = tile_dir
        self.q = q
        self.m = mapnik.Map(256, 256)
        self.printLock = printLock
        # Load style XML
        mapnik.load_map(self.m, mapfile, True)
        # Obtain <Map> projection
        self.prj = mapnik.Projection(self.m.srs)
        # Projects between tile pixel co-ordinates and LatLong (EPSG:4326)
        self.tileproj = GoogleProjection(maxZoom+1)


    def render_tile(self, tile_uri, x, y, z):

        # Calculate pixel positions of bottom-left & top-right
        p0 = (x * 256, (y + 1) * 256)
        p1 = ((x + 1) * 256, y * 256)
        
        # Convert to LatLong (EPSG:4326)
        l0 = self.tileproj.fromPixelToLL(p0, z);
        l1 = self.tileproj.fromPixelToLL(p1, z);
        
        # Convert to map projection (e.g. mercator co-ords EPSG:900913)
        c0 = self.prj.forward(mapnik.Coord(l0[0],l0[1]))
        c1 = self.prj.forward(mapnik.Coord(l1[0],l1[1]))
        # Bounding box for the tile
        if hasattr(mapnik,'mapnik_version') and mapnik.mapnik_version() >= 800:
            bbox = mapnik.Box2d(c0.x,c0.y, c1.x,c1.y)
        else:
            bbox = mapnik.Envelope(c0.x,c0.y, c1.x,c1.y)
        render_size = 256
        self.m.resize(render_size, render_size)
        self.m.zoom_to_box(bbox)
        if(self.m.buffer_size < 128):
            self.m.buffer_size = 128

        # Render image with default Agg renderer
        im = mapnik.Image(render_size, render_size)
        mapnik.render(self.m, im)
        im.save(tile_uri, 'png256')


    def loop(self):
        while True:
            #Fetch a tile from the queue and render it
            r = self.q.get()
            if (r == None):
                self.q.task_done()
                break
            else:
                (name, tile_uri, x, y, z) = r

            exists= ""
            if os.path.isfile(tile_uri):
                exists= "exists"
            else:
                self.render_tile(tile_uri, x, y, z)
            bytes=os.stat(tile_uri)[6]            
            empty= ''
            if bytes == 103:
                empty = " Empty Tile "
                ##and we put this key in the dictionary
                EmptyTiles[(z)][(x)][(y)]=1;
                ##and we suppress the tile.
                os.remove(tile_uri)
            self.printLock.acquire()
            print name, ":", z, x, y, exists, empty
            self.printLock.release()
            self.q.task_done()



def render_tiles(bbox, mapfile, tile_dir, minZoom=1,maxZoom=34, name="unknown", num_threads=NUM_THREADS, tms_scheme=False):
    print "render_tiles(",bbox, mapfile, tile_dir, minZoom,maxZoom, name,")"

    # Launch rendering threads
    queue = Queue(32)
    printLock = threading.Lock()
    renderers = {}
    for i in range(num_threads):
        renderer = RenderThread(tile_dir, mapfile, queue, printLock, maxZoom)
        render_thread = threading.Thread(target=renderer.loop)
        render_thread.start()
        #print "Started render thread %s" % render_thread.getName()
        renderers[i] = render_thread

    if not os.path.isdir(tile_dir):
         os.mkdir(tile_dir)

    gprj = GoogleProjection(maxZoom+1) 

    ll0 = (bbox[0],bbox[3])
    ll1 = (bbox[2],bbox[1])

    for z in range(minZoom,maxZoom + 1):
        px0 = gprj.fromLLtoPixel(ll0,z)
        px1 = gprj.fromLLtoPixel(ll1,z)
        EmptyTiles[(z)]={}
        # check if we have directories in place
        zoom = "%s" % z
        if not os.path.isdir(tile_dir + zoom):
            os.mkdir(tile_dir + zoom)
        for x in range(int(px0[0]/256.0),int(px1[0]/256.0)+1):
            # Validate x co-ordinate
            if (x < 0) or (x >= 2**z):
                continue
            EmptyTiles[(z)][(x)]={}
            # check if we have directories in place
            str_x = "%s" % x
            if not os.path.isdir(tile_dir + zoom + '/' + str_x):
                os.mkdir(tile_dir + zoom + '/' + str_x)
            for y in range(int(px0[1]/256.0),int(px1[1]/256.0)+1):
                # Validate x co-ordinate
                if (y < 0) or (y >= 2**z):
                    continue
                EmptyTiles[(z)][(x)][(y)]=0 ##initialize the EmptyTiles Triple Dict to 0 for this combination. Will be 1 if combination is empty.
                # flip y to match OSGEO TMS spec
                if tms_scheme:
                    str_y = "%s" % ((2**z-1) - y)
                else:
                    str_y = "%s" % y
                tile_uri = tile_dir + zoom + '/' + str_x + '/' + str_y + '.png'
                ##HERE WE COULD CHECK WHTHER THE DAD TILE OF THE TILE WAS EMPTY.
                ##IF YES? NO NEED TO GENERATE THIS ONE. AN WE NEED TO ADD IT TO
                ##THE LIST OF EMPTY TILES FOR NEXT ZOOM LEVEL
                # Calculate pixel positions of bottom-left & top-right for this xy combination
                p0 = (x * 256, (y + 1) * 256)
                p1 = ((x + 1) * 256, y * 256)
                # Convert to LatLong (EPSG:4326)
                l0 = gprj.fromPixelToLL(p0, z);
                l1 = gprj.fromPixelToLL(p1, z);
                #FROM THIS WE TAKE THE CENTER OF THE TILE IN LAT LONG
                #AND WE CAN THUS COMPUTE X AND Y TILE NAMES OF THE DAD TILE!
                centerlon=(l0[0]+l1[0])/2
                centerlat=(l0[1]+l1[1])/2
                zdad = z-2 #dad zoom level
                #It may be safert o consider the grand dad and not dad because of threading that may be somehow late...                
                n = 2 ** zdad
                xdad = n * ((centerlon + 180) / 360)
                ydad = n * (1 - (np.log(np.tan(DEG_TO_RAD * centerlat) + 1/(np.cos(DEG_TO_RAD * centerlat))) / np.pi)) / 2
                xdad=int(xdad)
                ydad=int(ydad)
                ##print "TILES: %d %d  ---  %d %d " % (x,y,xdad,ydad)
                ##NEW WE CHECK IN THE Dictionary of empty tile if this dad has a value of 1
                ##If yes, no need to generate its descendants

                if (EmptyTiles[(zdad)][(xdad)][(ydad)]==1): ##The dad was empty
                    ##no need to compute the sonand we add this descendant to the list of empty tiles
                    EmptyTiles[(z)][(x)][(y)]=1 ##no need to go further
                    print "Not rendered: %d %d %d " % (z,x,y)
                else:
                    # Submit tile to be rendered into the queue
                    t = (name, tile_uri, x, y, z)
                    try:
                        queue.put(t)
                    except KeyboardInterrupt:
                        raise SystemExit("Ctrl-c detected, exiting...")
                    
    # Signal render threads to exit by sending empty request to queue
    for i in range(num_threads):
        queue.put(None)
    # wait for pending rendering jobs to complete
    queue.join()
    for i in range(num_threads):
        renderers[i].join()



if __name__ == "__main__":
    home = os.environ['HOME']
    try:
        mapfile = os.environ['MAPNIK_MAP_FILE']
    except KeyError:
        mapfile = home + "/Documents/WORK/VISUAL.TREE.OF.LIFE/mapnik.dat/osm.xml"
    try:
        tile_dir = os.environ['MAPNIK_TILE_DIR']
    except KeyError:
        tile_dir = home + "/Documents/WORK/VISUAL.TREE.OF.LIFE/mapnik.dat/tiles/"

    if not tile_dir.endswith('/'):
        tile_dir = tile_dir + '/'

    #-------------------------------------------------------------------------
    #
    # Change the following for different bounding boxes and zoom levels
    #
    # Start with an overview
    # World
    bbox = (-3,-3, 3,3)
    
    render_tiles(bbox, mapfile, tile_dir, 6,8, "Test")
