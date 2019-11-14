Reading and assessing LIDAR data
================
P.Norman
18th July 2019

This mark down script is for the reading and analysis of LIDAR data, in
particular the assessment of large trees.

Lidar data can be sourced from a range of sources including some
government databases, such as <https://earthexplorer.usgs.gov/> ,
<http://elevation.fsdf.org.au/> or <http://opentopo.sdsc.edu/lidar>

For todays and tomorrows labs we will be going through the use of the
Australian governments commissioned LIDAR data sets, via
<http://elevation.fsdf.org.au/>

We will be primarily using the lidR package. To check it out, this is
the Github page <https://github.com/Jean-Romain/lidR> It is a very easy
package to use and allows for some fun graphics
:)

``` r
install.packages(c('lidR','raster','rgdal','graphics','sp','rgeos','mapview'))
```

``` r
library(lidR)
library(raster)
library(rgdal)
library(graphics)
library(sp)
library(rgeos)
library(mapview)
```

``` r
setwd("your_work_directory_path")
```

After getting your LIDAR data, you then have to extract all the files
that you want to assess into a single folder. Making sure that your end
product is either a .las or .laz file/s. First you need to get the
directory to your folder that contains the point cloud
data.

``` r
ctg <- catalog("your_path_to_the_folder_with_newly_downloaded_las/laz_files")
```

Once you have this sorted, it is time to load it into R and get a single
las file. If your data contains errors such as incorrect z values, you
can filter them out when reading catalog.

``` r
#reading in the las or laz file
lidar <- readLAS(ctg, filter = "-keep_first -drop_z_below 0 -drop_z_above 1000")
```

    ## Warning: Invalid data: 1 points with a number of returns equal to 0 found.

    ## Warning: Invalid data: 1 points with a 'return number' greater than the
    ## 'number of returns'.

If you return an error after reading don’t worry as the filter takes
care of this before you analyse the data.

Once you have read in the file, it is time to have a quick look at the
extent of the las/laz files on a map. This uses a neat package called
mapview, that allows you to plot things directly onto basemaps or
satellite images.I’m yet to find how to project this map to something
other than WGS84
    though.

``` r
plot(ctg, map = TRUE)
```

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

![Alt text](D:/GIS_data/Lidar_data/mapview_plot.jpg) With this map you
can also change the basemap to a satellite, if you click the layer stack
object in the top left corner.

I you want, you can then plot the point cloud. This may take a while and
we will trim this down to plot it in later steps.

``` r
plot(lidar,colorPalette=terrain.colors(100))
```

It should look something like this.

![Alt text](D:/GIS_data/Lidar_data/Large_scale_point_plot.png)

After inspecting you point cloud, you’ll notice that it is centred
around a single point and doesn’t necessarily get you in close to areas
that you want to inspect.

So in order to get into a group of tree or a point of interest, it is
easy to then create a subset area within the LAS file.

But before we can do this we want to know the location of the point we
are interested in. To do this we can create a polygon around the extent
of the lidar object.

``` r
#creating a polygon around the extent of the las file
ex <- extent(lidar)
las_extent <- as(ex,'SpatialPolygons')
proj_ex <- '+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
data = data.frame(f=99.9)
las_ex <- SpatialPolygonsDataFrame(las_extent, data = data, match.ID = TRUE)
proj4string(las_ex) <- CRS(proj_ex)
#This returns a shapefile if you then want to map the extent using QGIS, ArcGIS or other software.
writeOGR(obj=las_ex, dsn="Las_extent",layer="las_ex", driver="ESRI Shapefile", overwrite_layer = TRUE)
#plotting the shapefile of the extent
c1 <- gCentroid(las_ex)
plot(las_ex)
plot(c1, col='blue', add=TRUE)
```

![](Markdown_for_lidar_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Woohoo an empty box. So lets get some perspective.

An extent polygon shapefile should now be in the work directory. Now you
can open up the shapefile in ArcGIS and add a satellite basemap, then
using the identify tool to get the location of your point of interest.

If you can’t remember, here are the steps: -Open ArcGIS -Add data and
navigate to your work directory and there should now be a new folder
with a .shp file in there. Now put this on the map. -Once the polygon
has loaded, add a satellite basemap (similar location to adding data)
-Make sure the projection is one used by the polygon -Use the identify
tool to find your point of interest ![ArcGIS
method](D:/GIS_data/Lidar_data/ArcGIS_method_LI.jpg)

Once you have you coordinates, copy and paste them into the script
below. This is the x and y coordinate of the centre of subsection that
you want to assess. The radius is the size of your clip, generally in
metres but it goes on the coordinate system the original lidar data
used. So you are really just taking a cookie cut of the whole point
cloud.

``` r
pnt_of_interest <- lasclipCircle(lidar,493523.680,6867465.617, radius= 100)
```

Now you have a nice little clip that you can manipulate more easily
with. If you plot it now you can clearly see objects in your point
clouds as well as zoom around them. A bit like magic.

``` r
plot(pnt_of_interest)
```

It should look like a nice cookie cut of your large point cloud. For
example ![Subset cloud](D:/GIS_data/Lidar_data/Small_cloud.png)

After this you can find out some things about our area of interest.
First create a digital terrain model using the kriging method.

``` r
dtm <- grid_terrain(pnt_of_interest, algorithm = tin())
plot(dtm)
```

![](Markdown_for_lidar_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

After we have created the dtm, we can now 3D map what the terrain would
look like if you removed all the vegetation.

``` r
plot_dtm3d(dtm)
```

Next you can create a digital surface model so you can eye off the
trees/objects in you area.

``` r
dsm <- grid_canopy(pnt_of_interest, res = 1, dsmtin())
plot(dsm)
```

![](Markdown_for_lidar_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

So if you take the terrain model from you surface model, you can
generate tree heights :)

``` r
tree_heights <- dsm - dtm
cellStats(tree_heights, stat='max')#returns the tallest tree height value
```

    ## [1] 34.94356

Just watch if the terrain is steep, the max tree height may be
inaccurate due to a steep downslope.

``` r
plot(tree_heights)
```

![](Markdown_for_lidar_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

You can then create a slope
layer.

``` r
proj4string(dtm) <- CRS("+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
slope <- terrain(dtm, filename='', out=c('slope', 'aspect'), unit='degrees', 
                 neighbors=8)
plot(slope)
```

![](Markdown_for_lidar_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

It is then possible to then export all the layers that you have made as
geotiffs. Sometimes the projection doesn’t come through or if you would
like to change projections this is the code.

``` r
#reprojection the tree height layer
proj4string(tree_heights) <- CRS("+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#reprojecting the digital surface model
proj4string(dsm) <- CRS("+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#reprojecting the slope
proj4string(slope) <- CRS("+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
```

Now to export the geotiffs

``` r
#write the dsm
writeRaster(dsm, filename="dsmtin_1m.tif", format="GTiff", overwrite=TRUE)
#writing tree heights
writeRaster(tree_heights, filename="treeheights_raster.tif", format="GTiff", overwrite=TRUE)
#writing the dtm
writeRaster(dtm, filename="dtm_kriging.tif", format="GTiff", overwrite=TRUE)
#write the dtm slope
writeRaster(slope, filename = "slope.tif", format="GTiff", overwrite=TRUE)
```

I hope you have found this useful and a bit of fun. The thing I like to
do the most is to explore the canopies of the local rainforests and to
find tall trees.
